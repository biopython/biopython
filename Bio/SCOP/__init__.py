# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# Modifications Copyright 2004/2005 James Casbon. All rights Reserved.
# Modifications Copyright 2010 Jeffrey Finkelstein. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Changes made by James Casbon:
# - New Astral class
# - SQL functionality for both Scop and Astral classes
# - All sunids are int not strings
#
# Code written by Jeffrey Chang to access SCOP over the internet, which
# was previously in Bio.WWW.SCOP, has now been merged into this module.


"""SCOP: Structural Classification of Proteins.

The SCOP database aims to provide a manually constructed classification of
all know protein structures into a hierarchy, the main levels of which
are family, superfamily and fold.

* "SCOP":http://scop.mrc-lmb.cam.ac.uk/legacy/
* "Introduction":http://scop.mrc-lmb.cam.ac.uk/legacy/intro.html
* "SCOP parsable files":http://scop.mrc-lmb.cam.ac.uk/legacy/parse/

The Scop object in this module represents the entire SCOP classification. It
can be built from the three SCOP parsable files, modified is so desired, and
converted back to the same file formats. A single SCOP domain (represented
by the Domain class) can be obtained from Scop using the domain's SCOP
identifier (sid).

- nodeCodeDict -- A mapping between known 2 letter node codes and a longer
                  description. The known node types are 'cl' (class), 'cf'
                  (fold), 'sf' (superfamily), 'fa' (family), 'dm' (domain),
                  'sp' (species), 'px' (domain). Additional node types may
                  be added in the future.

This module also provides code to access SCOP over the WWW.

Functions:
 - search        -- Access the main CGI script.
 - _open         -- Internally used function.

"""


import os
import re

from urllib.parse import urlencode
from urllib.request import urlopen

from . import Des
from . import Cla
from . import Hie
from . import Residues
from Bio import SeqIO
from Bio.Seq import Seq

# Turn black code style off
# fmt: off
nodeCodeDict = {"cl": "class", "cf": "fold", "sf": "superfamily",
                "fa": "family", "dm": "protein", "sp": "species", "px": "domain"}


_nodetype_to_code = {"class": "cl", "fold": "cf", "superfamily": "sf",
                     "family": "fa", "protein": "dm", "species": "sp", "domain": "px"}

nodeCodeOrder = ["ro", "cl", "cf", "sf", "fa", "dm", "sp", "px"]

astralBibIds = [10, 20, 25, 30, 35, 40, 50, 70, 90, 95, 100]

astralEvs = [10, 5, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 1e-4, 1e-5, 1e-10, 1e-15,
             1e-20, 1e-25, 1e-50]

astralEv_to_file = {10: "e+1", 5: "e+0,7", 1: "e+0", 0.5: "e-0,3", 0.1: "e-1",
                    0.05: "e-1,3", 0.01: "e-2", 0.005: "e-2,3", 0.001: "e-3",
                    1e-4: "e-4", 1e-5: "e-5", 1e-10: "e-10", 1e-15: "e-15",
                    1e-20: "e-20", 1e-25: "e-25", 1e-50: "e-50"}

astralEv_to_sql = {10: "e1", 5: "e0_7", 1: "e0", 0.5: "e_0_3", 0.1: "e_1",
                   0.05: "e_1_3", 0.01: "e_2", 0.005: "e_2_3", 0.001: "e_3",
                   1e-4: "e_4", 1e-5: "e_5", 1e-10: "e_10", 1e-15: "e_15",
                   1e-20: "e_20", 1e-25: "e_25", 1e-50: "e_50"}
# Turn black code style on
# fmt: on


def cmp_sccs(sccs1, sccs2):
    """Order SCOP concise classification strings (sccs).

    a.4.5.1 < a.4.5.11 < b.1.1.1

    A sccs (e.g. a.4.5.11) compactly represents a domain's classification.
    The letter represents the class, and the numbers are the fold,
    superfamily, and family, respectively.
    """
    s1 = sccs1.split(".")
    s2 = sccs2.split(".")

    c1, c2 = s1[0], s2[0]
    if c1 < c2:
        return -1
    if c1 > c2:
        return +1

    for c1, c2 in zip(s1[1:], s2[1:]):
        i1 = int(c1)
        i2 = int(c2)
        if i1 < i2:
            return -1
        if i1 > i2:
            return +1

    n1 = len(s1)
    n2 = len(s2)
    if n1 < n2:
        return -1
    if n1 > n2:
        return +1

    return 0


_domain_re = re.compile(r">?([\w_\.]*)\s+([\w\.]*)\s+\(([^)]*)\) (.*)")


def parse_domain(term):
    """Convert an ASTRAL header string into a Scop domain.

    An ASTRAL (http://astral.stanford.edu/) header contains a concise
    description of a SCOP domain. A very similar format is used when a
    Domain object is converted into a string.  The Domain returned by this
    method contains most of the SCOP information, but it will not be located
    within the SCOP hierarchy (i.e. The parent node will be None). The
    description is composed of the SCOP protein and species descriptions.

    A typical ASTRAL header looks like --
    >d1tpt_1 a.46.2.1 (1-70) Thymidine phosphorylase {Escherichia coli}
    """
    m = _domain_re.match(term)
    if not m:
        raise ValueError("Domain: " + term)

    dom = Domain()
    dom.sid = m.group(1)
    dom.sccs = m.group(2)
    dom.residues = Residues.Residues(m.group(3))
    if not dom.residues.pdbid:
        dom.residues.pdbid = dom.sid[1:5]
    dom.description = m.group(4).strip()

    return dom


def _open_scop_file(scop_dir_path, version, filetype):
    filename = "dir.%s.scop.txt_%s" % (filetype, version)
    handle = open(os.path.join(scop_dir_path, filename))
    return handle


class Scop:
    """The entire SCOP hierarchy.

    root -- The root node of the hierarchy
    """

    def __init__(
        self,
        cla_handle=None,
        des_handle=None,
        hie_handle=None,
        dir_path=None,
        db_handle=None,
        version=None,
    ):
        """Build the SCOP hierarchy from the SCOP parsable files, or a sql backend.

        If no file handles are given, then a Scop object with a single
        empty root node is returned.

        If a directory and version are given (with dir_path=.., version=...) or
        file handles for each file, the whole scop tree will be built in memory.

        If a MySQLdb database handle is given, the tree will be built as needed,
        minimising construction times.  To build the SQL database to the methods
        write_xxx_sql to create the tables.

        """
        self._sidDict = {}
        self._sunidDict = {}

        if all(
            h is None for h in [cla_handle, des_handle, hie_handle, dir_path, db_handle]
        ):
            return

        if dir_path is None and db_handle is None:
            if cla_handle is None or des_handle is None or hie_handle is None:
                raise RuntimeError("Need CLA, DES and HIE files to build SCOP")

        sunidDict = {}

        self.db_handle = db_handle
        try:
            if db_handle:
                # do nothing if we have a db handle, we'll do it all on the fly
                pass
            else:
                # open SCOP parseable files
                if dir_path:
                    if not version:
                        raise RuntimeError(
                            "Need SCOP version to find parsable files in directory"
                        )
                    if cla_handle or des_handle or hie_handle:
                        raise RuntimeError(
                            "Cannot specify SCOP directory and specific files"
                        )

                    cla_handle = _open_scop_file(dir_path, version, "cla")
                    des_handle = _open_scop_file(dir_path, version, "des")
                    hie_handle = _open_scop_file(dir_path, version, "hie")

                root = Node()
                domains = []
                root.sunid = 0
                root.type = "ro"
                sunidDict[root.sunid] = root
                self.root = root
                root.description = "SCOP Root"

                # Build the rest of the nodes using the DES file
                records = Des.parse(des_handle)
                for record in records:
                    if record.nodetype == "px":
                        n = Domain()
                        n.sid = record.name
                        domains.append(n)
                    else:
                        n = Node()
                    n.sunid = record.sunid
                    n.type = record.nodetype
                    n.sccs = record.sccs
                    n.description = record.description

                    sunidDict[n.sunid] = n

                # Glue all of the Nodes together using the HIE file
                records = Hie.parse(hie_handle)
                for record in records:
                    if record.sunid not in sunidDict:
                        print(record.sunid)

                    n = sunidDict[record.sunid]

                    if record.parent != "":  # Not root node

                        if record.parent not in sunidDict:
                            raise ValueError("Incomplete data?")

                        n.parent = sunidDict[record.parent]

                    for c in record.children:
                        if c not in sunidDict:
                            raise ValueError("Incomplete data?")
                        n.children.append(sunidDict[c])

                # Fill in the gaps with information from the CLA file
                sidDict = {}
                records = Cla.parse(cla_handle)
                for record in records:
                    n = sunidDict[record.sunid]
                    assert n.sccs == record.sccs
                    assert n.sid == record.sid
                    n.residues = record.residues
                    sidDict[n.sid] = n

                # Clean up
                self._sunidDict = sunidDict
                self._sidDict = sidDict
                self._domains = tuple(domains)

        finally:
            if dir_path:
                # If we opened the files, we close the files
                if cla_handle:
                    cla_handle.close()
                if des_handle:
                    des_handle.close()
                if hie_handle:
                    hie_handle.close()

    def getRoot(self):
        """Get root node."""
        return self.getNodeBySunid(0)

    def getDomainBySid(self, sid):
        """Return a domain from its sid."""
        if sid in self._sidDict:
            return self._sidDict[sid]
        if self.db_handle:
            self.getDomainFromSQL(sid=sid)
            if sid in self._sidDict:
                return self._sidDict[sid]
        else:
            return None

    def getNodeBySunid(self, sunid):
        """Return a node from its sunid."""
        if sunid in self._sunidDict:
            return self._sunidDict[sunid]
        if self.db_handle:
            self.getDomainFromSQL(sunid=sunid)
            if sunid in self._sunidDict:
                return self._sunidDict[sunid]
        else:
            return None

    def getDomains(self):
        """Return an ordered tuple of all SCOP Domains."""
        if self.db_handle:
            return self.getRoot().getDescendents("px")
        else:
            return self._domains

    def write_hie(self, handle):
        """Build an HIE SCOP parsable file from this object."""
        # We order nodes to ease comparison with original file
        for n in sorted(self._sunidDict.values(), key=lambda n: n.sunid):
            handle.write(str(n.toHieRecord()))

    def write_des(self, handle):
        """Build a DES SCOP parsable file from this object."""
        # Original SCOP file is not ordered?
        for n in sorted(self._sunidDict.values(), key=lambda n: n.sunid):
            if n != self.root:
                handle.write(str(n.toDesRecord()))

    def write_cla(self, handle):
        """Build a CLA SCOP parsable file from this object."""
        # We order nodes to ease comparison with original file
        for n in sorted(self._sidDict.values(), key=lambda n: n.sunid):
            handle.write(str(n.toClaRecord()))

    def getDomainFromSQL(self, sunid=None, sid=None):
        """Load a node from the SQL backend using sunid or sid."""
        if sunid is None and sid is None:
            return None

        cur = self.db_handle.cursor()

        if sid:
            cur.execute("SELECT sunid FROM cla WHERE sid=%s", sid)
            res = cur.fetchone()
            if res is None:
                return None
            sunid = res[0]

        cur.execute("SELECT * FROM des WHERE sunid=%s", sunid)
        data = cur.fetchone()

        if data is not None:
            n = None

            # determine if Node or Domain
            if data[1] != "px":
                n = Node(scop=self)

                cur.execute("SELECT child FROM hie WHERE parent=%s", sunid)
                children = []
                for c in cur.fetchall():
                    children.append(c[0])
                n.children = children
            else:
                n = Domain(scop=self)
                cur.execute(
                    "select sid, residues, pdbid from cla where sunid=%s", sunid
                )

                n.sid, n.residues, pdbid = cur.fetchone()
                n.residues = Residues.Residues(n.residues)
                n.residues.pdbid = pdbid
                self._sidDict[n.sid] = n

            n.sunid, n.type, n.sccs, n.description = data

            if data[1] != "ro":
                cur.execute("SELECT parent FROM hie WHERE child=%s", sunid)
                n.parent = cur.fetchone()[0]

            n.sunid = int(n.sunid)

            self._sunidDict[n.sunid] = n

    def getAscendentFromSQL(self, node, type):
        """Get ascendents using SQL backend."""
        if nodeCodeOrder.index(type) >= nodeCodeOrder.index(node.type):
            return None

        cur = self.db_handle.cursor()
        cur.execute(
            "SELECT " + type + " from cla WHERE " + node.type + "=%s", (node.sunid)
        )
        result = cur.fetchone()
        if result is not None:
            return self.getNodeBySunid(result[0])
        else:
            return None

    def getDescendentsFromSQL(self, node, type):
        """Get descendents of a node using the database backend.

        This avoids repeated iteration of SQL calls and is therefore much
        quicker than repeatedly calling node.getChildren().
        """
        if nodeCodeOrder.index(type) <= nodeCodeOrder.index(node.type):
            return []

        des_list = []

        # SQL cla table knows nothing about 'ro'
        if node.type == "ro":
            for c in node.getChildren():
                for d in self.getDescendentsFromSQL(c, type):
                    des_list.append(d)
            return des_list

        cur = self.db_handle.cursor()

        if type != "px":
            cur.execute(
                "SELECT DISTINCT des.sunid,des.type,des.sccs,description FROM "
                "cla,des WHERE cla." + node.type + "=%s AND cla." + type + "=des.sunid",
                (node.sunid),
            )
            data = cur.fetchall()
            for d in data:
                if int(d[0]) not in self._sunidDict:
                    n = Node(scop=self)
                    n.sunid, n.type, n.sccs, n.description = d
                    n.sunid = int(n.sunid)
                    self._sunidDict[n.sunid] = n

                    cur.execute("SELECT parent FROM hie WHERE child=%s", n.sunid)
                    n.parent = cur.fetchone()[0]

                    cur.execute("SELECT child FROM hie WHERE parent=%s", n.sunid)
                    children = []
                    for c in cur.fetchall():
                        children.append(c[0])
                    n.children = children

                des_list.append(self._sunidDict[int(d[0])])

        else:
            cur.execute(
                "SELECT cla.sunid,sid,pdbid,residues,cla.sccs,type,description,sp "
                "FROM cla,des where cla.sunid=des.sunid and cla." + node.type + "=%s",
                node.sunid,
            )

            data = cur.fetchall()
            for d in data:
                if int(d[0]) not in self._sunidDict:
                    n = Domain(scop=self)
                    (
                        n.sunid,
                        n.sid,
                        pdbid,
                        n.residues,
                        n.sccs,
                        n.type,
                        n.description,
                        n.parent,
                    ) = d[0:8]
                    n.residues = Residues.Residues(n.residues)
                    n.residues.pdbid = pdbid
                    n.sunid = int(n.sunid)
                    self._sunidDict[n.sunid] = n
                    self._sidDict[n.sid] = n

                des_list.append(self._sunidDict[int(d[0])])

        return des_list

    def write_hie_sql(self, handle):
        """Write HIE data to SQL database."""
        cur = handle.cursor()

        cur.execute("DROP TABLE IF EXISTS hie")
        cur.execute(
            "CREATE TABLE hie (parent INT, child INT, PRIMARY KEY (child), "
            "INDEX (parent) )"
        )

        for p in self._sunidDict.values():
            for c in p.children:
                cur.execute("INSERT INTO hie VALUES (%s,%s)" % (p.sunid, c.sunid))

    def write_cla_sql(self, handle):
        """Write CLA data to SQL database."""
        cur = handle.cursor()

        cur.execute("DROP TABLE IF EXISTS cla")
        cur.execute(
            "CREATE TABLE cla (sunid INT, sid CHAR(8), pdbid CHAR(4), "
            "residues VARCHAR(50), sccs CHAR(10), cl INT, cf INT, sf INT, fa INT, "
            "dm INT, sp INT, px INT, PRIMARY KEY (sunid), INDEX (SID) )"
        )

        for n in self._sidDict.values():
            c = n.toClaRecord()
            cur.execute(
                "INSERT INTO cla VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
                (
                    n.sunid,
                    n.sid,
                    c.residues.pdbid,
                    c.residues,
                    n.sccs,
                    n.getAscendent("cl").sunid,
                    n.getAscendent("cf").sunid,
                    n.getAscendent("sf").sunid,
                    n.getAscendent("fa").sunid,
                    n.getAscendent("dm").sunid,
                    n.getAscendent("sp").sunid,
                    n.sunid,
                ),
            )

    def write_des_sql(self, handle):
        """Write DES data to SQL database."""
        cur = handle.cursor()

        cur.execute("DROP TABLE IF EXISTS des")
        cur.execute(
            "CREATE TABLE des (sunid INT, type CHAR(2), sccs CHAR(10), "
            "description VARCHAR(255), PRIMARY KEY (sunid) )"
        )

        for n in self._sunidDict.values():
            cur.execute(
                "INSERT INTO des VALUES (%s,%s,%s,%s)",
                (n.sunid, n.type, n.sccs, n.description),
            )


class Node:
    """A node in the Scop hierarchy.

    Attributes:
     - sunid  -- SCOP unique identifiers. e.g. '14986'
     - parent -- The parent node
     - children -- A list of child nodes
     - sccs     -- SCOP concise classification string. e.g. 'a.1.1.2'
     - type     -- A 2 letter node type code. e.g. 'px' for domains
     - description -- Description text.

    """

    def __init__(self, scop=None):
        """Initialize a Node in the scop hierarchy.

        If a Scop instance is provided to the constructor, this will be used
        to lookup related references using the SQL methods.  If no instance
        is provided, it is assumed the whole tree exists and is connected.
        """
        self.sunid = ""
        self.parent = None
        self.children = []
        self.sccs = ""
        self.type = ""
        self.description = ""
        self.scop = scop

    def __str__(self):
        """Represent the node as a string."""
        s = []
        s.append(str(self.sunid))
        s.append(self.sccs)
        s.append(self.type)
        s.append(self.description)

        return " ".join(s)

    def toHieRecord(self):
        """Return an Hie.Record."""
        rec = Hie.Record()
        rec.sunid = str(self.sunid)
        if self.getParent():  # Not root node
            rec.parent = str(self.getParent().sunid)
        else:
            rec.parent = "-"
        for c in self.getChildren():
            rec.children.append(str(c.sunid))
        return rec

    def toDesRecord(self):
        """Return a Des.Record."""
        rec = Des.Record()
        rec.sunid = str(self.sunid)
        rec.nodetype = self.type
        rec.sccs = self.sccs
        rec.description = self.description
        return rec

    def getChildren(self):
        """Return a list of children of this Node."""
        if self.scop is None:
            return self.children
        else:
            return [self.scop.getNodeBySunid(x) for x in self.children]

    def getParent(self):
        """Return the parent of this Node."""
        if self.scop is None:
            return self.parent
        else:
            return self.scop.getNodeBySunid(self.parent)

    def getDescendents(self, node_type):
        """Return a list of all descendant nodes of the given type.

        Node type can be a two letter code or longer description,
        e.g. 'fa' or 'family'.
        """
        if node_type in _nodetype_to_code:
            node_type = _nodetype_to_code[node_type]

        nodes = [self]
        if self.scop:
            return self.scop.getDescendentsFromSQL(self, node_type)
        while nodes[0].type != node_type:
            if nodes[0].type == "px":
                return []  # Fell of the bottom of the hierarchy
            child_list = []
            for n in nodes:
                for child in n.getChildren():
                    child_list.append(child)
                nodes = child_list

        return nodes

    def getAscendent(self, node_type):
        """Return the ancenstor node of the given type, or None.

        Node type can be a two letter code or longer description,
        e.g. 'fa' or 'family'.
        """
        if node_type in _nodetype_to_code:
            node_type = _nodetype_to_code[node_type]

        if self.scop:
            return self.scop.getAscendentFromSQL(self, node_type)
        else:
            n = self
            if n.type == node_type:
                return None

            while n.type != node_type:
                if n.type == "ro":
                    return None  # Fell of the top of the hierarchy
                n = n.getParent()

            return n


class Domain(Node):
    """A SCOP domain. A leaf node in the Scop hierarchy.

    Attributes:
        - sid - The SCOP domain identifier. e.g. ``"d5hbib_"``
        - residues - A Residue object. It defines the collection of PDB
          atoms that make up this domain.

    """

    def __init__(self, scop=None):
        """Initialize a SCOP Domain object."""
        Node.__init__(self, scop=scop)
        self.sid = ""
        self.residues = None

    def __str__(self):
        """Represent the SCOP Domain as a string."""
        s = []
        s.append(self.sid)
        s.append(self.sccs)
        s.append("(" + str(self.residues) + ")")

        if not self.getParent():
            s.append(self.description)
        else:
            sp = self.getParent()
            dm = sp.getParent()
            s.append(dm.description)
            s.append("{" + sp.description + "}")

        return " ".join(s)

    def toDesRecord(self):
        """Return a Des.Record."""
        rec = Node.toDesRecord(self)
        rec.name = self.sid
        return rec

    def toClaRecord(self):
        """Return a Cla.Record."""
        rec = Cla.Record()
        rec.sid = self.sid
        rec.residues = self.residues
        rec.sccs = self.sccs
        rec.sunid = self.sunid

        n = self
        while n.sunid != 0:  # Not root node
            rec.hierarchy[n.type] = str(n.sunid)
            n = n.getParent()

        # Order does not matter in the hierarchy field. For more info, see
        # http://scop.mrc-lmb.cam.ac.uk/legacy/release-notes.html
        # rec.hierarchy.reverse()

        return rec


class Astral:
    """Representation of the ASTRAL database.

    Abstraction of the ASTRAL database, which has sequences for all the SCOP domains,
    as well as clusterings by percent id or evalue.
    """

    def __init__(
        self, dir_path=None, version=None, scop=None, astral_file=None, db_handle=None
    ):
        """Initialize the astral database.

        You must provide either a directory of SCOP files:
            - dir_path - string, the path to location of the scopseq-x.xx directory
                       (not the directory itself), and
            - version   -a version number.

        or, a FASTA file:
            - astral_file - string, a path to a fasta file (which will be loaded in memory)

        or, a MYSQL database:
            - db_handle - a database handle for a MYSQL database containing a table
              'astral' with the astral data in it.  This can be created
              using writeToSQL.

        """
        if astral_file is None and dir_path is None and db_handle is None:
            raise RuntimeError(
                "Need either file handle, or (dir_path + version), "
                "or database handle to construct Astral"
            )
        if not scop:
            raise RuntimeError("Must provide a Scop instance to construct")

        self.scop = scop
        self.db_handle = db_handle

        if not astral_file and not db_handle:
            if dir_path is None or version is None:
                raise RuntimeError("must provide dir_path and version")

            self.version = version
            self.path = os.path.join(dir_path, "scopseq-%s" % version)
            astral_file = "astral-scopdom-seqres-all-%s.fa" % self.version
            astral_file = os.path.join(self.path, astral_file)

        if astral_file:
            # Build a dictionary of SeqRecord objects in the FASTA file, IN MEMORY
            self.fasta_dict = SeqIO.to_dict(SeqIO.parse(astral_file, "fasta"))

        self.astral_file = astral_file
        self.EvDatasets = {}
        self.EvDatahash = {}
        self.IdDatasets = {}
        self.IdDatahash = {}

    def domainsClusteredByEv(self, id):
        """Get domains clustered by evalue."""
        if id not in self.EvDatasets:
            if self.db_handle:
                self.EvDatasets[id] = self.getAstralDomainsFromSQL(astralEv_to_sql[id])

            else:
                if not self.path:
                    raise RuntimeError("No scopseq directory specified")

                file_prefix = "astral-scopdom-seqres-sel-gs"
                filename = "%s-e100m-%s-%s.id" % (
                    file_prefix,
                    astralEv_to_file[id],
                    self.version,
                )
                filename = os.path.join(self.path, filename)
                self.EvDatasets[id] = self.getAstralDomainsFromFile(filename)
        return self.EvDatasets[id]

    def domainsClusteredById(self, id):
        """Get domains clustered by percentage identity."""
        if id not in self.IdDatasets:
            if self.db_handle:
                self.IdDatasets[id] = self.getAstralDomainsFromSQL("id" + str(id))
            else:
                if not self.path:
                    raise RuntimeError("No scopseq directory specified")

                file_prefix = "astral-scopdom-seqres-sel-gs"
                filename = "%s-bib-%s-%s.id" % (file_prefix, id, self.version)
                filename = os.path.join(self.path, filename)
                self.IdDatasets[id] = self.getAstralDomainsFromFile(filename)
        return self.IdDatasets[id]

    def getAstralDomainsFromFile(self, filename=None, file_handle=None):
        """Get the scop domains from a file containing a list of sids."""
        if file_handle is None and filename is None:
            raise RuntimeError("You must provide a filename or handle")
        if not file_handle:
            file_handle = open(filename)
        doms = []
        while True:
            line = file_handle.readline()
            if not line:
                break
            line = line.rstrip()
            doms.append(line)
        if filename:
            file_handle.close()

        doms = [a for a in doms if a[0] == "d"]
        doms = [self.scop.getDomainBySid(x) for x in doms]
        return doms

    def getAstralDomainsFromSQL(self, column):
        """Load ASTRAL domains from the MySQL database.

        Load a set of astral domains from a column in the astral table of a MYSQL
        database (which can be created with writeToSQL(...).
        """
        cur = self.db_handle.cursor()
        cur.execute("SELECT sid FROM astral WHERE " + column + "=1")
        data = cur.fetchall()
        data = [self.scop.getDomainBySid(x[0]) for x in data]

        return data

    def getSeqBySid(self, domain):
        """Get the seq record of a given domain from its sid."""
        if self.db_handle is None:
            return self.fasta_dict[domain].seq
        else:
            cur = self.db_handle.cursor()
            cur.execute("SELECT seq FROM astral WHERE sid=%s", domain)
            return Seq(cur.fetchone()[0])

    def getSeq(self, domain):
        """Return seq associated with domain."""
        return self.getSeqBySid(domain.sid)

    def hashedDomainsById(self, id):
        """Get domains clustered by sequence identity in a dict."""
        if id not in self.IdDatahash:
            self.IdDatahash[id] = {}
            for d in self.domainsClusteredById(id):
                self.IdDatahash[id][d] = 1
        return self.IdDatahash[id]

    def hashedDomainsByEv(self, id):
        """Get domains clustered by evalue in a dict."""
        if id not in self.EvDatahash:
            self.EvDatahash[id] = {}
            for d in self.domainsClusteredByEv(id):
                self.EvDatahash[id][d] = 1
        return self.EvDatahash[id]

    def isDomainInId(self, dom, id):
        """Return true if the domain is in the astral clusters for percent ID."""
        return dom in self.hashedDomainsById(id)

    def isDomainInEv(self, dom, id):
        """Return true if the domain is in the ASTRAL clusters for evalues."""
        return dom in self.hashedDomainsByEv(id)

    def writeToSQL(self, db_handle):
        """Write the ASTRAL database to a MYSQL database."""
        cur = db_handle.cursor()

        cur.execute("DROP TABLE IF EXISTS astral")
        cur.execute("CREATE TABLE astral (sid CHAR(8), seq TEXT, PRIMARY KEY (sid))")

        for dom in self.fasta_dict:
            cur.execute(
                "INSERT INTO astral (sid,seq) values (%s,%s)",
                (dom, self.fasta_dict[dom].seq.data),
            )

        for i in astralBibIds:
            cur.execute("ALTER TABLE astral ADD (id" + str(i) + " TINYINT)")

            for d in self.domainsClusteredById(i):
                cur.execute("UPDATE astral SET id" + str(i) + "=1  WHERE sid=%s", d.sid)

        for ev in astralEvs:
            cur.execute("ALTER TABLE astral ADD (" + astralEv_to_sql[ev] + " TINYINT)")

            for d in self.domainsClusteredByEv(ev):
                cur.execute(
                    "UPDATE astral SET " + astralEv_to_sql[ev] + "=1  WHERE sid=%s",
                    d.sid,
                )


def search(
    pdb=None,
    key=None,
    sid=None,
    disp=None,
    dir=None,
    loc=None,
    cgi="http://scop.mrc-lmb.cam.ac.uk/legacy/search.cgi",
    **keywds
):
    """Access SCOP search and return a handle to the results.

    Access search.cgi and return a handle to the results.  See the
    online help file for an explanation of the parameters:
    http://scop.mrc-lmb.cam.ac.uk/legacy/help.html

    Raises an IOError if there's a network error.

    """
    params = {"pdb": pdb, "key": key, "sid": sid, "disp": disp, "dir": dir, "loc": loc}
    variables = {}
    for k, v in params.items():
        if v is not None:
            variables[k] = v
    variables.update(keywds)
    return _open(cgi, variables)


def _open(cgi, params=None, get=1):
    """Open a handle to SCOP and return it (PRIVATE).

    Open a handle to SCOP.  cgi is the URL for the cgi script to access.
    params is a dictionary with the options to pass to it.  get is a boolean
    that describes whether a GET should be used.

    """
    # Open a handle to SCOP.
    if params is None:
        params = {}
    options = urlencode(params)
    if get:  # do a GET
        if options:
            cgi += "?" + options
        handle = urlopen(cgi)
    else:  # do a POST
        handle = urlopen(cgi, data=options)
    return handle
