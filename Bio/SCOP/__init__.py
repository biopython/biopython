# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


""" SCOP: Structural Classification of Proteins.

The SCOP database aims to provide a manually constructed classification of
all know protein structures into a hierarchy, the main levels of which
are family, superfamily and fold.

* "SCOP":http://scop.mrc-lmb.cam.ac.uk/scop/

* "Introduction":http://scop.mrc-lmb.cam.ac.uk/scop/intro.html

* "SCOP parsable files":http://scop.mrc-lmb.cam.ac.uk/scop/parse/

The Scop object in this module represents the entire SCOP classification. It
can be built from the three SCOP parsable files, modified is so desired, and
converted back to the same file formats. A single SCOP domain (represented
by the Domain class) can be obtained from Scop using the domain's SCOP
identifier (sid).


nodeCodeDict  -- A mapping between known 2 letter node codes and a longer
                  description. The known node types are 'cl' (class), 'cf'
                  (fold), 'sf' (superfamily), 'fa' (family), 'dm' (domain), 
                  'sp' (species), 'px' (domain). Additional node types may
                  be added in the future.
"""


__all__ = [
    'Dom',
    'Cla',
    'Hie',
    'Des',
    'FileIndex',
    'Residues',
    'Raf',
    'parse_domain',
    'cmp_sccs',
    'Scop'
    ]



from types import *

import Des
import Cla
import Hie
from Residues import * 


nodeCodeDict = { 'cl':'class', 'cf':'fold', 'sf':'superfamily',
                 'fa':'family', 'dm':'protein', 'sp':'species', 'px':'domain'}


def cmp_sccs(sccs1, sccs2) :
    """Order SCOP concise classification strings (sccs).

    a.4.5.1 < a.4.5.11 < b.1.1.1 

    A sccs (e.g. a.4.5.11) compactly represents a domain's classification.
    The letter represents the class, and the numbers are the fold,
    superfamily, and family, respectively.

    """

    s1 = sccs1.split(".")
    s2 = sccs2.split(".")

    if s1[0] != s2[0]: return cmp(s1[0], s2[0])

    s1 = map(int, s1[1:])
    s2 = map(int, s2[1:])

    return cmp(s1,s2)


_domain_re = re.compile(r">?([\w_\.]*)\s+([\w\.]*)\s+\(([^)]*)\) (.*)")
def parse_domain(str) :
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

    m = _domain_re.match(str)
    if (not m) : raise SyntaxError, "Domain: "+ str

    dom = Domain()
    dom.sid = m.group(1)
    dom.sccs = m.group(2)
    dom.residues = Residues(m.group(3))
    if not dom.residues.pdbid :
        dom.residues.pdbid= dom.sid[1:5]
    dom.description = m.group(4).strip()

    return dom



class Scop:
    """The entire SCOP hierarchy.

    root -- The root node of the hierarchy 
    """
    def __init__(self, cla_handle=None, des_handle=None, hie_handle=None):
        """Build the SCOP hierarchy from the SCOP parsable files.

        If no file handles are given, then a Scop object with a single
        empty root node is returned.
        """
        self.root = Node()
        self._sidDict = {}
        self._sunidDict = {}

        if cla_handle==des_handle==hie_handle==None: return 

        if cla_handle == None or des_handle==None or hie_handle==None:
             raise RuntimeError,"Need CLA, DES and HIE files to build SCOP" 

        sunidDict = {}

        #Build the root node
        root = Node()
        domains = []
        root.sunid='0'
        sunidDict[root.sunid] = root
        self.root = root
        self.description = 'SCOP'

        # Build the rest of the nodes using the DES file
        i = Des.Iterator(des_handle, Des.Parser())
        while 1 :
            rec = i.next() 
            if rec is None : break
            if rec.nodetype =='px' :
                n = Domain()
                n.sid = rec.name
                domains.append(n)
            else : 
                n = Node()
            n.sunid = rec.sunid
            n.type = rec.nodetype
            n.sccs = rec.sccs
            n.description = rec.description
            
            sunidDict[n.sunid] = n
 
        # Glue all of the Nodes together using the HIE file
        i = Hie.Iterator(hie_handle, Hie.Parser())
        while 1 :
            rec = i.next()
            if rec is None : break
            if not sunidDict.has_key(rec.sunid) :
                print rec.sunid
    
            n = sunidDict[rec.sunid]

            if rec.parent : # Not root node
                if not sunidDict.has_key(rec.parent):
                    raise SyntaxError, "Incomplete data?"
                                       
                n.parent = sunidDict[rec.parent]
                
            for c in rec.children:
                if not sunidDict.has_key(c) :
                    raise SyntaxError, "Incomplete data?"
                n.children.append(sunidDict[c])


        # Fill in the gaps with information from the CLA file
        sidDict = {}
        i = Cla.Iterator(cla_handle, Cla.Parser())
        while 1 :
            rec = i.next()
            if rec is None : break
            n = sunidDict[rec.sunid]
            assert n.sccs == rec.sccs
            assert n.sid == rec.sid
            n.residues = rec.residues
            sidDict[n.sid] = n

        # Clean up
        self._sunidDict = sunidDict
        self._sidDict = sidDict
        self._domains = tuple(domains)




    def getDomainBySid(self, sid) :
        return self._sidDict[sid]


    def getDomains(self) :
        """Returns an ordered tuple of all SCOP Domains"""
        return self._domains






    def write_hie(self, handle) :
        """Build an HIE SCOP parsable file from this object"""
        nodes = self._sunidDict.values()
        # We order nodes to ease comparison with original file
        nodes.sort(lambda n1,n2: cmp(n1.sunid, n2.sunid))

        for n in nodes :
            handle.write(str(n.toHieRecord()))


    def write_des(self, handle) :
        """Build a DES SCOP parsable file from this object""" 
        nodes = self._sunidDict.values()
        # Origional SCOP file is not ordered?
        nodes.sort(lambda n1,n2: cmp(n1.sunid, n2.sunid))

        for n in nodes :
            if n != self.root :
                handle.write(str(n.toDesRecord()))


    def write_cla(self, handle) :
        """Build a CLA SCOP parsable file from this object"""                
        nodes = self._sidDict.values()
        # We order nodes to ease comparison with original file
        nodes.sort(lambda n1,n2: cmp(n1.sunid, n2.sunid))

        for n in nodes :
            handle.write(str(n.toClaRecord()))


  
class Node :
    """ A node in the Scop hierarchy

    sunid  -- SCOP unique identifiers. e.g. '14986'

    parent -- The parent node

    children -- A list of child nodes

    sccs     -- SCOP concise classification string. e.g. 'a.1.1.2'

    type     -- A 2 letter node type code. e.g. 'px' for domains

    description -- 
        
    """
    def __init__(self) :
        self.sunid=''    
        self.parent = None
        self.children=[]
        self.sccs = ''   
        self.type =''    
        self.description =''

    def __str__(self) :
        s = []
        s.append(self.sunid)
        s.append(self.sccs)
        s.append(self.type)
        s.append(self.description)

        return " ".join(s)

    def toHieRecord(self):
        """Return an Hie.Record"""
        rec = Hie.Record()
        rec.sunid =self.sunid
        if self.parent : #Not root node
            rec.parent = self.parent.sunid
        for c in self.children :
            rec.children.append(c.sunid)
        return rec
    
    def toDesRecord(self):
        """Return a Des.Record"""        
        rec = Des.Record()
        rec.sunid = self.sunid
        rec.nodetype = self.type
        rec.sccs = self.sccs
        rec.description = self.description
        return rec


class Domain(Node) :
    """ A SCOP domain. A leaf node in the Scop hierarchy.

    sid      -- The SCOP domain identifier. e.g. 'd5hbib_'

    residues -- A Residue object. It defines the collection
                  of PDB atoms that make up this domain.
    """
    def __init__(self) :
        Node.__init__(self)
        self.sid = ''         
        self.residues = None

    def __str__(self) :
        s = []
        s.append(self.sid)
        s.append(self.sccs)
        s.append("("+str(self.residues)+")")

        if not self.parent :
            s.append(self.description)
        else :
            sp = self.parent
            dm = sp.parent
            s.append(dm.description)
            s.append("{"+sp.description+"}")

        return " ".join(s)

    def toDesRecord(self):
        """Return a Des.Record"""
        rec = Node.toDesRecord(self)
        rec.name = self.sid
        return rec

    def toClaRecord(self) :
        """Return a Cla.Record"""        
        rec = Cla.Record()
        rec.sid = self.sid
        rec.residues = self.residues
        rec.sccs = self.sccs
        rec.sunid = self.sunid
        
        n = self
        while n.sunid !='0': #Not root node
            rec.hierarchy.append( (n.type, n.sunid) )
            n = n.parent
        rec.hierarchy.reverse()
       
        return rec
            
















