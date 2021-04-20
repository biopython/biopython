# Copyright 2013 by David Arenillas and Anthony Mathelier. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
"""Provides read access to a JASPAR5 formatted database.

This modules requires MySQLdb to be installed.

Example, substitute the your database credentials as
appropriate::

        from Bio.motifs.jaspar.db import JASPAR5
        JASPAR_DB_HOST = "hostname.example.org"
        JASPAR_DB_NAME = "JASPAR2018"
        JASPAR_DB_USER = "guest"
        JASPAR_DB_PASS = "guest"

        jdb = JASPAR5(
            host=JASPAR_DB_HOST,
            name=JASPAR_DB_NAME,
            user=JASPAR_DB_USER,
            password=JASPAR_DB_PASS
        )
        ets1 = jdb.fetch_motif_by_id('MA0098')
        print(ets1)
    TF name ETS1
    Matrix ID   MA0098.3
    Collection  CORE
    TF class    Tryptophan cluster factors
    TF family   Ets-related factors
    Species 9606
    Taxonomic group vertebrates
    Accession   ['P14921']
    Data type used  HT-SELEX
    Medline 20517297
    PAZAR ID    TF0000070
    Comments    Data is from Taipale HTSELEX DBD (2013)
    Matrix:
            0      1      2      3      4      5      6      7      8      9
    A: 2683.00 180.00 425.00   0.00   0.00 2683.00 2683.00 1102.00  89.00 803.00
    C: 210.00 2683.00 2683.00  21.00   0.00   0.00   9.00  21.00 712.00 401.00
    G: 640.00 297.00   7.00 2683.00 2683.00   0.00  31.00 1580.00 124.00 1083.00
    T: 241.00  22.00   0.00   0.00  12.00   0.00 909.00  12.00 1970.00 396.00

        motifs = jdb.fetch_motifs(
            collection = 'CORE',
            tax_group = ['vertebrates', 'insects'],
            tf_class = 'Homeo domain factors',
            tf_family = ['TALE-type homeo domain factors', 'POU domain factors'],
            min_ic = 12
        )
        for motif in motifs:
            pass # do something with the motif
"""


import warnings
from Bio import BiopythonWarning
from Bio import MissingPythonDependencyError

try:
    import MySQLdb as mdb
except ImportError:
    raise MissingPythonDependencyError(
        "Install MySQLdb if you want to use Bio.motifs.jaspar.db"
    )

from Bio.motifs import jaspar, matrix


JASPAR_DFLT_COLLECTION = "CORE"


class JASPAR5:
    """Class representing a JASPAR5 database.

    Class representing a JASPAR5 DB. The methods within are loosely based
    on the perl TFBS::DB::JASPAR5 module.

    Note: We will only implement reading of JASPAR motifs from the DB.
    Unlike the perl module, we will not attempt to implement any methods to
    store JASPAR motifs or create a new DB at this time.
    """

    def __init__(self, host=None, name=None, user=None, password=None):
        """Construct a JASPAR5 instance and connect to specified DB.

        Arguments:
         - host - host name of the the JASPAR DB server
         - name - name of the JASPAR database
         - user - user name to connect to the JASPAR DB
         - password - JASPAR DB password

        """
        self.name = name
        self.host = host
        self.user = user
        self.password = password

        self.dbh = mdb.connect(host, user, password, name)

    def __str__(self):
        """Return a string represention of the JASPAR5 DB connection."""
        return r"%s\@%s:%s" % (self.user, self.host, self.name)

    def fetch_motif_by_id(self, id):
        """Fetch a single JASPAR motif from the DB by its JASPAR matrix ID.

        Example id 'MA0001.1'.

        Arguments:
         - id - JASPAR matrix ID. This may be a fully specified ID including
                the version number (e.g. MA0049.2) or just the base ID (e.g.
                MA0049). If only a base ID is provided, the latest version is
                returned.

        Returns:
         - A Bio.motifs.jaspar.Motif object

        **NOTE:** The perl TFBS module allows you to specify the type of matrix
        to return (PFM, PWM, ICM) but matrices are always stored in JASPAR as
        PFMs so this does not really belong here. Once a PFM is fetched the
        pwm() and pssm() methods can be called to return the normalized and
        log-odds matrices.

        """
        # separate stable ID and version number
        (base_id, version) = jaspar.split_jaspar_id(id)
        if not version:
            # if ID contains no version portion, fetch the latest version
            version = self._fetch_latest_version(base_id)

        # fetch internal JASPAR matrix ID - also a check for validity
        int_id = None
        if version:
            int_id = self._fetch_internal_id(base_id, version)

        # fetch JASPAR motif using internal ID
        motif = None
        if int_id:
            motif = self._fetch_motif_by_internal_id(int_id)

        return motif

    def fetch_motifs_by_name(self, name):
        """Fetch a list of JASPAR motifs from a JASPAR DB by the given TF name(s).

        Arguments:
        name - a single name or list of names
        Returns:
        A list of Bio.motifs.jaspar.Motif objects

        Notes:
        Names are not guaranteed to be unique. There may be more than one
        motif with the same name. Therefore even if name specifies a single
        name, a list of motifs is returned. This just calls
        self.fetch_motifs(collection = None, tf_name = name).

        This behaviour is different from the TFBS perl module's
        get_Matrix_by_name() method which always returns a single matrix,
        issuing a warning message and returning the first matrix retrieved
        in the case where multiple matrices have the same name.

        """
        return self.fetch_motifs(collection=None, tf_name=name)

    def fetch_motifs(
        self,
        collection=JASPAR_DFLT_COLLECTION,
        tf_name=None,
        tf_class=None,
        tf_family=None,
        matrix_id=None,
        tax_group=None,
        species=None,
        pazar_id=None,
        data_type=None,
        medline=None,
        min_ic=0,
        min_length=0,
        min_sites=0,
        all=False,
        all_versions=False,
    ):
        """Fetch jaspar.Record (list) of motifs using selection criteria.

        Arguments::

            Except where obvious, all selection criteria arguments may be
            specified as a single value or a list of values. Motifs must
            meet ALL the specified selection criteria to be returned with
            the precedent exceptions noted below.

            all         - Takes precedent of all other selection criteria.
                          Every motif is returned. If 'all_versions' is also
                          specified, all versions of every motif are returned,
                          otherwise just the latest version of every motif is
                          returned.
            matrix_id   - Takes precedence over all other selection criteria
                          except 'all'.  Only motifs with the given JASPAR
                          matrix ID(s) are returned. A matrix ID may be
                          specified as just a base ID or full JASPAR IDs
                          including version number. If only a base ID is
                          provided for specific motif(s), then just the latest
                          version of those motif(s) are returned unless
                          'all_versions' is also specified.
            collection  - Only motifs from the specified JASPAR collection(s)
                          are returned. NOTE - if not specified, the collection
                          defaults to CORE for all other selection criteria
                          except 'all' and 'matrix_id'. To apply the other
                          selection criteria across all JASPAR collections,
                          explicitly set collection=None.
            tf_name     - Only motifs with the given name(s) are returned.
            tf_class    - Only motifs of the given TF class(es) are returned.
            tf_family   - Only motifs from the given TF families are returned.
            tax_group   - Only motifs belonging to the given taxonomic
                          supergroups are returned (e.g. 'vertebrates',
                          'insects', 'nematodes' etc.)
            species     - Only motifs derived from the given species are
                          returned.  Species are specified as taxonomy IDs.
            data_type   - Only motifs generated with the given data type (e.g.
                          ('ChIP-seq', 'PBM', 'SELEX' etc.) are returned.
                          NOTE - must match exactly as stored in the database.
            pazar_id    - Only motifs with the given PAZAR TF ID are returned.
            medline     - Only motifs with the given medline (PubmMed IDs) are
                          returned.
            min_ic      - Only motifs whose profile matrices have at least this
                          information content (specificty) are returned.
            min_length  - Only motifs whose profiles are of at least this
                          length are returned.
            min_sites   - Only motifs compiled from at least these many binding
                          sites are returned.
            all_versions- Unless specified, just the latest version of motifs
                          determined by the other selection criteria are
                          returned. Otherwise all versions of the selected
                          motifs are returned.

        Returns:
            - A Bio.motifs.jaspar.Record (list) of motifs.

        """
        # Fetch the internal IDs of the motifs using the criteria provided
        int_ids = self._fetch_internal_id_list(
            collection=collection,
            tf_name=tf_name,
            tf_class=tf_class,
            tf_family=tf_family,
            matrix_id=matrix_id,
            tax_group=tax_group,
            species=species,
            pazar_id=pazar_id,
            data_type=data_type,
            medline=medline,
            all=all,
            all_versions=all_versions,
        )

        record = jaspar.Record()

        """
        Now further filter motifs returned above based on any specified
        matrix specific criteria.
        """
        for int_id in int_ids:
            motif = self._fetch_motif_by_internal_id(int_id)

            # Filter motifs to those with matrix IC greater than min_ic
            if min_ic:
                if motif.pssm.mean() < min_ic:
                    continue

            # Filter motifs to those with minimum length of min_length
            if min_length:
                if motif.length < min_length:
                    continue

            # XXX We could also supply a max_length filter.

            """
            Filter motifs to those composed of at least this many sites.
            The perl TFBS module assumes column sums may be different but
            this should be strictly enforced here we will ignore this and
            just use the first column sum.
            """
            if min_sites:
                num_sites = sum(motif.counts[nt][0] for nt in motif.alphabet)
                if num_sites < min_sites:
                    continue

            record.append(motif)

        return record

    def _fetch_latest_version(self, base_id):
        """Get the latest version number for the given base_id (PRIVATE)."""
        cur = self.dbh.cursor()
        cur.execute(
            "select VERSION from MATRIX where BASE_id = %s order by VERSION"
            " desc limit 1",
            (base_id,),
        )

        row = cur.fetchone()

        latest = None
        if row:
            latest = row[0]
        else:
            warnings.warn(
                "Failed to fetch latest version number for JASPAR motif"
                f" with base ID '{base_id}'. No JASPAR motif with this"
                " base ID appears to exist in the database.",
                BiopythonWarning,
            )

        return latest

    def _fetch_internal_id(self, base_id, version):
        """Fetch the internal id for a base id + version (PRIVATE).

        Also checks if this combo exists or not.
        """
        cur = self.dbh.cursor()
        cur.execute(
            "select id from MATRIX where BASE_id = %s and VERSION = %s",
            (base_id, version),
        )

        row = cur.fetchone()

        int_id = None
        if row:
            int_id = row[0]
        else:
            warnings.warn(
                "Failed to fetch internal database ID for JASPAR motif"
                f" with matrix ID '{base_id}.{version}'. No JASPAR motif"
                " with this matrix ID appears to exist.",
                BiopythonWarning,
            )

        return int_id

    def _fetch_motif_by_internal_id(self, int_id):
        """Fetch basic motif information (PRIVATE)."""
        cur = self.dbh.cursor()
        cur.execute(
            "select BASE_ID, VERSION, COLLECTION, NAME from MATRIX where id = %s",
            (int_id,),
        )

        row = cur.fetchone()

        # This should never happen as it is an internal method. If it does
        # we should probably raise an exception
        if not row:
            warnings.warn(
                f"Could not fetch JASPAR motif with internal ID = {int_id}",
                BiopythonWarning,
            )
            return None

        base_id = row[0]
        version = row[1]
        collection = row[2]
        name = row[3]

        matrix_id = "".join([base_id, ".", str(version)])

        # fetch the counts matrix
        counts = self._fetch_counts_matrix(int_id)

        # Create new JASPAR motif
        motif = jaspar.Motif(matrix_id, name, collection=collection, counts=counts)

        # fetch species
        cur.execute("select TAX_ID from MATRIX_SPECIES where id = %s", (int_id,))
        tax_ids = []
        rows = cur.fetchall()
        for row in rows:
            tax_ids.append(row[0])

        # Many JASPAR motifs (especially those not in the CORE collection)
        # do not have taxonomy IDs. So this warning would get annoying.
        # if not tax_ids:
        #     warnings.warn("Could not fetch any taxonomy IDs for JASPAR motif"
        #                   " {0}".format(motif.matrix_id), BiopythonWarning)

        motif.species = tax_ids

        # fetch protein accession numbers
        cur.execute("select ACC FROM MATRIX_PROTEIN where id = %s", (int_id,))
        accs = []
        rows = cur.fetchall()
        for row in rows:
            accs.append(row[0])

        # Similarly as for taxonomy IDs, it would get annoying to print
        # warnings for JASPAR motifs which do not have accession numbers.

        motif.acc = accs

        # fetch remaining annotation as tags from the ANNOTATION table
        cur.execute("select TAG, VAL from MATRIX_ANNOTATION where id = %s", (int_id,))
        rows = cur.fetchall()
        for row in rows:
            attr = row[0]
            val = row[1]
            if attr == "class":
                motif.tf_class = val
            elif attr == "family":
                motif.tf_family = val
            elif attr == "tax_group":
                motif.tax_group = val
            elif attr == "type":
                motif.data_type = val
            elif attr == "pazar_tf_id":
                motif.pazar_id = val
            elif attr == "medline":
                motif.medline = val
            elif attr == "comment":
                motif.comment = val
            else:
                # TODO If we were to implement additional abitrary tags
                # motif.tag(attr, val)
                pass

        return motif

    def _fetch_counts_matrix(self, int_id):
        """Fetch the counts matrix from the JASPAR DB by the internal ID (PRIVATE).

        Returns a Bio.motifs.matrix.GenericPositionMatrix
        """
        counts = {}
        cur = self.dbh.cursor()

        for base in "ACGT":
            base_counts = []

            cur.execute(
                "select val from MATRIX_DATA where ID = %s and row = %s order by col",
                (int_id, base),
            )

            rows = cur.fetchall()
            for row in rows:
                base_counts.append(row[0])

            counts[base] = [float(x) for x in base_counts]

        return matrix.GenericPositionMatrix("ACGT", counts)

    def _fetch_internal_id_list(
        self,
        collection=JASPAR_DFLT_COLLECTION,
        tf_name=None,
        tf_class=None,
        tf_family=None,
        matrix_id=None,
        tax_group=None,
        species=None,
        pazar_id=None,
        data_type=None,
        medline=None,
        all=False,
        all_versions=False,
    ):
        """Fetch list of internal JASPAR motif IDs.

        Fetch a list of internal JASPAR motif IDs based on various passed
        parameters which may then be used to fetch the rest of the motif data.

        Caller:
            fetch_motifs()

        Arguments:
            See arguments sections of fetch_motifs()

        Returns:
            A list of internal JASPAR motif IDs which match the given
            selection criteria arguments.


        Build an SQL query based on the selection arguments provided.

        1: First add table joins and sub-clauses for criteria corresponding to
           named fields from the MATRIX and MATRIX_SPECIES tables such as
           collection, matrix ID, name, species etc.

        2: Then add joins/sub-clauses for tag/value parameters from the
           MATRIX_ANNOTATION table.

        For the surviving matrices, the responsibility to do matrix-based
        feature filtering such as ic, number of sites etc, fall on the
        calling fetch_motifs() method.

        """
        int_ids = []

        cur = self.dbh.cursor()

        """
        Special case 1: fetch ALL motifs. Highest priority.
        Ignore all other selection arguments.
        """
        if all:
            cur.execute("select ID from MATRIX")
            rows = cur.fetchall()

            for row in rows:
                int_ids.append(row[0])

            return int_ids

        """
        Special case 2: fetch specific motifs by their JASPAR IDs. This
        has higher priority than any other except the above 'all' case.
        Ignore all other selection arguments.
        """
        if matrix_id:
            """
            These might be either stable IDs or stable_ID.version.
            If just stable ID and if all_versions == 1, return all versions,
            otherwise just the latest
            """
            if all_versions:
                for id in matrix_id:
                    # ignore vesion here, this is a stupidity filter
                    (base_id, version) = jaspar.split_jaspar_id(id)
                    cur.execute("select ID from MATRIX where BASE_ID = %s", (base_id,))

                    rows = cur.fetchall()
                    for row in rows:
                        int_ids.append(row[0])
            else:
                # only the lastest version, or the requested version
                for id in matrix_id:
                    (base_id, version) = jaspar.split_jaspar_id(id)

                    if not version:
                        version = self._fetch_latest_version(base_id)

                    int_id = None
                    if version:
                        int_id = self._fetch_internal_id(base_id, version)

                    if int_id:
                        int_ids.append(int_id)

            return int_ids

        tables = ["MATRIX m"]
        where_clauses = []

        # Select by MATRIX.COLLECTION
        if collection:
            if isinstance(collection, list):
                # Multiple collections passed in as a list
                clause = "m.COLLECTION in ('"
                clause = "".join([clause, "','".join(collection)])
                clause = "".join([clause, "')"])
            else:
                # A single collection - typical usage
                clause = "m.COLLECTION = '%s'" % collection

            where_clauses.append(clause)

        # Select by MATRIX.NAME
        if tf_name:
            if isinstance(tf_name, list):
                # Multiple names passed in as a list
                clause = "m.NAME in ('"
                clause = "".join([clause, "','".join(tf_name)])
                clause = "".join([clause, "')"])
            else:
                # A single name
                clause = "m.NAME = '%s'" % tf_name

            where_clauses.append(clause)

        # Select by MATRIX_SPECIES.TAX_ID
        if species:
            tables.append("MATRIX_SPECIES ms")
            where_clauses.append("m.ID = ms.ID")

            """
            NOTE: species are numeric taxonomy IDs but stored as varchars
            in the DB.
            """
            if isinstance(species, list):
                # Multiple tax IDs passed in as a list
                clause = "ms.TAX_ID in ('"
                clause = "".join([clause, "','".join(str(s) for s in species)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "ms.TAX_ID = '%s'" % species

            where_clauses.append(clause)

        """
        Tag based selection from MATRIX_ANNOTATION
        Differs from perl TFBS module in that the matrix class explicitly
        has a tag attribute corresponding to the tags in the database. This
        provides tremendous flexibility in adding new tags to the DB and
        being able to select based on those tags with out adding new code.
        In the JASPAR Motif class we have elected to use specific attributes
        for the most commonly used tags and here correspondingly only allow
        selection on these attributes.

        The attributes corresponding to the tags for which selection is
        provided are:

           Attribute   Tag
           tf_class    class
           tf_family   family
           pazar_id    pazar_tf_id
           medline     medline
           data_type   type
           tax_group   tax_group
        """

        # Select by TF class(es) (MATRIX_ANNOTATION.TAG="class")
        if tf_class:
            tables.append("MATRIX_ANNOTATION ma1")
            where_clauses.append("m.ID = ma1.ID")

            clause = "ma1.TAG = 'class'"
            if isinstance(tf_class, list):
                # A list of TF classes
                clause = "".join([clause, " and ma1.VAL in ('"])
                clause = "".join([clause, "','".join(tf_class)])
                clause = "".join([clause, "')"])
            else:
                # A single TF class
                clause = "".join([clause, " and ma1.VAL = '%s' " % tf_class])

            where_clauses.append(clause)

        # Select by TF families (MATRIX_ANNOTATION.TAG="family")
        if tf_family:
            tables.append("MATRIX_ANNOTATION ma2")
            where_clauses.append("m.ID = ma2.ID")

            clause = "ma2.TAG = 'family'"
            if isinstance(tf_family, list):
                # A list of TF families
                clause = "".join([clause, " and ma2.VAL in ('"])
                clause = "".join([clause, "','".join(tf_family)])
                clause = "".join([clause, "')"])
            else:
                # A single TF family
                clause = "".join([clause, " and ma2.VAL = '%s' " % tf_family])

            where_clauses.append(clause)

        # Select by PAZAR TF ID(s) (MATRIX_ANNOTATION.TAG="pazar_tf_id")
        if pazar_id:
            tables.append("MATRIX_ANNOTATION ma3")
            where_clauses.append("m.ID = ma3.ID")

            clause = "ma3.TAG = 'pazar_tf_id'"
            if isinstance(pazar_id, list):
                # A list of PAZAR IDs
                clause = "".join([clause, " and ma3.VAL in ('"])
                clause = "".join([clause, "','".join(pazar_id)])
                clause = "".join([clause, "')"])
            else:
                # A single PAZAR ID
                clause = "".join([" and ma3.VAL = '%s' " % pazar_id])

            where_clauses.append(clause)

        # Select by PubMed ID(s) (MATRIX_ANNOTATION.TAG="medline")
        if medline:
            tables.append("MATRIX_ANNOTATION ma4")
            where_clauses.append("m.ID = ma4.ID")

            clause = "ma4.TAG = 'medline'"
            if isinstance(medline, list):
                # A list of PubMed IDs
                clause = "".join([clause, " and ma4.VAL in ('"])
                clause = "".join([clause, "','".join(medline)])
                clause = "".join([clause, "')"])
            else:
                # A single PubMed ID
                clause = "".join([" and ma4.VAL = '%s' " % medline])

            where_clauses.append(clause)

        # Select by data type(s) used to compile the matrix
        # (MATRIX_ANNOTATION.TAG="type")
        if data_type:
            tables.append("MATRIX_ANNOTATION ma5")
            where_clauses.append("m.ID = ma5.ID")

            clause = "ma5.TAG = 'type'"
            if isinstance(data_type, list):
                # A list of data types
                clause = "".join([clause, " and ma5.VAL in ('"])
                clause = "".join([clause, "','".join(data_type)])
                clause = "".join([clause, "')"])
            else:
                # A single data type
                clause = "".join([" and ma5.VAL = '%s' " % data_type])

            where_clauses.append(clause)

        # Select by taxonomic supergroup(s) (MATRIX_ANNOTATION.TAG="tax_group")
        if tax_group:
            tables.append("MATRIX_ANNOTATION ma6")
            where_clauses.append("m.ID = ma6.ID")

            clause = "ma6.TAG = 'tax_group'"
            if isinstance(tax_group, list):
                # A list of tax IDs
                clause = "".join([clause, " and ma6.VAL in ('"])
                clause = "".join([clause, "','".join(tax_group)])
                clause = "".join([clause, "')"])
            else:
                # A single tax ID
                clause = "".join([clause, " and ma6.VAL = '%s' " % tax_group])

            where_clauses.append(clause)

        sql = "".join(["select distinct(m.ID) from ", ", ".join(tables)])

        if where_clauses:
            sql = "".join([sql, " where ", " and ".join(where_clauses)])

        # print("sql = %s" % sql)

        cur.execute(sql)
        rows = cur.fetchall()

        for row in rows:
            id = row[0]
            if all_versions:
                int_ids.append(id)
            else:
                # is the latest version?
                if self._is_latest_version(id):
                    int_ids.append(id)

        if len(int_ids) < 1:
            warnings.warn(
                "Zero motifs returned with current select critera", BiopythonWarning
            )

        return int_ids

    def _is_latest_version(self, int_id):
        """Check if the internal ID represents the latest JASPAR matrix (PRIVATE).

        Does this internal ID represent the latest version of the JASPAR
        matrix (collapse on base ids)
        """
        cur = self.dbh.cursor()

        cur.execute(
            "select count(*) from MATRIX where "
            "BASE_ID = (select BASE_ID from MATRIX where ID = %s) "
            "and VERSION > (select VERSION from MATRIX where ID = %s)",
            (int_id, int_id),
        )

        row = cur.fetchone()

        count = row[0]

        if count == 0:
            # no matrices with higher version ID and same base id
            return True

        return False
