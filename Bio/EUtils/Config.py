"""Configuration information about NCBI's databases"""

# Used to figure out if efetch supports the start, stop, strand, and
# complexity fields
PUBLICATION_TYPE = 0
SEQUENCE_TYPE = 1

# Map from database name to database type
class DatabaseInfo:
    """stores NCBI's name for the database and its type"""
    def __init__(self, db, dbtype):
        self.db = db
        self.dbtype = dbtype

class DatabaseDict(dict):
    """map from name to DatabaseInfo for that database name

    Entries are also available through attributes like PUBMED,
    OMIM, and NUCLEOTIDE.
    """
    def gettype(self, db, dbtype = None):
        """Given a database name and optional type, return the database type"""
        if dbtype not in (None, SEQUENCE_TYPE, PUBLICATION_TYPE):
            raise TypeError("Unknown database type: %r" % (dbtype,))
        if dbtype is None:
            dbtype = self[db].dbtype
        return dbtype

databases = DatabaseDict()

def _add_db(x):
    databases[x.db] = x
    return x.db

# XXX Try these
# <option value="structure">Structure</option>
# <option value="pmc">PMC</option>
# <option value="taxonomy">Taxonomy</option>
# <option value="books">Books</option>
# <option value="geo">ProbeSet</option>
# <option value="domains">3D Domains</option>
# <option value="UniSts">UniSTS</option>
# <option value="cdd">Domains</option>
# <option value="snp">SNP</option>
# <option value="popset">PopSet</option>

databases.PUBMED = _add_db(DatabaseInfo("pubmed", 0))
databases.OMIM = _add_db(DatabaseInfo("omim", 0))
databases.JOURNALS = _add_db(DatabaseInfo("journals", 0))
                
databases.GENOME = _add_db(DatabaseInfo("genome", 1))
databases.NUCLEOTIDE = _add_db(DatabaseInfo("nucleotide", 1))
databases.PROTEIN = _add_db(DatabaseInfo("protein", 1))
databases.POPSET = _add_db(DatabaseInfo("popset", 1))
databases.SEQUENCES = _add_db(DatabaseInfo("sequences", 1))
databases.UNIGENE = _add_db(DatabaseInfo("unigene", 1))
databases.GENE = _add_db(DatabaseInfo("gene", 1))


# Someday I want to make it easier to get a given format.  I would
# rather not have to specify the retmode/rettype pair, but I don't
# know what people want from this feature, so skip for now.  Plus,
# it's harder than I thought.

##class FormatInfo:
##    def __init__(self, name, retmode):
##        self.name = name
##        self.retmode = retmode
##        self.rettype = rettype
