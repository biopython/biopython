# BLASTN 2.0a19MP-WashU [05-Feb-1998] [Build decunix3.2 01:53:29 05-Feb-1998]
# BLASTP 2.0.4 [Feb-24-1998]
class Algorithm(object):
    def __init__(self, name, version, description = ""):
        self.name = name                 # 'blastx', 'blastn', etc.
        self.version = version           # '2.1.2' or '2.0a19MP-WashU'
        self.description = description   # '[05-Feb-1998] [Build dec ...1998]'

# Query=  YAL001C YAL001C, Chr I from 147596 to 147665, and 147756 to 151168,
#     reverse complement
#         (3483 letters)
class Query(object):
    def __init__(self, name, accession, description, length):
        self.name = name                # 'YAL001C'
        self.accession = accession      # or None if missing
        self.description = description  # 'YAL001C, Chr I from 147596 to ... '
        self.length = length            # 3483

# Database:  ArabidopsisN
#            66,211 sequences; 69,074,155 total letters.
class Database(object):
    def __init__(self, name, letters, entries):
        self.name = name            # ArabidopsisN
        self.letters = letters      # 69074155
        self.entries = entries      # 66211

class TableInfo(object):
    def __init__(self, full_description, info):
        self.__dict__.update(info)
        self.full_description = full_description


class Search(object):
    def __init__(self, algorithm, query, database, table, hits,
                 parameters, statistics):
        self.algorithm = algorithm
        self.query = query
        self.database = database
        self.table = table
        self.hits = hits
        self.parameters = parameters
        self.statistics = statistics

class Hit(object):
    def __init__(self, name, description, accession, length,
                 algorithm, hsps = None):
        self.name = name
        self.description = description
        self.accession = accession
        self.length = length
        self.algorithm = algorithm
        if hsps is None:
            hsps = []
        self.hsps = hsps

    def __len__(self):
        return self.length



# >GB_PL:ATF18F4 AL021637 Arabidopsis thaliana DNA chromosome 4, BAC clone 
#           F18F4 (ESSAII project). 2/98
#             Length = 93,646
#  
#   Minus Strand HSPs:
#  
#  Score = 226 (33.9 bits), Expect = 0.80, P = 0.55
#  Identities = 98/142 (69%), Positives = 98/142 (69%), Strand = Minus / Plus
#    [...lines deleted...]
# Query:  2486 ATATCAAGCAATTTGATAAGATCTAG 2461
#              A AT  A C ATT GA AAGATC AG
# Sbjct: 85387 AGATTTACCTATT-GAGAAGATCAAG 85411

# computed from the strings
class _SeqLength:
    def __init__(self, length, identical, positives, gaps):
        self.length = length
        self.identical = identical
        self.positives = positives
        self.gaps = gaps
    def __len__(self):
        return self.length
    def __getattr__(self, name):
        if name == "frac_identical":
            return float(self.identical) / self.length
        elif name == "frac_positives":
            return float(self.positives) / self.length
        raise AttributeError(name)


class HomologySeq(_SeqLength):
    def __init__(self, seq, identical, positives, gaps):
        _SeqLength.__init__(self, len(seq), identical, positives, gaps)
        self.seq = seq

class HSPSeq(_SeqLength):
    def __init__(self, name, seq, location, identical, positives, gaps):
        _SeqLength.__init__(self, len(seq), identical, positives, gaps)
        self.name = name
        self.seq = seq
        self.location = location
        

class HSP(_SeqLength):
    def __init__(self,
                 query_seq,    # ATATCAAGCAATTTGATAAGATCTAG
                 homology_seq, # A AT  A C ATT GA AAGATC AG
                 subject_seq,  # AGATTTACCTATT-GAGAAGATCAAG

                 query_location,   # (2486, 2461, negative strand)
                 subject_location, # (85387, 85411)

                 query_name,     # Query (or None)
                 subject_name,   # Sbjct (or None)

                 algorithm,  # an Algorithm
                 info,       # contains Key/value pairs
                 homology_gaps = None,  # Is this needed?
                 ):
        assert len(query_seq) == len(homology_seq) == len(subject_seq), \
               (query_seq, homology_seq, subject_seq)
        self.algorithm = algorithm

        query_gaps = query_seq.count("-")
        subject_gaps = subject_seq.count("-")
        if homology_gaps is None:
            homology_gaps = query_gaps + subject_gaps
        self.info = info

        identical = info["identical"]
        # bioperl calls this 'conserved'
        positives = info.get("positives", identical)
        
        _SeqLength.__init__(self, len(query_seq), identical,
                            positives, homology_gaps)

        self.query = HSPSeq(name = query_name,
                            seq = query_seq,
                            location = query_location,
                            identical = identical,
                            positives = positives,
                            gaps = query_gaps)

        self.subject = HSPSeq(name = subject_name,
                              seq = subject_seq,
                              location = subject_location,
                              identical = identical,
                              positives = positives,
                              gaps = subject_gaps)
        self.homology = HomologySeq(seq = homology_seq,
                                    identical = identical,
                                    positives = positives,
                                    gaps = homology_gaps)
