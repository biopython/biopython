# Stores data about the sequence

# I HAVEN'T REALLY THOUGHT ABOUT THIS DESIGN!!!

# There's only enough here for me to bootstrap reading FASTA
# sequences.

# There will probably be a restricted vocabulary on the features and
# annotations, and one of them will be a "vocabulary" to distinguish
# between differences in vocabulary. :)

from Bio import FormatIO

# Should this be in the module namespace or the record namespace?
io = FormatIO.FormatIO("SeqRecord",
                       default_input_format = "sequence",
                       default_output_format = "fasta")


class SeqRecord:                           
    def __init__(self, seq, id = "<unknown id>", name = "<unknown name>",
                 description = "<unknown description>", dbxrefs = None,
                 features = None):
        self.seq = seq
        self.id = id
        self.name = name
        self.description = description
        if dbxrefs is None:
            dbxrefs = []
        self.dbxrefs = dbxrefs
        # annotations about the whole sequence
        self.annotations = {}
        
        
        # annotations about parts of the sequence
        if features is None:
            features = []
        self.features = features
        
