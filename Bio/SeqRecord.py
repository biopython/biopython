# Stores data about the sequence

# I HAVEN'T REALLY THOUGHT ABOUT THIS DESIGN!!!

# There's only enough here for me to bootstrap reading FASTA
# sequences.

# There will probably be a restricted vocabulary on the features and
# annotations, and one of them will be a "vocabulary" to distinguish
# between differences in vocabulary. :)

try:
    from Bio import FormatIO

    # Should this be in the module namespace or the record namespace?
    io = FormatIO.FormatIO("SeqRecord",
                           default_input_format = "sequence",
                           default_output_format = "fasta")
except SyntaxError: # we will get an error right now without python2.2
    import warnings # requires python 2.1, so we aren't doing much better :-)
    warnings.warn("No IO support because your python is too old.")


class SeqRecord:                           
    def __init__(self, seq, id = "<unknown id>", name = "<unknown name>",
                 description = "<unknown description>"):
        self.seq = seq
        self.id = id
        self.name = name
        self.description = description
        # annotations about the whole sequence
        self.annotations = {}
        
        # annotations about parts of the sequence
        self.features = []
        
