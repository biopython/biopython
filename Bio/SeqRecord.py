# Stores data about the sequence

# I HAVEN'T REALLY THOUGHT ABOUT THIS DESIGN!!!

# There's only enough here for me to bootstrap reading FASTA
# sequences.

# There will probably be a restricted vocabulary on the features and
# annotations, and one of them will be a "vocabulary" to distinguish
# between differences in vocabulary. :)

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
        
