# Stores data about the sequence

# NEEDS TO BE SYNCH WITH THE REST OF BIOPYTHON AND BIOPERL

import sys
if getattr(sys, "version_info", (1, 5))[:2] >= (2,1):
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
        
