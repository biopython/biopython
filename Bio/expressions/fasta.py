# Read a FASTA description

import operator
from Martel import *
from Martel import RecordReader
from Bio import Std

### Parse dbxrefs given the NCBI|descr|line as explained
### in ftp://ncbi.nlm.nih.gov/blast/db/README and augmented
### by experience

def make_2id(s, dbname, primary_name, secondary_name):
    assert secondary_name is not None
    if primary_name is None:
        return Str(s + "||") + \
               Std.dbxref_dbid(UntilSep(sep = "| "), {"dbname": dbname,
                                                      "type": secondary_name})

    return Str(s + "|") + \
           Std.dbxref_dbid(UntilSep(sep = "|"), {"dbname": dbname,
                                                 "type": primary_name}) + \
           Str("|") + \
           Std.dbxref_dbid(UntilSep(sep = "| "), {"dbname": dbname,
                                                  "type": secondary_name})

def make_1id(s, dbname, name):
    return Str(s + "|") + \
           Std.dbxref_dbid(UntilSep(sep = "| "), {"dbname": dbname,
                                                  "type": name})
                                    
ids = []
# gene identifier              gi|id  # This isn't in the README
ids.append(make_1id("gi", "x-gi", "primary"))

# GenBank                      gb|accession|locus
# gb|U37104|APU37104
ids.append(make_2id("gb", "gb", "primary", "secondary"))

# EMBL Data Library            emb|accession|locus
# emb|F19596|HSPD04201
ids.append(make_2id("emb", "embl", "primary", "secondary"))

# DDBJ, DNA Database of Japan  dbj|accession|locus
ids.append(make_2id("dbj", "ddbj", "primary", "secondary"))

# NBRF PIR                     pir||entry
ids.append(make_2id("pir", "pir", None, "primary"))

# Protein Research Foundation  prf||name
ids.append(make_2id("prf", "x-prf", None, "primary"))

# SWISS-PROT                   sp|accession|entry name
ids.append(make_2id("sp", "sp", "primary", "secondary"))

# Brookhaven Protein Data Bank pdb|entry|chain
ids.append(make_2id("pdb", "x-pdb", "primary", "secondary"))  # XXX not correct

# Patents                      pat|country|number 
ids.append(make_2id("pat", "x-pat", "primary", "secondary"))  # XXX not correct

# GenInfo Backbone Id          bbs|number 
ids.append(make_1id("bbs", "x-bbs", "primary"))

# General database identifier  gnl|database|identifier
gnl_id = Str("gnl|") + \
         Std.dbxref_dbname(UntilSep(sep = "| ")) + \
         Str("|") + \
         Std.dbxref_dbid(UntilSep(sep = "| "))


# NCBI Reference Sequence      ref|accession|locus
ids.append(make_2id("ref", "x-ref", "primary", "secondary"))

# Local Sequence identifier    lcl|identifier
ids.append(make_1id("lcl", "local", "primary"))

# "|" them all together
ncbi_word = Std.dbxref(reduce(operator.or_, ids))

#ncbi_term = Assert(Re("[^ \R]+\|")) + \
ncbi_term =  ncbi_word + Rep(Str("|") + ncbi_word)

# Anything else
generic_term = Std.dbxref(
                 Std.dbxref_dbid(UntilSep(sep = " "), {"dbname": "local"})
               )
id_term = ncbi_term | generic_term
###########################################################

comment = Str(">") + Std.description_line(id_term + UntilEol()) + AnyEol()
seqline = AssertNot(Str(">")) + Std.sequence(UntilEol()) + AnyEol()

sequence = Std.sequence_block(Rep(seqline))

record = Std.record(comment + sequence + Rep(AnyEol()))

format = ParseRecords("dataset", {"format": "fasta"}, record,
                      RecordReader.StartsWith, (">",) )
