# This was tested against BLOCKS-12.0, June 2000

import warnings
warnings.warn("Bio.expressions was deprecated, as it does not work with recent versions of mxTextTools. If you want to continue to use this module, please get in contact with the Biopython developers at biopython-dev@biopython.org to avoid permanent removal of this module from Biopython", DeprecationWarning)


from Martel import *
from Martel import RecordReader
from Bio import Std

# Header goes up to the line starting with "ID"
header = Rep(AssertNot(Str("ID   ")) + \
             ToEol())

# ID   kringle; BLOCK
# ID   14-3-3; BLOCK
#  but not!
# IDSA_METJA|Q58270  (  46) GGKRIRPYLTV  11
ID = Str("ID   ") + Std.dbid(ToSep(sep = ";"), {"type": "primary"}) + \
     Str(" BLOCK") + AnyEol()

# AC   IPB000001A; distance from previous block=(10,266)
AC = Str("AC   ") + Std.dbid(ToSep(sep = ";"), {"type": "accession"}) + \
     Str(" distance from previous block=(") + \
     Integer("dist1") + Str(",") + Integer("dist2") + \
     Str(")") + AnyEol()


# DE   Kringle domain
#  If the DE line is long, it doen't fold .. it's all on one line
DE = Str("DE   ") + ToEol("description")


# BL   CCY;  width=14; seqs=44; 99.5%=717; strength=1059
BL = Str("BL   ") + ToSep("protomat_id", ";") + \
     Str("  width=") + Digits("width") + \
     Str("; seqs=") + Digits("numseqs") + \
     Str("; 99.5%=") + Digits("protomat_count") + \
     Str("; strength=") + Digits("strength") + \
     AnyEol()


# PLMN_BOVIN|P06868  (  60) CEEETDFVCRAFQY  26
# ^^^^^^^^^^^^^^^^^
#                     ^^^^-- number of segments
#                           ^^^^^^^^^^^^^^-- matching sequence
#                                           ^^-- weight
#
identifier = (Std.dbxref_dbid(UntilSep(sep = "|."), 
                          {"dbname": "swissprot", "type": "primary"}) + \
              Str("|") + \
              Std.dbxref_dbid(UntilSep(sep = " "), 
                              {"dbname": "swissprot", "type": "accession"})) |\
              Std.dbxref_dbid(UntilSep(sep = " "))
                              
segment = AssertNot(Re(r".. ")) + \
          identifier + \
          Re(r" *\( *") + \
          Integer("position") + \
          Re(r"\) *") + \
          Word("matching_sequence") + Spaces() + \
          Digits("weight") + AnyEol()

segment_block = Rep1(segment | AnyEol())

end = Str("//") + AnyEol()

record = Group("record",
               ID + AC + DE + BL + segment_block + end)

format_expression = header + Rep1(record)
format = HeaderFooter("dataset", {"format": "blocks/12"},
                      header, RecordReader.Until, ("ID ",),
                      record, RecordReader.EndsWith, ("//\n",),
                      None, None, None)

if __name__ == "__main__":
    import os
    from xml.sax import saxutils
    filename = "/home/dalke/ftps/databases/blocks/unix/blocks-12.0/blocks.dat.Z"
    infile = os.popen("zcat " + filename)
    parser = format.make_parser(debug_level = 0)
    parser.setContentHandler(saxutils.XMLGenerator())
    parser.parseFile(infile)
