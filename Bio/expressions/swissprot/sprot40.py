import Martel
from Martel import Time
import sprot38

# HAS2_CHICK has a DT line like this
# DT   30-MAY-2000 (REL. 39, Created)
#                   ^^^ Note the upper-case "REL" instead of "Rel" !
DT_created_exp = (Martel.Str("DT   ") +
                  Time.make_expression("%(DD)-%(Jan)-%(YYYY)") + \
                  Martel.Re(" \(R[Ee][Ll]. (?P<release>\d\d), Created\)\R"))


OX_start = (Martel.Str("OX   NCBI_TaxID=") +
            Martel.Rep1(Martel.Digits("ncbi_taxid") +
                        Martel.Re("[,; ]+")) +
            Martel.AnyEol())
OX_cont = (Martel.Str("OX   ") +
           Martel.Rep1(Martel.Digits("ncbi_taxid") +
                       Martel.Re("[,; ]+")) +
           Martel.AnyEol())

OX_exp = OX_start + Martel.Rep(OX_cont)

# 0 or 1
# in 40 the line changed to look like this
#  RX   MEDLINE=93305731; PubMed=7916637;
#  RX   PubMed=11001938;
bib = (Martel.Word("bibliographic_database_name") + Martel.Str("=") +
       Martel.ToSep("bibliographic_identifier", ";")
       )
RX_exp = (Martel.Str("RX   ") + bib +
          Martel.Opt(Martel.Str(" ") + bib) +
          Martel.AnyEol())

# Here's the neq SQ line format -- uses a CRC64
# SQ   SEQUENCE   889 AA;  100368 MW;  ABD7E3CD53961B78 CRC64;
SQ_exp = Martel.Re("SQ   SEQUENCE +(?P<sequence_length>\d+) AA;" \
                   " +(?P<molecular_weight>\d+) MW;" \
                   " +(?P<crc?type=64>\w+) CRC64;\R")

replacements = [
    ("DT_created", DT_created_exp),
    ("OX_block", OX_exp),
    ("RX", RX_exp),
    ("SQ", SQ_exp),
    ]
record = Martel.replace_groups(sprot38.record, replacements)


format_expression = Martel.replace_groups(
    sprot38.format_expression, replacements)


format = Martel.replace_groups(sprot38.format, replacements)

if __name__ == "__main__":
    parser = format.make_parser()
    filename = "/home/dalke/ftps/databases/swiss-prot/release_compressed/sprot40.dat"
##    import os
##    infile = os.popen("zcat " + filename)
    infile = open(filename)
    infile.seek(107976062)
    parser.parseFile(infile)
    
