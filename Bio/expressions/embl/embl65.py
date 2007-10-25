"""Format from EMBL Nucleotide Sequence Database Release 65, December 2000

"""

import warnings
warnings.warn("Bio.expressions was deprecated, as it does not work with recent versions of mxTextTools. If you want to continue to use this module, please get in contact with the Biopython developers at biopython-dev@biopython.org to avoid permanent removal of this module from Biopython", DeprecationWarning)



import Martel
from Martel import RecordReader, Time
from Bio import Std

from Bio.expressions.swissprot import sprot38

whitespace = Martel.Spaces()

## ID - identification             (begins each entry; 1 per entry)
# ID   entryname  dataclass; molecule; division; sequencelength BP.

divisions = Martel.Re("EST|PHG|FUN|GSS|HTC|HTG|HUM|INV|ORG|MAM|VRT|PLN|" + \
                      "PRO|ROD|SYN|STS|UNC|VRL|[A-Z]{3}")

# XXX is found in S40706
ID_line = Martel.Str("ID   ") + \
          Std.dbid(Martel.UntilSep("entry_name", " "), {"type": "primary",
                                                        "dbname": "embl"}) + \
          whitespace + \
          Martel.ToSep("dataclass", ";") + \
          whitespace + \
          Martel.Group("molecule",
                       Std.alphabet(Martel.Str("DNA", "circular DNA"),
                                    {"alphabet": "iupac-ambiguous-dna"}) |
                       Std.alphabet(Martel.Str("RNA", "circular RNA"),
                                    {"alphabet": "iupac-ambiguous-rna"}) |
                       Std.alphabet(Martel.Str("XXX"),
                                    {"alphabet": "nucleotide"})) + \
          Martel.Str("; ") + \
          Martel.Group("division", divisions) + \
          Martel.Str("; ") + \
          Martel.Digits("length") + \
          Martel.Str(" BP.") + \
          Martel.AnyEol()


## AC - accession number           (>=1 per entry)
accession = Std.dbid(Martel.UntilSep("accession", ";"),
                     {"type": "accession",
                      "dbname": "embl"}) + Martel.Str(";")
AC_line = Martel.Str("AC   ") + \
          accession + Martel.Rep(Martel.Str(" ") + accession) + \
          Martel.AnyEol()
          
AC_block = Martel.Rep1(AC_line)

## SV - sequence version           (1 per entry)
SV_line = Martel.Str("SV   ") + \
          Martel.Group("sequence_version",
                       Martel.ToSep("accession", ".") + \
                       Martel.Digits("version")) + \
          Martel.AnyEol()


## DT - date                       (2 per entry)
date = Time.make_expression("%(day)-%(Jan)-%(year)")

DT_created_line = Martel.Str("DT   ") + \
                  Martel.Group("date_created", date) + \
                  Martel.Str(" (Rel. ") + \
                  Martel.Digits("release_created") + \
                  Martel.Str(", Created)") + \
                  Martel.AnyEol()

DT_updated_line = Martel.Str("DT   ") + \
                  Martel.Group("date_updated", date) + \
                  Martel.Str(" (Rel. ") + \
                  Martel.Digits("release_updated") + \
                  Martel.Str(", Last updated, Version ") + \
                  Martel.Digits("version_number") + \
                  Martel.Str(")") + \
                  Martel.AnyEol()

DT_block = DT_created_line + DT_updated_line

## DE - description                (>=1 per entry)
DE_line = Martel.Str("DE   ") + \
          Std.description(Martel.UntilEol("description")) + \
          Martel.AnyEol()

DE_block = Std.description_block(Martel.Group("description_block",
                                              Martel.Rep1(DE_line)))

## KW - keyword                    (>=1 per entry)
KW_line = Martel.Str("KW   ") + \
          Martel.ToEol("keyword_data")
KW_block = Martel.Rep1(KW_line)

## OS - organism species           (>=1 per entry)
OS_block = sprot38.OS_block

## OC - organism classification    (>=1 per entry)
OC_block = sprot38.OC_block

## OG - organelle                  (0 or 1 per entry)
OG_block = sprot38.OG_block

organism = Martel.Group("organism",
                        OS_block + \
                        OC_block + \
                        Martel.Opt(OG_block))

## RN - reference number           (>=1 per entry)
## RC - reference comment          (>=0 per entry)
## RP - reference positions        (>=1 per entry)
## RX - reference cross-reference  (>=0 per entry)
## RA - reference author(s)        (>=1 per entry)
## RT - reference title            (>=1 per entry)
## RL - reference location         (>=1 per entry)
RN_line = sprot38.RN
RC_block = sprot38.RC_block
RP_line = sprot38.RP

RX_line = sprot38.RX
RX_block = Martel.Group("RX_block", Martel.Rep1(RX_line))

RA_block = sprot38.RA_block
RT_block = sprot38.RT_block
RL_block = sprot38.RL_block

reference = Martel.Group("reference",
                         RN_line + \
                         Martel.Opt(RC_block) + \
                         Martel.Opt(RP_line) + \
                         Martel.Opt(RX_block) + \
                         RA_block + \
                         RT_block + \
                         RL_block)

## DR - database cross-reference   (>=0 per entry)
DR_block = sprot38.DR_block

## FH - feature table header       (0 or 2 per entry)
FH_block = Martel.Str("FH   Key             Location/Qualifiers") + \
           Martel.AnyEol() + \
           Martel.Str("FH") + \
           Martel.AnyEol()

## FT - feature table data         (>=0 per entry)
##FT_line = Martel.Str("FT   ") + \
##          Martel.ToEol("ft_data")
##FT_block = Martel.Rep1(FT_line)

fq_dbxref = Std.feature_qualifier_name(Martel.Str("db_xref")) + \
            Martel.Str('=') + \
            Std.feature_qualifier_description(
                Martel.Str('"') + \
                Std.dbxref(Std.dbxref_dbname(Martel.UntilSep(None, ":")) + \
                           Martel.Str(":") + \
                           Std.dbxref_dbid(Martel.UntilSep(None, '"'))) + \
                Martel.Str('"')) + \
            Martel.AnyEol()
                       

fq_generic = \
           Martel.Assert(Martel.Word() + Martel.Str("=")) + \
           Std.feature_qualifier_name(Martel.Word()) + \
           Martel.Str("=") + \
           Std.feature_qualifier_description(Martel.UntilEol()) + \
           Martel.AnyEol() + \
           Martel.Rep(
               Martel.Str("FT                   ") + \
               (Martel.AssertNot(Martel.Str("/")) |
               Martel.AssertNot(Martel.Re(r"/\w+="))) + \
           Std.feature_qualifier_description(Martel.UntilEol()) + \
               Martel.AnyEol())

feature_qualifier = Std.feature_qualifier(
    Martel.Str("FT                   /") + \
    (fq_dbxref | fq_generic))

feature = Std.feature(
    Martel.Str("FT   ") + \
    Std.feature_name(Martel.UntilSep(sep = " ")) + \
    whitespace + \
    Std.feature_location(Martel.UntilEol()) + \
    Martel.AnyEol() + \
    Martel.Rep(Martel.Str("FT                   ") + \
               Martel.AssertNot(Martel.Str("/")) + \
               Std.feature_location(Martel.UntilEol()) + \
               Martel.AnyEol()
               ) + \
    Martel.Rep(feature_qualifier)
    )
    
FT_block = Std.feature_block(Martel.Rep(feature),
                             {"location-style": "genbank"})

                         
    

## CC - comments or notes          (>=0 per entry)
CC_line = Martel.Str("CC   ") + \
          Martel.ToEol("comment")
CC_block = Martel.Rep1(CC_line)

## XX - spacer line                (many per entry)
XX = Martel.Str("XX") + Martel.AnyEol()

## SQ - sequence header            (1 per entry)
SQ_line = Martel.Str("SQ   Sequence ") + \
          Martel.Digits("num_BP") + \
          Martel.Str(" BP; ") + \
          Martel.Digits("num_A") + \
          Martel.Str(" A; ") + \
          Martel.Digits("num_C") + \
          Martel.Str(" C; ") + \
          Martel.Digits("num_G") + \
          Martel.Str(" G; ") + \
          Martel.Digits("num_T") + \
          Martel.Str(" T; ") + \
          Martel.Digits("num_other") + \
          Martel.Str(" other;") + \
          Martel.AnyEol()

## bb - (blanks) sequence data     (>=1 per entry)
SQ_data = Martel.Str("     ") + \
          Std.sequence(Martel.Re(".{65}")) + \
          whitespace + \
          Martel.Digits("end_position") + \
          Martel.AnyEol()

SQ_block = Std.sequence_block(SQ_line + Martel.Rep1(SQ_data))

## // - termination line           (ends each entry; 1 per entry)
end = Martel.Str("//") + Martel.AnyEol()

record = Martel.Group("record", \
                      ID_line + \
                      Martel.Opt(XX) + \
                      AC_block + \
                      Martel.Opt(XX) + \
                      SV_line + \
                      Martel.Opt(XX) + \
                      DT_block + \
                      Martel.Opt(XX) + \
                      DE_block + \
                      Martel.Opt(XX) + \
                      KW_block + \
                      Martel.Opt(XX) + \
                      Martel.Rep1(organism + Martel.Opt(XX)) + \
                      Martel.Rep(reference + Martel.Opt(XX)) + \
                      Martel.Opt(DR_block + \
                                 Martel.Opt(XX)) + \
                      Martel.Rep(CC_block + \
                                 Martel.Opt(XX)) + \
                      FH_block + \
                      FT_block + \
                      Martel.Opt(XX) + \
                      SQ_block + \
                      end,
                      {"format": "embl/65"})

format_expression = Martel.Group("dataset", Martel.Rep1(record),
                                 {"format": "embl/65"})

format = Martel.ParseRecords("dataset", {"format": "embl/65"},
                             record, RecordReader.EndsWith, ("//\n",) )
