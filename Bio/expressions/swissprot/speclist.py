"""Parser for the SPECLIST.TXT file in SWISS-PROT.

You probably want to use the variables 'record' (for a single record)
and 'format' (for a set of records).

"""
from Martel import *


def SkipLinesToNoEOL(expr):
    return Rep(AssertNot(expr) + ToEol()) + expr
    


DESCRIPTION = Str("Description:") + Spaces() + UntilEol("description") + AnyEol()
NAME = Str("Name:") + Spaces() + UntilEol("name") + AnyEol()
RELEASE = Str("Release:") + Spaces() + UntilEol("release") + AnyEol()

TOTAL_CODES = Str("Total number of organism identification codes currently defined:") + Spaces() + Integer("num_organism_codes") + Str(".") + AnyEol()

TABLE_HEADER = Group("table_header",
                     Str("Code") + Spaces() + Str("Taxon") + Spaces() +
                     Str("N=Official name") + ToEol() +
                     Spaces() + Str("Node") + Spaces() +
                     Str("C=Common name") + ToEol() +
                     Spaces() + Str("S=Synonym") + ToEol() +
                     Rep1(Any(" _")) + AnyEol()
                     )

_dash_line = Rep1(Str("-")) + AnyEol()
COPYRIGHT = Group("copyright",
                  _dash_line +
                  Str("SWISS-PROT is copyright.") + ToEol() +
                  SkipLinesToNoEOL(_dash_line) +
                  Rep(AnyEol())
                  )

_code_line = Group("code", Re(r"[A-Z0-9]{1,5}")) + \
             Spaces() + \
             Group("kingdom", Any("ABEV")) + \
             Spaces() + \
             Group("taxon_node", Digits()) + Str(":") + \
             Spaces() + \
             Str("N=") + Group("official_name", UntilEol()) + \
             AnyEol()
_common_line = Spaces() + \
               Str("C=") + Group("common_name", UntilEol()) + \
               AnyEol()
_synonym_line = Spaces() + \
                Str("S=") + Group("synonym", UntilEol()) + \
                AnyEol()
record = Group("record",
               _code_line + Rep(_common_line) + Rep(_synonym_line)
               )

format = Group("format",
               SkipLinesToNoEOL(DESCRIPTION) +
               NAME +
               RELEASE +
               SkipLinesToNoEOL(TOTAL_CODES) +
               SkipLinesToNoEOL(TABLE_HEADER) +
               Rep1(record) +
               AnyEol() +
               COPYRIGHT
               )
