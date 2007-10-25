"""A Martel format to parse the output from transfac.

Formats:
format             Format for a whole file.

"""

import warnings
warnings.warn("Bio.expressions was deprecated, as it does not work with recent versions of mxTextTools. If you want to continue to use this module, please get in contact with the Biopython developers at biopython-dev@biopython.org to avoid permanent removal of this module from Biopython", DeprecationWarning)


import sys

from Martel import *
from Martel import RecordReader

blank_line = Opt(Spaces()) + AnyEol()

MATRIX_LINE = Str("Search for sites by WeightMatrix library:") + Spaces() + \
              UntilEol("matrix_file") + AnyEol()
SEQUENCE_LINE = Str("Sequence file:") + Spaces() + \
                UntilEol("sequence_file") + AnyEol()
PROFILE_LINE = Str("Site selection profile:") + Spaces() + \
               UntilSep("profile_file", sep=" ") + Spaces() + \
               UntilEol("profile_description") + AnyEol()

TITLE_LINE = Str("Inspecting sequence ID") + Spaces() + \
             UntilSep("entryname", sep=" ") + Spaces() + \
             UntilSep("dataclass", sep=";") + Str(";") + Spaces() + \
             UntilSep("molecule", sep=";") + Str(";") + Spaces() + \
             UntilSep("division", sep=";") + Str(";") + Spaces() + \
             UntilSep("sequencelength", sep=" ") + Spaces() + Str("BP") + \
             UntilEol() + AnyEol()

def SS(exp):  # expression surrounded by optional spaces.
    return Opt(Spaces()) + exp + Opt(Spaces())


DATA_LINE = \
    SS(UntilSep("matrix_identifier", sep=" |")) + \
    Str("|") + \
    SS(UntilSep("position", sep=" ")) + \
    SS(Str("(") + Group("strand", Any("+-")) + Str(")")) + \
    Str("|") + \
    SS(Float("core_match")) + \
    Str("|") + \
    SS(Float("matrix_match")) + \
    Str("|") + \
    Opt(Spaces()) + UntilEol("sequence") + AnyEol()

SEQUENCES_LENGTH_LINE = \
    Spaces() + Str("Total sequences length=") + Integer("sequences_length") + \
    AnyEol()

FOUND_SITES_LINE = \
    Spaces() + Str("Total number of found sites=") + Integer("found_sites") + \
    AnyEol()

SITE_FREQUENCY_LINE = \
    Spaces() + Str("Frequency of sites per nucleotide=") + \
    Float("sites_per_nucleotide") + AnyEol()
    
format = MATRIX_LINE + \
         SEQUENCE_LINE + \
         PROFILE_LINE + \
         blank_line + \
         TITLE_LINE + \
         blank_line + \
         Rep(DATA_LINE) + \
         blank_line + \
         SEQUENCES_LENGTH_LINE + \
         blank_line + \
         FOUND_SITES_LINE + \
         blank_line + \
         SITE_FREQUENCY_LINE 
