# soon-to-be deprecated. Use Bio.Data.SCOPData instead!

import warnings

from Bio.Data.SCOPData import protein_letters_3to1 as to_one_letter_code
from Bio import BiopythonDeprecationWarning


warnings.warn("The Bio.SCOP.three_to_one_dict has been moved to "
        "Bio.Data.SCOPData and renamed to protein_letters_3to1. "
        "Bio.SCOP.three_to_one_dict will be deprecated and removed in "
        "future versions of Biopython.", BiopythonDeprecationWarning)
