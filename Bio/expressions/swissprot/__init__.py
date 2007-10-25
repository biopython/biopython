# This is a Python module.

import warnings
warnings.warn("Bio.expressions was deprecated, as it does not work with recent versions of mxTextTools. If you want to continue to use this module, please get in contact with the Biopython developers at biopython-dev@biopython.org to avoid permanent removal of this module from Biopython", DeprecationWarning)


import Martel
from Martel import RecordReader

# Very simple filter
filter = Martel.SimpleRecordFilter(Martel.Str("ID  "),
                                   RecordReader.CountLines, (1,))
