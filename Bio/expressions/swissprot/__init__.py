# This is a Python module.

import Martel
from Martel import RecordReader

# Very simple filter
filter = Martel.SimpleRecordFilter(Martel.Str("ID  "),
                                   RecordReader.CountLines, (1,))
