"""A client-side library for the Entrez databases at NCBI (DEPRECATED).

Please use Bio.Entrez instead.
"""
import re
import warnings
warnings.warn("Bio.UEtils has been deprecated. Please use Bio.Entrez instead,"\
              " which is described in the tutorial.", DeprecationWarning)

__version__ = "1.0p1"
__authors__ = ["Andrew Dalke, Dalke Scientific Software"]

import Datatypes
from Datatypes import DateRange, WithinNDays, EUtilsError, \
     EUtilsSearchError

from Config import databases

def DBIds(db, ids):
    if isinstance(ids, str):
        ids = [ids]
    db = db.lower()
    if db not in databases:
        raise TypeError("Unknown database: %s" % (db,))
    return Datatypes.DBIds(db, ids)

PLUS_STRAND = 1
MINUS_STRAND = 2
