# This is a Python module.
"""This module is DEPRECATED.

Andrew Dalke is no longer maintaining Martel or Bio.Mindy, and these modules
and associate ones like Bio.builders are now deprecated.  They are no longer
used in any of the current Biopython parsers, and are likely to be removed
in a future release.
"""

import warnings
warnings.warn("Martel and those parts of Biopython depending on it" \
              +" directly (such as Bio.Mindy and Bio.builders) are now" \
              +" deprecated, and will be removed in a future release of"\
              +" Biopython.  If you want to continue to use this code,"\
              +" please get in contact with the Biopython developers via"\
              +" the mailing lists to avoid its permanent removal from"\
              +" Biopython.", \
              DeprecationWarning)
