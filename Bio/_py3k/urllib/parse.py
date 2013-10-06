# Copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Python 3 urllib compatibility module (PRIVATE)."""

import sys
if sys.version_info[0] >= 3:
    from urllib.parse import quote, quote_plus, unquote, unquote_plus, urlencode
else:
    from urllib import quote, quote_plus, unquote, unquote_plus, urlencode


