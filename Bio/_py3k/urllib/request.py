# Copyright 2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Python 3 urllib compatibility module (PRIVATE)."""

import sys
if sys.version_info[0] >= 3:
    from urllib.request import urlretrieve, urlcleanup, pathname2url, url2pathname
    from urllib.request import getproxies, URLopener, FancyURLopener
    from urllib.request import urlopen, install_opener, build_opener
    from urllib.request import Request, OpenerDirector, BaseHandler
    #...
    from urllib.request import HTTPHandler, FileHandler, FTPHandler
    from urllib.request import CacheFTPHandler, UnknownHandler
else:
    from urllib import urlretrieve, urlcleanup, pathname2url, url2pathname
    from urllib import getproxies, URLopener, FancyURLopener
    from urllib2 import urlopen, install_opener, build_opener
    from urllib2 import Request, OpenerDirector, BaseHandler
    #...
    from urllib2 import HTTPHandler, FileHandler,FTPHandler
    from urllib2 import CacheFTPHandler,UnknownHandler

