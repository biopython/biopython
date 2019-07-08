#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
"""Update the Rebase EMBOSS files.

The Rebase EMBOSS files are used by ``ranacompiler.py`` to build the updated
``Restriction_Dictionary.py`` module for ``Bio.Restriction``.

"""

from __future__ import print_function

import os
from datetime import date

try:
    # Python 2
    from urllib import urlretrieve, urlcleanup
except ImportError:
    # Python 3
    from urllib.request import urlretrieve, urlcleanup


# Rebase ftp location, do not modify these addresses:
ftp_Rebase = 'ftp://ftp.neb.com/'
ftp_emb_e = ftp_Rebase + 'pub/rebase/emboss_e.###'
ftp_emb_s = ftp_Rebase + 'pub/rebase/emboss_s.###'
ftp_emb_r = ftp_Rebase + 'pub/rebase/emboss_r.###'

# Generate 'time stamp' of type ymm to add to Rebase file names.
# This is the 3 digit number REBASE release number (e.g. 312).
# The first digit is the last digit of the year (e.g. 3 for 2013)
# and the two last the month (e.g. 12 for December)
release_number = date.today().strftime("%y%m")[1:]

# Replace '###' with the 'time stamp'
files = [x.replace('###', release_number)
         for x in [ftp_emb_e, ftp_emb_s, ftp_emb_r]]


def get_files():
    """Download Rebase files."""
    for file in files:
        print('copying %s' % file)
        fn = os.path.basename(file)
        filename = os.path.join(os.getcwd(), fn)
        print('to %s' % filename)
        try:
            urlretrieve(file, filename)
            urlcleanup()
        except IOError as e:
            print(e)
            print('Download of Rebase files failed. Please download the files '
                  '"emboss_e.{0}", "emboss_s.{0}" and "emboss_r.{0}" manually '
                  'from: ftp://ftp.neb.com/pub/rebase.'
                  .format(release_number))
            return
    return


if __name__ == "__main__":
    get_files()
