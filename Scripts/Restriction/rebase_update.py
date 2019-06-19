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

The Rebase EMBOSS files are used by `ranacompiler.py` to build the updated
`Restriction_Dictionary.py` module for  `Bio.Restriction`.
"""

from __future__ import print_function

import optparse
import os
import sys
import time
from ftplib import FTP

from Bio.Restriction.RanaConfig import ftp_proxy, ftp_Rebase, Rebase_name
from Bio.Restriction.RanaConfig import ftp_emb_e, ftp_emb_s, ftp_emb_r


def localtime():
    """Generate 'time stamp' of type ymm to add to Rebase file names."""
    t = time.gmtime()
    year = str(t.tm_year)[-1]
    month = str(t.tm_mon)
    if len(month) == 1:
        month = '0' + month
    return year + month


class RebaseUpdate(FTP):
    """A class to fetch the Rebase EMBOSS files."""

    def __init__(self, ftpproxy=''):
        """RebaseUpdate([ftpproxy]]) -> new RebaseUpdate instance.

        if ftpproxy is not given RebaseUpdate uses the corresponding
        variable from RanaConfig.

        ftpproxy is the proxy to use if any.
        """
        proxy = {'ftp': ftpproxy or ftp_proxy}
        if not proxy['ftp']:
            FTP.__init__(self, host=ftp_Rebase)
        else:
            FTP.__init__(self, host=proxy)

    def getfiles(self, *files):
        """Download Rebase files."""
        for file in self.update(*files):
            print('copying %s' % file)
            fn = os.path.basename(file)
            filename = os.path.join(os.getcwd(), fn)
            print('to %s' % filename)
            try:
                self.retrbinary('RETR ' + file, open(filename, 'wb').write)
            except IOError as e:
                print(e)
                self.close()
                return
        self.close()
        return

    def update(self, *files):
        """Update filenames to recent versions (indicated by 'time stamp')."""
        if not files:
            files = [ftp_emb_e, ftp_emb_s, ftp_emb_r]
        return [x.replace('###', localtime()) for x in files]


if __name__ == '__main__':
    parser = optparse.OptionParser()
    add = parser.add_option

    add('-p', '--proxy',
        action="store",
        dest='ftp_proxy',
        default='',
        help="set the proxy to be used by the ftp connection.")

    (option, args) = parser.parse_args()

    Getfiles = RebaseUpdate(option.ftp_proxy)
    Getfiles.openRebase()
    Getfiles.getfiles()
    Getfiles.close()
    sys.exit()
