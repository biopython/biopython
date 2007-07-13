#!/usr/bin/env python
#
#      Restriction Analysis Libraries.
#      Copyright (C) 2004. Frederic Sohm.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import sys
import optparse
from Bio.Restriction._Update.Update import *

        
if __name__ == '__main__' :
    parser = optparse.OptionParser()
    add = parser.add_option

    add('-m', '--e-mail',
        action  = "store",
        dest    = 'rebase_password',
        default = '',
        help    = "set the e-mail address to be used as password for the"
        "anonymous ftp connection to Rebase.")
    add('-p', '--proxy',
        action  = "store",
        dest    = 'ftp_proxy',
        default = '',
        help    = "set the proxy to be used by the ftp connection.")
    
    (option, args) = parser.parse_args()
    
    Getfiles = RebaseUpdate(option.rebase_password, option.ftp_proxy)
    Getfiles.openRebase()
    Getfiles.getfiles()
    Getfiles.close()
    sys.exit()
