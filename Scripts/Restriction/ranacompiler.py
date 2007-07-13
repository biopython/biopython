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
from Bio.Restriction._Update.RestrictionCompiler import DictionaryBuilder


def standalone() :
    parser = optparse.OptionParser() 
    add = parser.add_option

    add('-i', '--install',
        action  = "store_true",
        dest    = 'i',
        default = False,
        help    = "compile and install the newly created file. "
        "default behaviour (without switch): "
        "Compile the enzymes and store them in the Updates folder")
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
    options, args = parser.parse_args()
    return options, args

if __name__ == '__main__' :
    options, args = standalone()
    Builder = DictionaryBuilder(options.rebase_password, options.ftp_proxy)
    Builder.build_dict()
    if options.i :
        Builder.install_dict()
    else :
        Builder.no_install()
    sys.exit()
