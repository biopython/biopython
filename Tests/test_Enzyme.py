# Copyright 1999 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
import Bio.File
from Bio import ParserSupport
from Bio import Enzyme


tests = [ 'lipoprotein.txt', 'proline.txt', 'valine.txt', 'lactate.txt' ]

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Enzymes', test )
    enzyme_handle = open( datafile )
    enzyme_undo_handle = Bio.File.UndoHandle( enzyme_handle )
    enzyme_consumer = ParserSupport.TaggingConsumer()
    enzyme_scanner = Enzyme._Scanner()
    enzyme_scanner.feed( enzyme_undo_handle, enzyme_consumer )

