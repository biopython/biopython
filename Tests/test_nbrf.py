"""Tests the basic functionality of the KEGG parsers.
"""

import os

import Bio.NBRF
from Bio.File import SGMLHandle

testfiles   = [  'B_nuc.pir', 'Cw_prot.pir', 'DMA_nuc.pir', 'DMB_prot.pir',
                 'clustalw.pir']

for file in testfiles:
    fh = open(os.path.join("NBRF", file))
    print "Testing Bio.NBRF on " + file + "\n\n"
    records = Bio.NBRF.Iterator(fh, Bio.NBRF.RecordParser(debug_level=0))
    while 1:
#    for j in range( 0, 1 ):
        record = records.next()
        if record is not None:
            print record
        else:
            break
    print "\n"





