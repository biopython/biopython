"""Tests the basic functionality of the KEGG parsers.
"""

import os

import Bio.NBRF
from Bio.RecordFile import RecordFile
from Bio.File import SGMLHandle

testfiles   = [  'B_nuc.pir', 'Cw_prot.pir', 'DMA_nuc.pir', 'DMB_prot.pir'
                            ]

for file in testfiles:
    fh = open(os.path.join("NBRF", file))
    record_file = RecordFile( fh, '>', chr( 0x2a )  )
    print "Testing Bio.NBRF on " + file + "\n\n"
    records = Bio.NBRF.Iterator(record_file, Bio.NBRF.RecordParser(debug_level=0))
#    records = Bio.NBRF.Iterator(record_file, None)
    while 1:
#    for j in range( 0, 1 ):
        record = records.next()
        if record is not None:
            print record
        else:
            break
    print "\n"





