# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import string
from types import *
from Bio import File
from Bio import ParserSupport
from Bio import Fasta
from Bio import Alphabet

def title_to_ids(title):
    """Function to convert a title into the id, name, and description.

    This is just a quick-n-dirty implementation, and is definately not meant
    to handle every FASTA title line case.
    """
    # first split the id information from the description
    # the first item is the id info block, the rest is the description
    all_info = string.split(title, " ")
    id_info = all_info[0]
    rest = all_info[1:]
    descr = string.join(rest, " ")

    # now extract the ids from the id block
    # gi|5690369|gb|AF158246.1|AF158246
    id_info_items = string.split(id_info, "|")
    id = id_info_items[3] # the id with version info
    name = id_info_items[4] # the id without version info

    return id, name, descr

tests = [ 'lupine.nu', 'elderberry.nu', 'phlox.nu', 'centaurea.nu', \
    'wisteria.nu', 'sweetpea.nu', 'lavender.nu' ]
record_parser = Fasta.RecordParser()
sequence_parser = Fasta.SequenceParser(Alphabet.generic_dna, title_to_ids)

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Nucleic', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Nucleic', test )
    src_handle = open( datafile )
    data = sequence_parser.parse( src_handle )
    print data.id
    print data.name
    print data.description
    print data.seq

tests = [ 'aster.pro', 'rosemary.pro', 'rose.pro', 'loveliesbleeding.pro' ]
sequence_parser = Fasta.SequenceParser(Alphabet.generic_protein, title_to_ids)
for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Amino', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    print data

for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Amino', test )
    src_handle = open( datafile )
    data = sequence_parser.parse( src_handle )
    print data.id
    print data.name
    print data.description
    print data.seq[ :60]
