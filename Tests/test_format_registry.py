# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys

import Bio


# Formats in Biopython:
#'blast', 'blastn', 'blastp', 'blastx', 'embl', 'embl/65', 'empty',
#'fasta', 'genbank', 'genbank-records', 'genbank-release',
#'ncbi-blastn', 'ncbi-blastp', 'ncbi-blastx', 'ncbi-tblastn',
#'ncbi-tblastx', 'search', 'sequence', 'swissprot', 'swissprot/38',
#'swissprot/40', 'tblastn', 'tblastx', 'wu-blastn', 'wu-blastp',
#'wu-blastx'


def print_title(name):
    linelen = (60 - len(name)+2)/2
    line = "=" * linelen
    print "%s %s %s" % (line, name, line)

def _open(filename, *args):
    filename = os.path.join('Registry', filename)
    return open(filename, *args)


print_title("Testing Format Identifier")
x = Bio.formats['swissprot/38'].identifyFile(_open('EDD_RAT.dat'))
print x.name   # swissprot/38
x = Bio.formats['swissprot'].identifyFile(_open('EDD_RAT.dat'))
print x.name   # swissprot/38
x = Bio.formats['sequence'].identifyFile(_open('EDD_RAT.dat'))
print x.name   # swissprot/38

x = Bio.formats['search'].identifyFile(_open('bt001'))
print x.name   # ncbi-blastp


print_title("Testing Format Parser")
from Bio import SeqRecord
format = Bio.formats['sequence'].identifyFile(_open('EDD_RAT.dat'))
parser = format.make_parser()
builder = Bio.formats.find_builder(Bio.formats['sequence'], SeqRecord.io)
parser.setContentHandler(builder)
parser.parseFile(_open('EDD_RAT.dat'))
print builder.document.id    # EDD_RAT
print builder.document.seq   # Seq('ARR...')

from Bio import Search
r = Search.io.readFile(_open("bt001"))
print r.algorithm.name       # BLASTP


print_title("Testing FormatIO")
r = SeqRecord.io.readFile(_open("EDD_RAT.dat"))
print r.__class__
for seq in r:
    print seq.id              # EDD_RAT
    print seq.description     # Ubiquitin ...
    print seq.seq             # Seq('ARR...')
    
SeqRecord.io.convert(infile=_open("EDD_RAT.dat"))   # FASTA-format

print_title("Testing Reading Multiple Records")
SeqRecord.io.convert(infile=_open("seqs.fasta"))    # 2 FASTA-format

