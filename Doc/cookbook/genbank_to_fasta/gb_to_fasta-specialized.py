#!/usr/bin/env python
"""Convert GenBank into Fasta, with specialized control over titles.

Usage:
    gb_to_fasta-specialized.py <file to convert>
"""
import sys
import os

# get the file names
input_file = sys.argv[1]
basename, ext = os.path.splitext(input_file)
output_file = basename + ".fasta"
input_handle = open(input_file)
output_handle = open(output_file, "w")

# do the conversion
from Bio import GenBank
from Bio import Fasta

iterator = GenBank.Iterator(input_handle, GenBank.RecordParser())
for gb_rec in iterator:
    fasta_rec = Fasta.Record()
    fasta_rec.title = "%s|%s|%s %s" % \
     (gb_rec.gi, gb_rec.locus, gb_rec.version, gb_rec.definition)
    fasta_rec.sequence = gb_rec.sequence
    output_handle.write(str(fasta_rec) + "\n")
