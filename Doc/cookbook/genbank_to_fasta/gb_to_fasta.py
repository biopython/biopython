#!/usr/bin/env python
"""Convert a GenBank format file into a Fasta file, as simply as possible.

Usage:
    gb_to_fasta.py <file to convert>
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
from Bio import formats
from Bio.FormatIO import FormatIO

formatter = FormatIO("SeqRecord", formats["genbank"], formats["fasta"])
formatter.convert(input_handle, output_handle)

