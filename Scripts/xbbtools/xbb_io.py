#!/usr/bin/env python
# Created: Wed Jun 21 13:46:35 2000
# Last changed: Time-stamp: <00/12/02 14:18:23 thomas>
# Thomas.Sicheritz@molbio.uu.se, http://evolution.bmc.uu.se/~thomas
# File: xbb_io.py

import os, sys  # os.system, sys.argv
sys.path.insert(0, '.')
sys.path.insert(0, os.path.expanduser('~thomas/cbs/python/biopython'))

from Bio.ParserSupport import *
from Bio import Fasta

class xbb_io:
    def __init__(self):
        ""

    def error(self, str):
        print str
        
    def read_fasta_file(self, file):
        genes = []
        iter = Fasta.Iterator(handle = open(file), parser = Fasta.RecordParser())
        while 1:
            rec = iter.next()
            if not rec: break
            genes.append((rec.sequence, rec.title))

        return genes
    
            


