#!/usr/bin/env python
# Created: Wed Jun 21 10:26:53 2000
# Last changed: Time-stamp: <00/12/02 14:18:34 thomas>
# Thomas.Sicheritz@molbio.uu.se, http://evolution.bmc.uu.se/~thomas
# File: xbb_sequence.py

import sys
sys.path.insert(0, '.')
from Bio import Sequence


class xbb_sequence(Sequence):
    def __init__(self):
        ""



if __name__ == '__main__':
    test = xbb_sequence()
    
    
