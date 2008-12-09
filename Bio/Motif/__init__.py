# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Module containing different tools for sequence motif analysis.

it contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.
It also inlcudes functionality for parsing AlignACE and MEME programs
"""
from __future__ import generators
from Motif import Motif
from AlignAceParser import AlignAceParser, CompareAceParser
from MEMEParser import MEMEParser,MASTParser
from Thresholds import score_distribution
from AlignAceStandalone import AlignAce

_parsers={"AlignAce":AlignAceParser,
          "MEME":MEMEParser,
          }

def _from_pfm(handle):
    return Motif().from_jaspar_pfm(handle)

def _from_sites(handle):
    return Motif().from_jaspar_sites(handle)

_readers={"pfm": _from_pfm,
          "sites": _from_sites
          }

          
def parse(handle,format):
    """Parses an output file of motif finding programs.

    Currently supported formats:
    -AlignAce
    -MEME

    you can also use single-motif formats:
    -pfm
    -sites
    """
    try:
        parser=_parsers[format]
        
    except KeyError:
        try: #not a true parser, try reader formats
            reader=_readers[format]
        except:
            raise ValueError("Wrong parser format")
        else: #we have a proper reader 
            yield reader(handle)
    else: # we have a proper reader
        for m in parser().parse(handle).motifs:
            yield m

def read(handle,format):
    """Reads a motif from a handle using a specified file-format.

    Currently supported formats:
    -sites
    -pfm

    yau can also use parser formats to get the first motif
    """
    return parse(handle,format).next()


        
        
    
