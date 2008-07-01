# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read SAF formatted files.

This is a huge regular regular expression for Saf, built using
the 'regular expressiona on steroids' capabilities of Martel.

#http://www.embl-heidelberg.de/predictprotein/Dexa/optin_safDes.html


Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""

from Bio.Align.Generic import Alignment
import Bio.Alphabet



"""Hold SAF data in a straightforward format.

classes:
o Record - All of the information in an Saf record.
"""

class Record:
    """Hold Saf information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at Saf data.

    Attributes:
    alignment

    """
    def __init__(self):
        self.alignment = Alignment( Bio.Alphabet.generic_alphabet )

    def __str__( self ):
        output = ''
        sequences = self.alignment.get_all_seqs()
        for sequence_record in sequences:
            output = output + '%s\n' % sequence_record.description
            output = output + out_sequence( sequence_record.seq.data )
        return output

def out_sequence( seq ):
    output = ''
    for j in range( 0, len( seq ), 80 ):
        output = output + '%s\n'  % seq[ j: j + 80 ]
    output = output + '\n'
    return output










