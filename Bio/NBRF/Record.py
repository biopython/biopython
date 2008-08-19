# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read NBRF formatted files.

This is a huge regular regular expression for NBRF, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://www-nbrf.georgetown.edu/pirwww/pirhome.shtml


Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""


from Bio.Seq import Seq
from Bio.NBRF.ValSeq import valid_sequence_dict



"""Hold NBRF data in a straightforward format.

classes:
o Record - All of the information in an NBRF record.
"""

class Record:
    """Hold NBRF information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at NBRF data.

    Attributes:
    sequence_type
    sequence_name
    comment
    sequence

    """
    def __init__(self):
        self.sequence_type = ''
        self.sequence_name = ''
        self.comment = ''
        self.sequence = Seq('')

    def __str__( self ):
        sequence_type = valid_sequence_dict[ self.sequence_type ]
        output =  'Sequence type %s\n' % sequence_type
        output =  output + 'Sequence name %s\n' % self.sequence_name
        output = output + '%s\n' % self.comment
        output = output + out_sequence( self.sequence.data )
        return output

def out_sequence( seq ):
    output = ''
    for j in range( 0, len( seq ), 80 ):
        output = output + '%s\n'  % seq[ j: j + 80 ]
    output = output + '\n'
    return output










