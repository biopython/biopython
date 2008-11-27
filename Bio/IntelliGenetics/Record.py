# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read IntelliGenetics formatted files (DEPRECATED).

This module defines a Record class to hold a sequence from the IntelliGenetics/
MASE file format in a similar represenation to the original raw file.
"""


from Bio.Seq import Seq
"""Hold IntelliGenetics data in a straightforward format.

classes:
o Record - All of the information in an IntelliGenetics record.
"""

class Record:
    """Hold IntelliGenetics information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GenBank data.

    Attributes:
    comments
    title
    sequence
    """
    def __init__(self):
        self.comments = []
        self.title = ''
        self.sequence = Seq('')

    def __str__( self ):
        output =  'Title: %s\n' % self.title
        for comment in self.comments:
            output = output + '%s\n' % comment
        output = output + out_sequence( self.sequence.data )
        return output

def out_sequence( seq ):
    output = ''
    for j in range( 0, len( seq ), 80 ):
        output = output + '%s\n'  % seq[ j: j + 80 ]
    output = output + '\n'
    return output










