# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read IntelliGenetics formatted files.

This is a huge regular regular expression for IntelliGenetics, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://hiv-web.lanl.gov/ALIGN_97/HIV12SIV-index.html

Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""
# standard library
import string


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










