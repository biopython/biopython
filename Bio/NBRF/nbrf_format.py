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
# standard library
#http://www-nbrf.georgetown.edu/pirwww/pirhome.shtml
import string

# Martel
import Martel
from Martel import RecordReader
from Martel import Str
from Martel import AnyEol
from Martel import ToEol
from Martel import Group
from Martel import Alt
from Martel import Rep
from Martel import Rep1
from Martel import Any
from Martel import AnyBut


from Bio.NBRF.ValSeq import valid_sequence_dict



# --- first set up some helper constants and functions
# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


sequence_types = map( Str, valid_sequence_dict.keys() )
sequence_type = Group( "sequence_type", Alt( *sequence_types ) )
sequence_name = Group( "sequence_name", Rep1( Martel.Expression.Dot() ) )
name_line = Martel.Group( "name_line", \
    Str( ">" ) +
    sequence_type +
    Str( ";" ) +
    sequence_name +
    AnyEol() )
comment_line = Group( "comment_line", 
                      Rep1(AnyBut(' ')) +
                      ToEol( "comment" ) )

excluded_chars = chr( 0x2a ) + chr( 10 ) + chr( 13 )
sequence_text = Group( "sequence_text", \
    Martel.Rep1( AnyBut( excluded_chars ) ) )
sequence_final_text = Group( "sequence_final_text", \
    Martel.Rep1( AnyBut( excluded_chars ) ) )
sequence_final_line = Group( "sequence_final_line",
    sequence_final_text +
    Str( chr( 0x2a ) ) +
    AnyEol() )
sequence_line = Group( "sequence_line", sequence_text +
    AnyEol() )
sequence_block = Group( "sequence_block", Rep( sequence_line ) )
nbrf_record =  name_line + comment_line + sequence_block + sequence_final_line









