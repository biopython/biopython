# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read NBRF formatted files.

This is a huge regular regular expression for NBRF, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://www-nbrf.georgetown.edu/pirwww/pirhome.shtml
"""

# Martel
import Martel
from Martel import Str
from Martel import AnyEol, UntilEol
from Martel import Group
from Martel import Alt
from Martel import Rep
from Martel import Rep1
from Martel import AnyBut
from Martel import UntilSep

from Bio.NBRF.ValSeq import valid_sequence_dict

sequence_types = map( Str, valid_sequence_dict.keys() )
sequence_type = Group( "sequence_type", Alt( *sequence_types ) )
name_line = Martel.Group( "name_line", \
    Str( ">" ) +
    sequence_type +
    Str( ";" ) +
    UntilEol("sequence_name") +
    AnyEol() )

comment_line = UntilEol("comment") + AnyEol()

# 0x2a -- '*'
# 10 -- '\n', 13 -- '\r' newline endings
excluded_chars = chr(0x2a) + chr(10) + chr(13)
# sequence lines with only sequence
sequence_text = Group( "sequence_text", \
    Martel.Rep1( AnyBut( excluded_chars ) ) )
sequence_line = Group( "sequence_line", sequence_text +
    AnyEol())
# the final line, has a '*' and potentially some sequence
sequence_final_line = Group( "sequence_final_line",
        UntilSep("sequence_final_text", chr(0x2a)) + Str(chr(0x2a)) +
        Rep1(AnyEol()))

sequence_block = Group("sequence_block", Rep( sequence_line ))
nbrf_record =  name_line + comment_line + sequence_block + sequence_final_line
