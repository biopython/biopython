# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read SAF formatted files.

This is a huge regular regular expression for SAF, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://www.embl-heidelberg.de/predictprotein/Dexa/optin_safDes.html


Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""
#http://www.embl-heidelberg.de/predictprotein/Dexa/optin_safDes.html


# Martel
import Martel
from Martel import Str
from Martel import AnyEol
from Martel import ToEol
from Martel import Group
from Martel import Alt
from Martel import Rep
from Martel import Rep1
from Martel import Any
from Martel import Opt
from Martel import ToSep
from Martel.Expression import Assert



# --- first set up some helper constants and functions
# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


digits = "0123456789"
valid_sequence_characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-. \t'
white_space = "\t "
valid_residue_characters = digits + white_space + chr( 0x2e )
residue_number_line = Group( "residue_number_line", \
                      Rep1( Any( valid_residue_characters ) ) +
                      AnyEol())
comment_line = Group( "comment_line", \
               Str( "#" ) +
               ToEol() )
ignored_line = Group( "ignored_line", \
               Alt( comment_line, residue_number_line ) )
candidate_line = Group( "candidate_line", \
                 Assert( Str( "#" ), 1 ) +
                 Assert( Any( valid_residue_characters ), 1 ) +
                 ToSep( sep = ' ' ) +
                 Rep( Any( valid_sequence_characters ) ) +
                 ToEol() )
saf_record =  Group( "saf_record", \
    candidate_line + Rep( Alt( candidate_line, ignored_line ) ) + Opt( Str( "#" ) ) )









