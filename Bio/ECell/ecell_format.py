# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read ECell formatted files.

This is a huge regular regular expression for ECell, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://www.bioinformatics.org/ecell2/
Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""
# standard library
#http://www.bioinformatics.org/ecell2/
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
from Martel import Expression





# --- first set up some helper constants and functions
# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


excluded_chars = ' ' + chr( 0x09 ) + chr( 10 ) + chr( 13 )

block_type = Group( "block_type", Expression.NoCase( Str( "Type" ) ) )
header_line = Group( "header_line", \
    block_type + ToEol())
tab = Group( "tab", Str( '\t' ) )
system_tag = Group( "system_tag", Expression.NoCase( Str( "system" ) ) )
reactor_tag = Group( "reactor_tag", Expression.NoCase( Str( "Reactor" ) ) )
substance_tag = Group( "substance_tag", Expression.NoCase( Str( "Substance" ) ) )
system_line = Group( "system_line", system_tag + ToEol() )
reactor_line = Group( "reactor_line", reactor_tag + ToEol() )
substance_line = Group( "substance_line", substance_tag + ToEol() )
continuation_line = Group( "continuation_line", tab + ToEol() )
include_line = Group( "include_line", Str( 'include' ) + ToEol())

substance_multiline = Group( "substance_multiline", \
    substance_line +
    Rep( continuation_line ) )

reactor_multiline = Group( "reactor_multiline", \
    reactor_line +
    Rep( continuation_line ) )

system_block = Group( "system_block", \
               Rep1( system_line ) )
reactor_block = Group( "reactor_block", \
               Rep1( reactor_multiline ) )
substance_block = Group( "substance_block", \
               Rep1( substance_multiline ) )
valid_block = Group( "valid_block",
    header_line +
    Alt( system_block, reactor_block, substance_block ) )
valid_contents = Group( "valid_contents", Rep1( valid_block ) )
ecell_record =  valid_contents









