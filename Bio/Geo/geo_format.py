# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read GEO formatted files.

This is a huge regular regular expression for GEO.

#http://www.ncbi.nlm.nih.gov/geo/

Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""
# standard library
#http://www.ncbi.nlm.nih.gov/geo/
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
from Martel import RepN
from Martel import Opt
from Martel import ToSep
from Martel.Expression import Assert
from Martel.Expression import NoCase



# --- first set up some helper constants and functions
# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


digits = "0123456789"
valid_sequence_characters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ-. \t'
valid_entity = NoCase( Alt( Str( "PLATFORM" ), Str( "SAMPLE" ), Str( "SERIES" ) ) )
white_space = "\t "
entity_line = Group( "entity_line", \
               Str( "^" ) +
               valid_entity +
               ToEol() )
attribute_line = Group( "attribute_line", \
               Str( "!" ) +
               ToEol() )
col_heading_line = Group( "col_heading_line", \
               Str( "#" ) +
               ToEol() )
row_line = Group( "row_line", \
                 Assert( Str( "^" ), 1 ) +
                 Assert( Str( "!" ), 1 ) +
                 Assert( Str( "#" ), 1 ) +
                 ToEol() )
entity_block = Group( "entity_block", \
                       entity_line + \
                       Rep( attribute_line ) )
col_heading_block = Group( "col_heading_block", Rep( col_heading_line ) )
row_block = Group( "row_block", Rep( row_line ) )
geo_record =  Group( "geo_record", \
    entity_block + Opt( col_heading_block + row_block ) )









