# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read CDD formatted files.

This is a huge regular regular expression for CDD, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

# standard library
http://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml
"""
import string

# Martel
import Martel
from Martel import RecordReader
from Martel import Str
from Martel import AnyEol
from Martel import ToEol
from Martel import Group
from Martel import Alt
from Martel import Opt
from Martel import Rep
from Martel import Rep1
from Martel import Any
from Martel import AnyBut
from Martel import Assert
from Martel import AssertNot





# --- first set up some helper constants and functions
# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

upper_alpha = Any( "ABCDEFGHIJKLMNOPQRSTUVWXYZ" )
white_space = Any( "\t " )
eols = chr( 13 ) + chr( 10 )
white_spaces = Rep( white_space )
summary_line = Str( "CD summary" ) + ToEol()

cd_tag = Group( "cd_tag", Str( "CD:" ) )
description_tag = Group( "description_tag", Str( "Description:" ) )
status_tag = Group( "status_tag", Str( "CD status:" ) )
source_tag = Group( "source_tag", Str( "Source:" ) )
date_tag = Group( "date_tag", Str( "Created:" ) )
reference_tag = Group( "reference_tag", Str( "References:" ) )
taxonomy_tag =  Group( "taxonomy_tag", Str( "Taxonomy spanned:" ) )
aligned_tag = Group( "aligned_tag", Str( "Aligned sequences:" ) )
representative_tag = Group( "representative_tag", Str( "Representative:" ) )
range_tag = Group( "range_tag", Str( "Aligned range:" ) )
sequence_tag = Group( "sequence_tag", Str( "Sequence:" ) )
has_tag = Alt( cd_tag, description_tag, status_tag, source_tag, date_tag, \
    reference_tag, taxonomy_tag, aligned_tag, representative_tag, range_tag, sequence_tag )

cd_key_line = cd_tag + white_spaces + AnyEol()
description_key_line = description_tag + white_spaces + AnyEol()
status_key_line = status_tag + white_spaces + AnyEol()
source_key_line = source_tag + white_spaces + AnyEol()
date_key_line = date_tag + white_spaces + AnyEol()
reference_key_line = reference_tag + white_spaces + AnyEol()
taxonomy_key_line = taxonomy_tag + white_spaces + AnyEol()
aligned_key_line = aligned_tag + white_spaces + AnyEol()
representative_key_line = representative_tag + white_spaces + AnyEol()
range_key_line = range_tag + white_spaces + AnyEol()
sequence_key_line = sequence_tag + white_spaces + AnyEol()

cd_contents_line = Group( "cd_contents_line", AssertNot( has_tag ) + ToEol() )
description_contents_line = AssertNot( has_tag ) + ToEol()
status_contents_line = AssertNot( has_tag ) + ToEol()
source_contents_line = AssertNot( has_tag ) + ToEol()
date_contents_line = AssertNot( has_tag ) + ToEol()
reference_contents_line = AssertNot( has_tag ) + ToEol()
taxonomy_contents_line = AssertNot( has_tag ) + ToEol()
aligned_contents_line = AssertNot( has_tag ) + ToEol()
representative_contents_line = AssertNot( has_tag ) + ToEol()
range_contents_line = AssertNot( has_tag ) + ToEol()
sequence_contents_line = Group( "sequence_contents_line", \
                            white_spaces + Rep1( upper_alpha ) + white_spaces + AnyEol() )
sentinel_line = white_spaces + Str( "Definition" ) + white_spaces + AnyEol()
boiler_plate = AssertNot( sentinel_line ) + ToEol()
definition_line = Group( "definition_line", \
    Rep( AnyBut( eols + '[' ) ) + Str( '[CD]' ) + white_spaces + AnyEol() )
pdb_id_line = AssertNot( definition_line ) + ToEol()
pdb_id_multiline = Group( "pdb_id_multiline", Rep1( pdb_id_line ) )
table_entry = Group( "table_entry", \
    pdb_id_multiline + definition_line )
table = Group( "table", Rep1( table_entry ) )

cd_contents_multiline = Group( "cd_contents_multiline", \
    Rep( cd_contents_line ) )
description_contents_multiline = Group( "description_contents_multiline", \
    Rep( description_contents_line ) )
status_contents_multiline = Group( "status_contents_multiline", \
    Rep( status_contents_line ) )
source_contents_multiline = Group( "source_contents_multiline", \
    Rep( source_contents_line ) )
date_contents_multiline = Group( "date_contents_multiline", \
    Rep( date_contents_line ) )
reference_contents_multiline = Group( "reference_contents_multiline", \
    Rep( reference_contents_line ) )
taxonomy_contents_multiline = Group( "taxonomy_contents_multiline", \
    Rep( taxonomy_contents_line ) )
aligned_contents_multiline = Group( "aligned_contents_multiline", \
    Rep( aligned_contents_line ) )
representative_contents_multiline = Group( "representative_contents_multiline", \
    Rep( representative_contents_line ) )
range_contents_multiline = Group( "range_contents_multiline", \
    Rep( range_contents_line ) )
sequence_contents_multiline = Group( "sequence_contents_multiline", \
    Rep( sequence_contents_line ) )

cd_block = cd_key_line + cd_contents_multiline
description_block = description_key_line + description_contents_multiline
status_block = status_key_line + status_contents_multiline
source_block = source_key_line + source_contents_multiline
date_block = date_key_line + date_contents_multiline
reference_block = Assert(reference_tag ) + reference_key_line + \
    reference_contents_multiline
taxonomy_block = taxonomy_key_line + taxonomy_contents_multiline
aligned_block = aligned_key_line + aligned_contents_multiline
representative_block = representative_key_line + representative_contents_multiline
range_block = range_key_line + range_contents_multiline
sequence_block = sequence_key_line + sequence_contents_multiline
trailer_line = ToEol()

cdd_record = summary_line + cd_block + description_block + status_block + \
    source_block + date_block + Opt( reference_block ) + taxonomy_block + \
    aligned_block + representative_block + range_block + sequence_block + \
    Rep( boiler_plate ) + sentinel_line + table









