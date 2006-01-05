# Copyright 2002 by Katharine Lindner.  All rights reserved.
# Copyright 2005 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read GEO formatted files.

This is a huge regular regular expression for GEO.

http://www.ncbi.nlm.nih.gov/geo/

The file format is described here, last revised in 2005:

http://www.ncbi.nlm.nih.gov/projects/geo/info/soft2.html#SOFTformat

There are four basic types of lines in GEO SOFT files,
those that start with ^, ! and # and the fourth type
are the tab separated data lines (including the column header line)

This parser considers any ^ line to be the start of a new GEO record.
Each record look like this:

The ^ line
Optionally, many ! lines
Optionally, many # lines
Optionally, a data table

The data table, used to look like this:

ID_REF  VALUE   ABS_CALL        DETECTION P-VALUE
141200_at       36.6    A       0.818657
141201_at       41.5    A       0.703191
...
141219_at       223.5   P       0.007827

As of 2005, the NCBI added two extra ! lines, one before and one after:

!Sample_table_begin
ID_REF  VALUE   ABS_CALL        DETECTION P-VALUE
141200_at       36.6    A       0.818657
141201_at       41.5    A       0.703191
...
141219_at       223.5   P       0.007827
!Sample_table_end

The exact text of these markers will vary depending on the table type,
but it seems that it should match the ^entity type at the start of the
record.

Except for the platform files, where they do this instead for the table:

ID (tab) ...
!table_begin
AFFX-BioB-5_at (tab) ...
AFFX-BioC-3_at (tab) ...
...
SYNPBR322_tet_w_at (tab) ...
!table_end

Which is just plain awkward of them, but not the end of the world.
"""

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

#There have been a few new "ENTITIES" added to the file format in the last few years...
valid_entity = NoCase( Alt( Str("PLATFORM"), Str("SAMPLE"), Str("SERIES"), 
                            Str("DATABASE"), Str("DATASET"), Str("SUBSET"),
                            Str("ANNOTATION")))


#Calling lines starting ^ attribute lines
entity_line = Group( "entity_line", \
               Str( "^" ) +
               valid_entity +
               ToEol() )

#Calling lines starting ! entity lines
attribute_line = Group( "attribute_line", \
               Str( "!" ) +
               ToEol() )

#Calling lines starting # column heading lines
col_heading_line = Group( "col_heading_line", \
               Str( "#" ) +
               ToEol() )

#Calling the data rows (which don't start ^,! or #) row_line
row_line = Group( "row_line", \
                 AnyBut( "^!#" ) +
                 ToEol() )

table_begin_line = Group("table_begin", \
                          Str( "!" ) + valid_entity + Str("_table_begin") +
                          AnyEol() )
table_end_line = Group("table_end", \
                          Str( "!" ) + valid_entity + Str("_table_end") +
                          AnyEol() )

ann_table_begin_line = Group("ann_table_begin", \
                          Str( "!table_begin") +
                          AnyEol() )

ann_table_end_line = Group("ann_table_begin", \
                          Str( "!table_end") +
                          AnyEol() )


geo_record =  Group( "geo_record",
                     #Must have the ^ entity line:
                     entity_line +
                     #Can then have none or more ! attribute lines:
                     Rep(attribute_line) +
                     #Can then have none or more # column headings:
                     Rep(col_heading_line) +
                     #Can then have a table, in one of three forms:
                     # (a) New annotation table with nasty !table_begin placement
                     #     after the header row and before the genes.
                     # (b) New data table with !XXX_table_begin before both the
                     #     header row and the genes
                     # (c) Old style data table with no table begin/end lines
                     Alt (row_line + ann_table_begin_line + Rep(row_line) + ann_table_end_line,
                          table_begin_line + Rep(row_line) + table_end_line,
                          Rep(row_line)) +
                     #Finally, allow none or more blank lines (just in case the
                     #file has been edited by hand):
                     Rep(AnyEol()))
