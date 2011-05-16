# Copyright 2001 by Katharine Lindner.  All rights reserved.
# Copyright 2006 by PeterC.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Hold GEO data in a straightforward format.

classes:
o Record - All of the information in an GEO record.

See http://www.ncbi.nlm.nih.gov/geo/
"""

class Record(object):
    """Hold GEO information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GEO data.

    Attributes:
    entity_type
    entity_id
    entity_attributes
    col_defs
    table_rows

    """
    def __init__(self):
        self.entity_type = ''
        self.entity_id = ''
        self.entity_attributes = {}
        self.col_defs = {}
        self.table_rows = []

    def __str__( self ):
        output = ''
        output = output + 'GEO Type: %s\n' % self.entity_type
        output = output + 'GEO Id: %s\n' % self.entity_id
        att_keys = self.entity_attributes.keys()
        att_keys.sort()
        for key in att_keys:
            contents = self.entity_attributes[ key ]
            if( type( contents ) == type( [] ) ):
                for item in contents:
                    try:
                        output = output + '%s: %s\n' % ( key, item[ :40 ] )
                        output = output + out_block( item[ 40: ] )
                    except:
                        pass
            elif( type( contents ) == type( '' ) ):
                output = output + '%s: %s\n' % ( key, contents[ :40 ] )
                output = output + out_block( contents[ 40: ] )
            else:
                print contents
                output = output + '%s: %s\n' % ( key, val[ :40 ] )
                output = output + out_block( val[ 40: ] )
        col_keys = self.col_defs.keys()
        col_keys.sort()
        output = output + 'Column Header Definitions\n'
        for key in col_keys:
            val = self.col_defs[ key ]
            output = output + '    %s: %s\n' % ( key, val[ :40 ] )
            output = output + out_block( val[ 40: ], '    ' )
        #May have to display VERY large tables,
        #so only show the first 20 lines of data
        MAX_ROWS = 20+1 # include header in count
        for row in self.table_rows[0:MAX_ROWS]:
            output = output + '%s: ' % self.table_rows.index( row )
            for col in row:
                output = output + '%s\t' % col
            output = output + '\n'
        if len(self.table_rows) > MAX_ROWS:
            output = output + '...\n'
            row = self.table_rows[-1]
            output = output + '%s: ' % self.table_rows.index( row )
            for col in row:
                output = output + '%s\t' % col
            output = output + '\n'
            
        return output

def out_block( text, prefix = '' ):
    output = ''
    for j in range( 0, len( text ), 80 ):
        output = output + '%s%s\n'  % ( prefix, text[ j: j + 80 ] )
    output = output + '\n'
    return output
