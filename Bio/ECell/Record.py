# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read ECell formatted files.

This is a huge regular regular expression for Ecell, built using
the 'regular expressiona on steroids' capabilities of Martel.

#http://www.bioinformatics.org/ecell2/



Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""
# standard library
import string





"""Hold ECell data in a straightforward format.

classes:
o Record - All of the information in an ECell record.
"""

class Record:
    """Hold ECell information in a format similar to the original record.


    """
    def __init__(self):
        self.cell_dict = {}
        self._max_dict = {}
        self._duplicates = []
        self._ecell_version = '0.9.4.6'
        self._name = 'Ecell module'
        self._version = ''
        self.include_buf = ''
        self.num_systems = 0
        self.num_substances = 0
        self.num_reactors = 0
        self.include_block = ''
        self.contains_concentration = 0


    def __str__( self ):
        output = '# this is ecell rule file for E-CELL' + self._ecell_version + '\n'
        output = output + '# converted by %s %s\n\n' % ( self._name, self._version )
        if( self.contains_concentration ):
            output = output + 'include(qty.er)\n'
        output = output + self.include_buf
        system_output = self._print_systems()
        substance_output = self._print_substances()
        reactor_output = self._print_reactors( substance_output )
        output = output + system_output + reactor_output
        return output

    def _print_systems( self ):
        output = ''
        for system in range( 1, self.num_systems + 1 ):
            composite_key = 'system' + str( system ) + 'class0'
            output = output + '\nsystem %s' % self.cell_dict[ composite_key ]
            composite_key = 'system' + str( system ) + 'path0'
            output = output + '(%s:' % self.cell_dict[ composite_key ]
            composite_key = 'system' + str( system ) + 'id0'
            output = output + '%s,' % self.cell_dict[ composite_key ]
            composite_key = 'system' + str( system ) + 'name0'
            output = output +  '"%s")\n' % self.cell_dict[ composite_key ]
            output = output + '{\n'

            output = output + '\tStepper SlaveStepper;\n'
            composite_key = 'system' + str( system ) + 'inside0'
            if( self.cell_dict.has_key( composite_key ) ):
                output = output + '\tInside  %s;\n' % self.cell_dict[ composite_key ]
            composite_key = 'system' + str( system ) + 'outside0'
            if( self.cell_dict.has_key( composite_key ) ):
                output = output + '\tOutside %s;\n' % self.cell_dict[ composite_key ]
            composite_key = 'system' + str( system ) + 'volumeindex0'
            if( self.cell_dict.has_key( composite_key ) ):
                output = output + '\tVolumeIndex %s;\n' % self.cell_dict[ composite_key ]
            output = output + '}\n'
        return output

    def _print_substances( self ):
        output = ''
        for substance in range( 1, self.num_substances + 1 ):
            composite_key = 'substance' + str( substance ) + 'class0'
            if( self.cell_dict.has_key( composite_key ) ):
                output = output + 'substance %s' % self.cell_dict[ composite_key ]
                composite_key = 'substance' + str( substance ) + 'path0'
                output = output + '(%s:' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'id0'
                output = output + '%s,' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'name0'
                output = output + '"%s")\n' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'qty0'
                output = output + '{\n\tQuantity %s;\n}\n' % get_entry( self.cell_dict, composite_key )
            else:
                composite_key = 'substance' + str( substance ) + 'path0'
                output = output + 'substance %s:' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'id0'
                output = output + '%s ' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'name0'
                output = output + '"%s" ' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'qty0'
                output = output + '%s;\n' % get_entry( self.cell_dict, composite_key )
            qty_key = 'substance' + str( substance ) + 'qty1'
            conc_key = 'substance' + str( substance ) + 'conc1'
            if( get_entry( self.cell_dict, qty_key ).lower() == 'fix' ):
                composite_key = 'substance' + str( substance ) + 'path0'
                output = output + 'fix %s:' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'id0'
                output = output + '%s;\n' % get_entry( self.cell_dict, composite_key )
            elif( get_entry( self.cell_dict, conc_key ).lower() == 'fix' ):
                composite_key = 'substance' + str( substance ) + 'path0'
                output = output + 'fix %s:' % get_entry( self.cell_dict, composite_key )
                composite_key = 'substance' + str( substance ) + 'id0'
                output = output + '%s;\n' % get_entry( self.cell_dict, composite_key )
        return output

    def _print_reactors( self, output ):
        volume_buf = '\n'
        self._check_duplicates()
        for reactor in range( 1, self.num_reactors + 1 ):
            output = output + '\nreactor '
            prefix = 'reactor' + str( reactor )
            composite_key = prefix + 'class0'
            output = output + self.cell_dict[ composite_key ]
            composite_key = prefix + 'path0'
            output = output + '(%s:' % self.cell_dict[ composite_key ]
            composite_key = prefix + 'id0'
            output = output + '%s,' % self.cell_dict[ composite_key ]
            composite_key = prefix + 'name0'
            output = output + '"%s")\n' % self.cell_dict[ composite_key ]

            composite_key = 's_' + str( reactor )
            num_substrates = get_entry( self._max_dict, composite_key )
            output = output + '{\n'
            for substrate in range( 1, num_substrates + 1 ):
                output = output + '\tSubstrate '
                composite_key = prefix + 's_path' + str( substrate )
                output = output + get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 's_id' + str( substrate )
                output = output + ':%s ' % get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 's_coeff' + str( substrate )
                output = output + '%s' % get_entry( self.cell_dict, composite_key )
                output = output + ';\n'

            composite_key = 'p_' + str( reactor )
            num_products = get_entry( self._max_dict, composite_key )
            for product in range( 1, num_products + 1 ):
                output = output + '\tProduct '
                composite_key = prefix + 'p_path' + str( product )
                output = output + get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'p_id' + str( product )
                output = output + ':%s ' % get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'p_coeff' + str( product )
                output = output + '%s' % get_entry( self.cell_dict, composite_key )
                output = output + ';\n'

            composite_key = 'c_' + str( reactor )
            num_catalysts = get_entry( self._max_dict, composite_key )
            for catalyst in range( 1, num_catalysts + 1 ):
                output = output + '\tCatalyst '
                composite_key = prefix + 'c_path' + str( catalyst )
                output = output + get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'c_id' + str( catalyst )
                output = output + ':%s' % get_entry( self.cell_dict, composite_key )
                output = output + ';\n'

            composite_key = 'e_' + str( reactor )
            num_effectors = get_entry( self._max_dict, composite_key )
            for effector in range( 1, num_effectors + 1 ):
                output = output + '\tEffector '
                composite_key = prefix + 'e_path' + str( effector )
                output = output + get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'e_id' + str( effector )
                output = output + ':%s ' % get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'e_coeff' + str( effector )
                output = output + '%s;\n' % get_entry( self.cell_dict, composite_key )
                output = output + ';\n'

            composite_key = 'o_' + str( reactor )
            num_options = get_entry( self._max_dict, composite_key )
            for option in range( 1, num_options + 1 ):
                composite_key = prefix + 'o_type' + str( option )
                output = output + '\t%s ' % get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'o_path' + str( option )
                output = output + get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'o_id' + str( option )
                output = output + ':%s ' % get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'o_coeff' + str( option )
                output = output + '%s;\n' % get_entry( self.cell_dict, composite_key )
                output = output + ';\n'

            composite_key = 'arg_tag' + str( reactor )
            num_args = get_entry( self._max_dict, composite_key )
            for arg in range( 1, num_args + 1 ):
                composite_key = prefix + 'arg_tag' + str( arg )
                output = output + '\t%s ' % get_entry( self.cell_dict, composite_key )
                composite_key = prefix + 'arg_coeff' + str( arg )
                output = output + '%s;\n' % get_entry( self.cell_dict, composite_key )

            for system in range( 1, self.num_systems + 1 ):
                path_key = prefix + 'path0'
                id_key = prefix + 'id0'
                reactor_path = get_entry( self.cell_dict, path_key )
                reactor_id = get_entry( self.cell_dict, id_key )
                path_id = '%s:%s' % ( reactor_path, reactor_id )
                volume_key = 'system' + str( system ) + 'volumeindex0'
                volume_index = get_entry( self.cell_dict, volume_key )
                if( path_id == volume_index ):
                    output = output + '\tInitialActivity '
                    init_act_key = prefix + 'init_act0'
                    init_act0 = get_entry( self.cell_dict, init_act_key )
                    if( init_act0 == '' ):
                        init_act0 = get_entry( self.cell_dict, prefix + 'init_act1' )
                    output = output + '%s;\n' % init_act0
                    if( not ( system in self._duplicates )  ):
                        volume_buf = volume_buf + R'_SETVOLUME('
                        volume_buf = volume_buf + '%s,%s)\n' % ( reactor_path, init_act0 )
                        volume_buf = volume_buf.replace( r'//', '\/' )

            output = output + '}\n'

        if( self.contains_concentration ):
            output = volume_buf + output
        return output

    def _check_duplicates( self ):
        self._duplicates = []
        for system in range( 1, self.num_systems + 1 ):
            target_path_key = 'system' + str( system ) + 'path0'
            target_id_key = 'system' + str( system ) + 'id0'
            for other in range( 1, system ):
                match_path_key = 'system' + str( other ) + 'path0'
                match_id_key = 'system' + str( other ) + 'id0'
                if( self.cell_dict.has_key( target_path_key ) and \
                    self.cell_dict.has_key( target_id_key ) ):
                        if self._match_cell_dict_entry( target_path_key, match_path_key ):
                            if self._match_cell_dict_entry( target_id_key, match_id_key ):
                                self._duplicates.append( system )

    def _match_cell_dict_entry( self, key, other_key ):
        if( get_entry( self.cell_dict, key ) == get_entry( self.cell_dict, other_key ) ):
            return 1
        return 0


def get_entry( dict, key ):
    try:
        entry = dict[ key ]
    except KeyError:
        entry = ''
    return entry









