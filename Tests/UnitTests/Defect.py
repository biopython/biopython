# Defect
#
#
# MODIFICATION NOTES: See bottom of file.

# Copyright (c) 1999 Katharine Lindner
# This module is free software; you can redistribute it and/or modify
# it under the same terms as Python itself.

# IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
# SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF
# THIS CODE, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
#
# THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE.  THE CODE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
# AND THERE IS NO OBLIGATION WHATSOEVER TO PROVIDE MAINTENANCE,
# SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

import sys
import traceback
class Defect:




    #default for 0 only for DEGUGGIN!!!
    def __init__( self, message, trace = 0 ):
        self.message = message
        self.trace = trace

    def print_message( self, outfile = sys.stdout ):
        """
        Title     print_message
        Usage     print_message
                  print_message( file_handle )
        Function  Print a failure message

        Returns
        Argument  File handle
        """


	outfile.write( '%s\n' % ( self.message ) )

    # Title     : print_trace
    # Usage     : print_trace
    #           : print_trace( file_handle )
    # Function  : Print a failure trace
    #           :
    # Returns   :
    # Argument  : File handle

    def print_trace( self, outfile = sys.stdout ):
        trace = self.trace
        frame = trace.tb_frame
	code = frame.f_code
	outfile.write( '%s %s %d\n' % ( code.co_filename, code.co_name, code.co_firstlineno ) )
        if frame.f_back:
            frame = frame.f_back
            code = frame.f_code
            outfile.write( '%s %s %d\n' % ( code.co_filename, code.co_name, code.co_firstlineno ) )


