# UnitTestResults.py
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
import UnitTestExceptions
import Defect

class UnitTestResults:
# biopython test results object


    def __init__( self ):
        """
        Title     __initialize__
        Usage
        Function  initialize a new object.
        Example
        Returns
        Argument
        """


        self.failures = []
        self.errors = []
        self.warnings = []

    def get_failures( self ):
        """
        Title     get_failures()
        Usage
        Function  Return the list of assertion failures

        Returns   List of failures
        Argument  None
        """


        return ( self.failures )

    def get_errors( self ):
        """
        Title     get_errors()
        Usage
        Function  Return the list of errors

        Returns   List of errors
        Argument  None
        """


	return ( self.errors )

    def get_warnings( self ):
        """
        Title     get_warnings()
        Usage
        Function  Return the list of warnings

        Returns   List of warnings
        Argument  None
        """


        return ( self.warnings )

    def count_failures( self ):
        """
        Title     count_failures()
        Usage
        Function  Returns the number of assertion failures

        Returns   number of failures
        Argument  none
        """


        return len( self.get_failures() )

    def count_errors( self ):
        """
        Title     count_errors()
        Usage
        Function  Returns the number of errors

        Returns   number of errors
        Argument  : none
        """


        return len( self.get_errors() )

    def count_warnings( self ):
        """
        Title     count_warnings()
        Usage
        Function  Returns the number of warnings

        Returns   number of warnings
        Argument  none
        """


        return len( self.get_warnings() )

    def add_failure( self, failure ):
        """
        Title     add_failure()
        Usage
        Function  Add a failure to the list of assertion failures

        Returns
        Argument  failure
        """


	self.get_failures().append( failure )

    def add_error( self, error ):
        """
        Title     add_error()
        Usage
        Function  Add an error to the list of errors

        Returns
        Argument  error
        """


        self.get_errors().append( error )

    def add_warning( self, warning ):
        """
        Title     add_warning()
        Usage
        Function  Add a warning to the list of warnings

        Returns
        Argument  warning
        """


        self.get_warnings().append( warning )

    def was_successful( self ):
        """
        Title     was_successful()
        Usage
        Function  Returns false if any failures or errors

        Returns   Returns false if any failures or errors
        Argument  none
        """


        return  ( not( self.count_failures() ) and  \
	        not( self.count_errors() ) )

    def print_failures( self, outfile = sys.stdout ):
        """
        Title     print_failures()
        Usage     print_failures
                  print_failures( file_handle );
        Function  Prints assertion failure messages

        Returns
        Argument  file handle
        """


        for failure in self.get_failures():
            failure.print_message( outfile )
            failure.print_trace( outfile )

    def print_errors( self, outfile = sys.stdout ):
        """
        Title     print_errors()
        Usage     print_errors
                  print_errors( file_handle )
        Function  Prints error messages

        Returns
        Argument  file handle
        """


        for error in self.get_errors():
            error.print_message( outfile )
            error.print_trace( outfile )

    def print_warnings( self, outfile = sys.stdout ):
        """
        Title     print_warnings()
        Usage     print_warnings
                  print_warnings( file_handle )
        Function  Prints warning messages

        Returns
        Argument  file handle
        """


        for warning in self.get_warnings():
            outfile.write(  '%s\n' % ( warning ) )


