# UnitTestCaller.pm
#
#
# MODIFICATION NOTES: See bottom of file.

# Copyright (c) 1999 Katharine Lindner
# This module is free software; you can redistribute it and/or modify
# it under the same terms as Python

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
from types import *
import string
import UnitTestExceptions
import Defect
import UnitTestCase



class UnitTestCaller( UnitTestCase.UnitTestCase ):

    def __init__( self, name, method, test_object ):
        """
        Title      __init__
        Usage      unit_test_caller = UnitTestCaller->new( method, test_object )
        Function   initializes the UnitTestCaller object
        Example    See usage
        Returns
        Argument
        """
        try:
	    if ( type( name ) != StringType ):
                raise UnitTestExceptions.BadParameter( 'method name must be a string' )
            if not( type( method ) == MethodType ):
                raise UnitTestExceptions.BadParameter( 'method must be code' )
            if not( name[ :5 ] == 'test_' ):
                raise UnitTestExceptions.BadParameter( 'method name must begin with "test"' )
        except UnitTestExceptions.BadParameter, bad_parameter:
            print '%s' % ( bad_parameter.message )
        else:
            self.name = name
	    self.method = method.im_func
            self.test_object = test_object

        """
        Title     name( self, name = '' )
        Usage     unit_test_caller.name()
                  unit_test_caller.name( name )
        Function  Return/get the name of the test case

        Returns   Name of the test case
        Argument  None
        """

    def name( self, name = '' ):
        """
        Title     name( self, name = '' )
        Usage     unit_test_caller.name()
                  unit_test_caller.name( name )
        Function  Return/get the name of the test case

        Returns   Name of the test case
        Argument  None
        """


        if( name and ( type( name ) == type( '' ) ) ):
            self.name = name
        return( self.name )

    def setup( self ):
        """Title      setup
        Usage      unit_test_case.setup()
        Function   Prepare for a test case

        Returns
        Argument   None
       """

	test_object = self.test_object
        test_object.setup()

    def tear_down( self ):
        """Title      tear_down
        Usage      unit_test_case.tear_down()
        Function   Clean up after a test case

        Returns
        Argument   None
       """

	test_object = self.test_object
        test_object.tear_down()



    def run_test( self ):
        """
        Title     run_test( self )
        Usage
        Function  Run the test case.

        Returns   outcome of the test
        Argument  Failure/ Success
        """

	method_name = self.name
        method = self.method
	test_object = self.test_object
	print '%s' % ( method_name )
 	return( method( test_object ) )

