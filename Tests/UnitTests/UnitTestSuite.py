# UnitTestSuite.py
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
import string
import UnitTestExceptions
import Defect
import UnitTestResults
import UnitTest
import UnitTestCaller

class UnitTestSuite( UnitTest.UnitTest ):



    def __init__( self, name ):
        """
        Title     new
        Usage     unit_test_suite = UnitTestSuite.( "suite" );
        Function  Initializes a new UnitTestSuite object.
        Example   See usage
        Returns   Bio::UnitTests::UnitTestSuite object
        Argument
        """


        try:
            if not( type( name ) == type( '' ) ):
                raise UnitTestExceptions.BadParameter( 'method name must''\
		     be a string' )
        except UnitTestExceptions.BadParameter, bad_parameter:
            print '%s' % ( bad_parameter.message )
        else:
            self.name = name
            self.suite = []

    def add_test( self, caller ):
        """
        Title     add_test()
        Usage     unit_test_suite.add_test( self, caller )
        Function  Add a test to the suite

        Returns
        Argument  : test caller
        """

        try:
	    caller_name = caller.name
            if not( type( caller_name ) == type( '' ) ):
                raise UnitTestExceptions.BadParameter( 'method name must''\
		     be a string' )
            if not ( caller_name[ :5 ] == 'test_' ):
                raise UnitTestExceptions.BadParameter( 'method name must''\
		     begin with "test"' )

        except UnitTestExceptions.BadParameter, bad_parameter:
            print '%s' % ( bad_parameter.message )
        else:
            self.suite.append( caller )


    def run( self, test_results ):
        """
        Title     run()
        Usage     unit_test_suite.run( self, test_results )
        Function  Run the suite of tests

        Returns
        Argument  Test results
        """


	for caller in self.suite:
	    caller.run( test_results )

        """
 	Title     build_suite()
        Usage
        Function  create a suite of tests

        Returns
        Argument  : None
        """

    def build_suite( self, modname, test_object ):

        print 'build_suite'
        module = sys.__dict__[ 'modules' ][ modname ]
	for entry in module.__dict__[ modname ].__dict__.keys():
	    if( entry[ :5 ] == 'test_' ):
                method = eval( 'test_object.' + entry )
                caller = UnitTestCaller.UnitTestCaller( entry, method, test_object )
                self.add_test( caller )


