# UnitTestSuite.pm
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
import AssertionFailure
import UnitTestResults
import UnitTest
import UnitTestCaller
class UnitTestSuite( UnitTest.UnitTest ):



    #Title     : new
    #Usage     : unit_test_suite = UnitTestSuite.( "suite" );
    #Function  : Initializes a new UnitTestSuite object.
    #Example   : See usage
    #Returns   : Bio::UnitTests::UnitTestSuite object
    #Argument  :

    def __init__( self, name ):

        try:
            if not( type( name ) == type( '' ) ):
                raise UnitTestExceptions.BadParameter( 'method name must''\
		     be a string' )
        except UnitTestExceptions.BadParameter, bad_parameter:
            print '%s' % ( bad_parameter.message )
        else:
            self.name = name
            self.suite = []

    #Title     : add_test()
    #Usage     : unit_test_suite.add_test( self, caller )
    #Function  : Add a test to the suite
    #          :
    #Returns   :
    #Argument  : test caller

    def add_test( self, caller ):
        try:
	    name = caller.name
            if not( name[ :4 ] == 'test' ):
                raise UnitTestExceptions.BadParameter( 'method name must''\
		     begin with "test"' )

        except UnitTestExceptions.BadParameter, bad_parameter:
            print '%s' % ( bad_parameter.message )
        else:
            self.suite.append( caller )


     #Title     : run()
     #Usage     : unit_test_suite.run( self, test_results )
     #Function  : Run the suite of tests
     #          :
     #Returns   :
     #Argument  : Test results

    def run( self, test_results ):

        for caller in self.suite:
	    caller.run( test_results )



