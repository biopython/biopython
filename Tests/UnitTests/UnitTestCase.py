# UnitTestCase.
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
import string
import UnitTestExceptions
import Defect
import UnitTestResults
import UnitTest



class UnitTestCase( UnitTest.UnitTest ):



    def __init__( self, name ):
        """
        Title      __init__
        Usage      unit_test_case = UnitTestCase( name )
        Function   Initializes a new UnitTestCase object.
        Example    See usage
        Returns
        Argument
        """

        if( type( name ) == type( '' ) ):
	    self.name = name


    def assert_condition( self, condition, message ):
        """
        Title      assert_condition
        Usage      unit_test_case.assert_condition( condition, message )
        Function   Return failure if the condition is false

        Returns    Failure
        Argument  : condition, message
        """

        try:
            failure = 0
	    if not condition:
		raise UnitTestExceptions.FailedUnitTest( message )

        except UnitTestExceptions.FailedUnitTest, failed_unit_test:
	    trace_back = sys.exc_info()[ 2 ]
            trace = traceback.extract_tb( trace_back )
	    text = failed_unit_test.message
	    failure = Defect.Defect( text, trace_back )

	return( failure )

    def assert_equals( self, expected, actual ):
        """
        Title      assert_equals
        Usage      unit_test_case.assert_equals( self, expected, actual )
        Function   Return failure if expected does not equal actual

        Returns   reference to failure object
        Argument  expected, actual, message
        """


        try:
            failure = None
            if not( expected == actual ):
                message = 'expected is %s actual is %s' % ( str( expected ), \
                    str( actual ) )
	        raise UnitTestExceptions.FailedUnitTest( message )
        except UnitTestExceptions.FailedUnitTest, failed_unit_test:
            trace_back = sys.exc_traceback
	    text = failed_unit_test.message
	    failure = Defect.Defect( text, trace_back )

        return( failure )


    def name( self, name = '' ):
        """
        Title     name()
        Usage     unit_test_case.name()
                  unit_test_case.name( name )
        Function  Return/set the name of the test case

        Returns   Name of the test case
        Argument  None
        """


        if( name and ( type( name ) == type( '' ) ) ):
            self.name = name
	    return( self.name )


    def setup( self ):

        """
        Title      setup()
        Usage      unit_test_case.setup()
        Function   Prepare fortest case

        Returns
        Argument   None
        """
        pass

    def tear_down( self ):
        """
	Title      tear_down()
        Usage      unit_test_case.teardown()
        Function   Clean up after test case

        Returns
        Argument   None
        """
        pass

    def run_test( self ):
        """
        Title      run_test
        Usage
        Function   Run the test case.

        Returns    failure
        Argument   None
        """

        failure = None
        return failure


    def run( self, test_results ):
        """
        Title      run
        Usage      run( self, test_results )
        Function   Run the test case.

        Returns
        Argument   result
        """

        self.setup()
	try:
	    failure = self.run_test()
            if failure:
	        print 'failure' + '\n'
		test_results.add_failure( failure )
            else:
                print 'success' + '\n'
        except:
	    trace_back = sys.exc_info()[ 2 ]
	    error = Defect.Defect( sys.exc_type, trace_back )
	    test_results.add_error( error )
            print 'error'
	self.tear_down()

