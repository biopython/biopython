# JustOneTestCase
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
import Defect
import UnitTestResults
import UnitTestCase
import UnitTestCaller



class JustOneTestCase( UnitTestCase.UnitTestCase ):

    def setup( self ):
        self.a = 0

	print 'setup'

    def tear_down( self ):
        """
	Title      tear_down()
        Usage      unit_test_case.teardown()
        Function   Clean up after test case

        Returns
        Argument   None
        """
        print 'tear_down'

    #Title     : test_condition( self )
    #Usage     :
    #Function  : Test for a condition
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_condition( self ):

        return ( self.assert_condition( self.a, 'failed condition' ) )

# call_test
#
# Title     : call_test
# Usage     :
# Function  : call a single test
#           :
# Returns   :
# Argument  : None

    def call_test( self ):

        test_results = UnitTestResults.UnitTestResults( )
	caller = UnitTestCaller.UnitTestCaller( 'test_condition', self.test_condition, self )
	caller.run( test_results )
#    open  $fh, ">unit.log";
        test_results.print_failures( )
#    close $fh;

