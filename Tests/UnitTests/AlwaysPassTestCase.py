# AlwaysPassTestCase
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
import UnitTestSuite



class AlwaysPassTestCase( UnitTestCase.UnitTestCase ):

    def setup( self ):
        self.a = 1
	self.b = 1
	self.texta = 'python'
        self.textb = 'python'
	self.languages_a = [ 'perl', 'python', 'tcl' ]
	self.languages_b = [ 'perl', 'python', 'tcl' ]
        self.hash_a = { 'dog' : 'bark', 'bird' : 'chirp', 'cow': 'moo' }
        self.hash_b = { 'dog' : 'bark', 'bird' : 'chirp', 'cow': 'moo' }

	print 'setup'

    def tear_down( self ):
        """
	Title      tear_down()
        Usage      unit_test_case.teardown()
        Function   Clean up after test case

        Returns
        Argument   None
        """
        pass

    #Title     : test_condition( self )
    #Usage     :
    #Function  : Test for a condition
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_condition( self ):

        a = 0;
        return ( self.assert_condition( self.a, 'failed condition' ) )

    #Title     : test_equals()
    #Usage     :
    #Function  : Test for an equality between two values
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_equals( self ):

        return ( self.assert_equals( self.a, self.b ) )

    #Title     : test_list_equals()
    #Usage     :
    #Function  : Test for an equality between two lists
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_list_equals( self ):

        return ( self.assert_equals( self.languages_a, self.languages_b ) )

    #Title     : test_string_equals()
    #Usage     :
    #Function  : Test for an equality between two strings
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_string_equals( self ):

        return ( self.assert_equals( self.texta, self.textb ) )

    #Title     : test_hash_equals()
    #Usage     :
    #Function  : Test for an equality between two hashes
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_hash_equals( self ):

        return ( self.assert_equals( self.hash_a, self.hash_b ) )


#invoke_suite
#
# Title     : invoke_suite()
# Usage     :
# Function  : create a suite of tests
#           :
# Returns   :
# Argument  : None

    def invoke_suite( self, suite_name ):

        test_results = UnitTestResults.UnitTestResults( )
        suite = UnitTestSuite.UnitTestSuite( suite_name )
	suite.build_suite( suite_name, self )
        suite.run( test_results )
#    open  $fh, ">unit.log";
        test_results.print_failures( )
#    close $fh;

