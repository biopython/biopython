# NeverPassTestCase
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



class NeverPassTestCase( UnitTestCase.UnitTestCase ):

    #Title     : test_condition( self )
    #Usage     :
    #Function  : Test for a condition
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_condition( self ):

        a = 0;
        return ( self.assert_condition( a, 'failed condition' ) )

    #Title     : test_equals()
    #Usage     :
    #Function  : Test for an equality between two values
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_equals( self ):

        a = 0
        b = 1

        return ( self.assert_equals( a, b ) )

    #Title     : test_list_equals()
    #Usage     :
    #Function  : Test for an equality between two lists
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_list_equals( self ):

        a = ( 'ruby', 'sapphire', 'diamond' )
        b = ( 'ruby', 'sapphire', 'zircon' )

        return ( self.assert_equals( a, b ) )

    #Title     : test_string_equals()
    #Usage     :
    #Function  : Test for an equality between two strings
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_string_equals( self ):

        a = 'python'
        b = 'perl'

        return ( self.assert_equals( a, b ) )

    #Title     : test_hash_equals()
    #Usage     :
    #Function  : Test for an equality between two hashes
    #          :
    #Returns   : Results of the test
    #Argument  : None

    def test_hash_equals( self ):

        a = { 'dog' : 'bark', 'bird' : 'chirp', 'cow': 'moo' }
        b = { 'dog' : 'howl', 'bird' : 'chirp', 'cow': 'moo' }

        return ( self.assert_equals( a, b ) )


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

