# PatternTestCase
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


import string, re
import Bio.Prosite.Pattern
import UnitTests.UnitTestCase
import UnitTests.UnitTestSuite
import UnitTests.UnitTestResults



class PatternTestCase( UnitTests.UnitTestCase.UnitTestCase ):

    def setup( self ):

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

    """
    Title      test_conversion_1()
    Usage
    Function   Test a conversion of a prosite pattern to a regular expression.

    Returns    Results of the test
    Argument   None
    """

    def test_conversion_1( self ):

        pattern = '[LIV]-G-{P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-[GSTACLIVMFY]-' \
        + 'x(5,18)-[LIVMFYWCSTAR]-[AIVP]-[LIVMFAGCKR]-K.'
        print pattern
        actual = Bio.Prosite.Pattern.prosite_to_re( pattern )
        print actual
        expected = '[LIV]G[^P]G[^P][FYWMGSTNH][SGA][^PW][LIVCAT][^PD].[GSTACLIVMFY]' \
        + '.{5,18}[LIVMFYWCSTAR][AIVP][LIVMFAGCKR]K'
	return ( self.assert_equals( expected, actual  ) )

    """
    Title      test_conversion_2()
    Usage
    Function   Test a conversion of a prosite pattern to a regular expression.

    Returns    Results of the test
    Argument   None
    """

    def test_conversion_2( self ):

        pattern = '[IV]-x-D-S-[GAS]-[GASC]-[GAST]-[GA]-T.'
        print pattern
        actual = Bio.Prosite.Pattern.prosite_to_re( pattern )
        print actual
        expected = '[IV].DS[GAS][GASC][GAST][GA]T'
	return ( self.assert_equals( expected, actual  ) )


    """
    Title      test_conversion_3()
    Usage
    Function   Test a conversion of a prosite pattern to a regular expression.

    Returns    Results of the test
    Argument   None
    """

    def test_conversion_3( self ):

        pattern = 'G-[LIVM]-x(3)-E-[LIV]-T-[LF]-R.'
        print pattern
        actual = Bio.Prosite.Pattern.prosite_to_re( pattern )
        print actual
        expected = 'G[LIVM].{3}E[LIV]T[LF]R'
	return ( self.assert_equals( expected, actual  ) )

    """
    Title      test_conversion_4()
    Usage
    Function   Test a conversion of a prosite pattern to a regular expression.

    Returns    Results of the test
    Argument   None
    """

    def test_conversion_4( self ):

        pattern = '[DESH]-x(4,5)-[STVG]-x-[AS]-[FYI]-K-[DLIFSA]-[RVMF]-[GA]-[LIVMGA].'
        print pattern
        actual = Bio.Prosite.Pattern.prosite_to_re( pattern )
        print actual
        expected = '[DESH].{4,5}[STVG].[AS][FYI]K[DLIFSA][RVMF][GA][LIVMGA]'
	return ( self.assert_equals( expected, actual  ) )

    """
    Title      test_conversion_5()
    Usage
    Function   Test a conversion of a prosite pattern to a regular expression.

    Returns    Results of the test
    Argument   None
    """

    def test_conversion_5( self ):

        pattern = 'W-[IV]-[STA]-[RK]-x-[DE]-Y-[DNE]-[DE].'
        print pattern
        actual = Bio.Prosite.Pattern.prosite_to_re( pattern )
        print actual
        expected = 'W[IV][STA][RK].[DE]Y[DNE][DE]'
	return ( self.assert_equals( expected, actual  ) )

    """
    Title      test_verify_pattern_1()
    Usage
    Function   Test verification of a pattern

    Returns    Results of the test
    Argument   None
    """

    def test_verify_pattern_1( self ):

        pattern = 'W-[IV]-[STA]-[RK]-x-[DE]-Y-[DNE]-[DE].'
        print pattern
        is_ok = Bio.Prosite.Pattern.verify_pattern( pattern )
	return ( self.assert_condition( is_ok == 1, 'failed to accept a good pattern'  ) )

    """
    Title      test_verify_pattern_2()
    Usage
    Function   Test verification of a pattern

    Returns    Results of the test
    Argument   None
    """

    def test_verify_pattern_2( self ):

        pattern = 'W-[IV]-[STA*]-[RK]-x-[DE]-Y-[DNE]-[DE].'
        print pattern
        is_ok = Bio.Prosite.Pattern.verify_pattern( pattern )
	return ( self.assert_condition( is_ok != 1, 'failed to reject a bad pattern'  ) )

    """
    Title      test_verify_pattern_3()
    Usage
    Function   Test verification of a pattern

    Returns    Results of the test
    Argument   None
    """

    def test_verify_pattern_3( self ):

        pattern = 'W-[IV]-[STA-[RK]-x-[DE]-Y-[DNE]-[DE].'
        print pattern
        is_ok = Bio.Prosite.Pattern.verify_pattern( pattern )
	return ( self.assert_condition( is_ok != 1, 'failed to reject a bad pattern'  ) )

    """
    Title      test_verify_pattern_4()
    Usage
    Function   Test verification of a pattern

    Returns    Results of the test
    Argument   None
    """

    def test_verify_pattern_4( self ):

        pattern = '[LIV]-G-P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-[GSTACLIVMFY]-' \
        + 'x(5,18)-[LIVMFYWCSTAR]-[AIVP]-[LIVMFAGCKR]-K.'
        print pattern
        is_ok = Bio.Prosite.Pattern.verify_pattern( pattern )
	return ( self.assert_condition( is_ok != 1, 'failed to reject a bad pattern'  ) )

    """
    Title      test_verify_pattern_5()
    Usage
    Function   Test verification of a pattern

    Returns    Results of the test
    Argument   None
    """

    def test_verify_pattern_5( self ):

        pattern = '[LIV]-G-{P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-[GSTACLIVMFY]-' \
        + 'x(5,18)-[LIVMFYWCSTAR]-[AIVP]-[LIVMFAGCKR]-K.'
        print pattern
        is_ok = Bio.Prosite.Pattern.verify_pattern( pattern )
	return ( self.assert_condition( is_ok == 1, 'failed to accept a good pattern'  ) )

#invoke_suite
#
# Title     : invoke_suite()
# Usage     :
# Function  : create a suite of tests
#           :
# Returns   :
# Argument  : None

    def invoke_suite( self, suite_name ):

        test_results = UnitTests.UnitTestResults.UnitTestResults( )
        suite = UnitTests.UnitTestSuite.UnitTestSuite( suite_name )
	suite.build_suite( suite_name, self )
        suite.run( test_results )
#    open  $fh, ">unit.log";
        test_results.print_failures( )
        test_results.print_errors( )
#    close $fh;

