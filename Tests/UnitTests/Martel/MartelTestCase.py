# MartelTestCase
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


import string
import Bio.Martel
from Bio.Martel.Generate import generate
import TextTools
from TextTools.mxTextTools import tag
import UnitTests.UnitTestCase
import UnitTests.UnitTestSuite
import UnitTests.UnitTestResults



class MartelTestCase( UnitTests.UnitTestCase.UnitTestCase ):

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
    Title      test_a()
    Usage
    Function   Any

    Returns    Results of the test
    Argument   None
    """

    def test_a( self ):

        expression = Bio.Martel.Any( "ace" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "a", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )


    """
    Title      test_b()
    Usage
    Function   Any

    Returns    Results of the test
    Argument   None
    """

    def test_b( self ):

        expression = Bio.Martel.Any( "ace" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "c", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )


    """
    Title      test_c()
    Usage
    Function   Any

    Returns    Results of the test
    Argument   None
    """

    def test_c( self ):

        expression = Bio.Martel.Any( "ace" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "e", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_d()
    Usage
    Function   Any

    Returns    Results of the test
    Argument   None
    """

    def test_d( self ):

        expression = Bio.Martel.Any( "ace" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "b", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_e()
    Usage
    Function   AnyBut

    Returns    Results of the test
    Argument   None
    """

    def test_e( self ):

        expression = Bio.Martel.AnyBut( "ace" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "d", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_f()
    Usage
    Function   AnyBut

    Returns    Results of the test
    Argument   None
    """

    def test_f( self ):

        expression = Bio.Martel.AnyBut( "ace" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "a", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_g()
    Usage
    Function   Dot

    Returns    Results of the test
    Argument   None
    """

    def test_g( self ):

        expression = Bio.Martel.Expression.Dot()
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "d", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_h()
    Usage
    Function   Dot

    Returns    Results of the test
    Argument   None
    """

    def test_h( self ):

        expression = Bio.Martel.Expression.Dot()
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "2", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_i()
    Usage
    Function   Dot

    Returns    Results of the test
    Argument   None
    """

    def test_i( self ):

        expression = Bio.Martel.Expression.Dot()
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "&", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_j()
    Usage
    Function   Dot

    Returns    Results of the test
    Argument   None
    """

    def test_j( self ):

        expression = Bio.Martel.Expression.Dot()
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( ".", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_k()
    Usage
    Function   Dot

    Returns    Results of the test
    Argument   None
    """

    def test_k( self ):

        expression = Bio.Martel.Expression.Dot()
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "\n", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_l()
    Usage
    Function   AtEnd

    Returns    Results of the test
    Argument   None
    """

    def test_l( self ):

        expression = Bio.Martel.Expression.AtEnd()
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "a\n", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_m()
    Usage
    Function   AtEnd

    Returns    Results of the test
    Argument   None

    """

    def test_m( self ):

        expression = Bio.Martel.Expression.AtEnd()
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "ab\n", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_n()
    Usage
    Function   Str

    Returns    Results of the test
    Argument   None
    """

    def test_n( self ):

        expression = Bio.Martel.Str( "abcd" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "abcd", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_o()
    Usage
    Function   Str

    Returns    Results of the test
    Argument   None
    """

    def test_o( self ):

        expression = Bio.Martel.Str( "abcd" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "abcz", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_p()
    Usage
    Function   Str

    Returns    Results of the test
    Argument   None
    """

    def test_p( self ):

        expression = Bio.Martel.Str( "abcd" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "abc", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_q()
    Usage
    Function   Str

    Returns    Results of the test
    Argument   None
    """

    def test_q( self ):

        expression = Bio.Martel.Str( "abcd" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "Abcd", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_r()
    Usage
    Function   Literal

    Returns    Results of the test
    Argument   None
    """

    def test_r( self ):

        expression = Bio.Martel.Expression.Literal( "Z" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "Ze", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_s()
    Usage
    Function   Literal

    Returns    Results of the test
    Argument   None
    """

    def test_s( self ):

        expression = Bio.Martel.Expression.Literal( "Z" )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "Ye", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_t()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_t( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "l", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )


    """
    Title      test_u()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_u( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "m", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_v()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_v( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "n", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_w()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_w( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "K", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_x()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_x( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "VWX", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )


    """
    Title      test_y()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_y( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "L", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


    """
    Title      test_z()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_z( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "N", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_a1()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_a1( self ):

        alt1 = Bio.Martel.Expression.Literal( 'K' )
        alt2 = Bio.Martel.Str( 'VWX' )
        alt3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Alt( alt1, alt2, alt3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "UVWX", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_b1()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_b1( self ):

        exp1 = Bio.Martel.Expression.Literal( 'K' )
        exp2 = Bio.Martel.Str( 'VWX' )
        exp3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Seq( exp1, exp2, exp3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "KVWXl", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_c1()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_c1( self ):

        exp1 = Bio.Martel.Expression.Literal( 'K' )
        exp2 = Bio.Martel.Str( 'VWX' )
        exp3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Seq( exp1, exp2, exp3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "KVWXm", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_d1()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_d1( self ):

        exp1 = Bio.Martel.Expression.Literal( 'K' )
        exp2 = Bio.Martel.Str( 'VWX' )
        exp3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Seq( exp1, exp2, exp3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "KVWXn", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_e1()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_e1( self ):

        exp1 = Bio.Martel.Expression.Literal( 'K' )
        exp2 = Bio.Martel.Str( 'VWX' )
        exp3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Seq( exp1, exp2, exp3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "lVWXK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


    """
    Title      test_f1()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_f1( self ):

        exp1 = Bio.Martel.Expression.Literal( 'K' )
        exp2 = Bio.Martel.Str( 'VWX' )
        exp3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Seq( exp1, exp2, exp3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "VWXKl", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_g1()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_g1( self ):

        exp1 = Bio.Martel.Expression.Literal( 'K' )
        exp2 = Bio.Martel.Str( 'VWX' )
        exp3 = Bio.Martel.Any( 'lmn' )
        expression = Bio.Martel.Seq( exp1, exp2, exp3 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "KmVWX", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_h1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_h1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 2, 2 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzAzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_i1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_i1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 2, 2 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_j1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_j1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 2, 2 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzAzAzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


    """
    Title      test_k1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_k1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 2 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzAzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_l1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_l1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 2 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_m1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_m1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 2 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzAzAzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_n1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_n1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 0, 1 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "K", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_o1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_o1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 0, 1 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )


    """
    Title      test_p1()
    Usage
    Function   MaxRepeat

    Returns    Results of the test
    Argument   None
    """

    def test_p1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 0, 1 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzAzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_z1()
    Usage
    Function   Maxrepeat

    Returns    Results of the test
    Argument   None
    """

    def test_z1( self ):

        exp1 = Bio.Martel.Str( 'Az' )
        exp2 = Bio.Martel.Expression.Literal( 'K' )

        exp3 = Bio.Martel.MaxRepeat( exp1, 0, 0 )
        expression = Bio.Martel.Seq( exp3, exp2 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "AzK", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_a2()
    Usage
    Function   Assert

    Returns    Results of the test
    Argument   None
    """

    def test_a2( self ):

        exp1 = Bio.Martel.Str( 'bC' )
        exp2 = Bio.Martel.Expression.Literal( 'V' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Expression.Assert( exp3 )
        exp5 = Bio.Martel.Expression.Literal( 'b' )
        expression = Bio.Martel.Seq( exp4, exp5 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "bCV", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_b2()
    Usage
    Function   Assert

    Returns    Results of the test
    Argument   None
    """

    def test_b2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Str( 'ij' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Expression.Assert( exp3, 0 )
        exp5 = Bio.Martel.Str( 'ti' )
        expression = Bio.Martel.Seq( exp4, exp5 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "tij", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_c2()
    Usage
    Function   Assert

    Returns    Results of the test
    Argument   None
    """

    def test_c2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Str( 'ij' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Expression.Assert( exp3, 1 )
        exp5 = Bio.Martel.Str( 'ti' )
        expression = Bio.Martel.Seq( exp4, exp5 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "tij", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_d2()
    Usage
    Function   Assert

    Returns    Results of the test
    Argument   None
    """

    def test_d2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Str( 'ij' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Expression.Assert( exp3, 0 )
        exp5 = Bio.Martel.Str( 'ti' )
        expression = Bio.Martel.Seq( exp4, exp5 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "tiJ", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


    """
    Title      test_e2()
    Usage
    Function   Assert

    Returns    Results of the test
    Argument   None
    """

    def test_e2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Str( 'ij' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Expression.Assert( exp3, 2 )
        exp5 = Bio.Martel.Str( 'ti' )
        expression = Bio.Martel.Seq( exp4, exp5 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "tij", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_f2()
    Usage
    Function   Assert

    Returns    Results of the test
    Argument   None
    """

    def test_f2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Str( 'ij' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Expression.Assert( exp3, 3 )
        exp5 = Bio.Martel.Str( 'ti' )
        expression = Bio.Martel.Seq( exp4, exp5 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "tij", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_g2()
    Usage
    Function   Group, GroupRef

    Returns    Results of the test
    Argument   None
    """

    def test_g2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Str( 'pQrsT' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Group( 'item', exp3 )
        exp5 = Bio.Martel.Str( 'ab' )
        expression = Bio.Martel.Seq( exp4, exp5, \
            Bio.Martel.Expression.GroupRef( 'item' ) )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "tpQrsTabtpQrsT", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_h2()
    Usage
    Function   Group, GroupRef

    Returns    Results of the test
    Argument   None
    """

    def test_h2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Str( 'pQrsT' )

        exp3 = Bio.Martel.Seq( exp1, exp2 )
        exp4 = Bio.Martel.Group( 'item', exp3 )
        exp5 = Bio.Martel.Str( 'ab' )
        expression = Bio.Martel.Seq( exp4, exp5, \
            Bio.Martel.Expression.GroupRef( 'item' ) )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "tpQrsTabtPqRSt", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_i2()
    Usage
    Function   Group, GroupRef

    Returns    Results of the test
    Argument   None
    """

    def test_i2( self ):

        exp1 = Bio.Martel.Expression.Literal( 't' )
        exp2 = Bio.Martel.Expression.Literal( 'T' )
        exp3 = Bio.Martel.Alt( exp1, exp2 )
        exp4 = Bio.Martel.Str( '01234' )

        exp5 = Bio.Martel.Seq( exp3, exp4 )
        exp6 = Bio.Martel.Group( 'item', exp5 )
        exp7 = Bio.Martel.Str( '56789' )
        exp8 = Bio.Martel.Expression.GroupRef( 'item' )
        expression = Bio.Martel.Seq( exp6, exp7, exp8, exp7, exp8, exp7 )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "t0123456789T0123456789t0123456789", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_j2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_j2( self ):

        exp1 = Bio.Martel.Expression.Literal( 'P' )
        exp2 = Bio.Martel.Any( 'ghijk' )
        exp3 = Bio.Martel.Str( '!!!!!' )
        expression = exp1 + exp2 + exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "Pk!!!!!", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_k2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_k2( self ):

        exp1 = Bio.Martel.Any( 'GHIjk' )
        exp2 = Bio.Martel.Str( ':) :)' )
        expression = exp1 + exp2
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "G:) :)", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_l2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_l2( self ):

        exp1 = Bio.Martel.Any( 'GHIjk' )
        exp2 = Bio.Martel.Str( ':) :)?' )
        expression = exp1 + exp2
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "K:) :)?", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_m2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_m2( self ):

        exp1 = Bio.Martel.Any( 'GHIjk' )
        exp2 = Bio.Martel.Str( ':) :)' )
        expression = exp1 + exp2
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "G:):)", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_n2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_n2( self ):

        exp1 = Bio.Martel.Str( '$%@#/\\' )
        exp2 = Bio.Martel.Any( 'GHIjk' )
        exp3 = Bio.Martel.Str( '(*+_-=' )
        expression = exp1 + exp2 + exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( '$%@#/\\I*+_-=', tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_o2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_o2( self ):

        exp1 = Bio.Martel.Str( '$%@#/\\' )
        exp2 = Bio.Martel.Any( 'GHIjk' )
        exp3 = Bio.Martel.Str( '(*+_-=' )
        expression = exp1 + exp2 + exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "$%@#\\/I*+_-=", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


    """
    Title      test_p2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_p2( self ):

        exp1 = Bio.Martel.Str( '$%@#/\\' )
        exp2 = Bio.Martel.Any( 'GHIjk' )
        exp3 = Bio.Martel.Str( '(*+_-=' )
        expression = exp1 + exp2 + exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "$%@#/\\I*+-_=", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_q2()
    Usage
    Function   Overloaded or

    Returns    Results of the test
    Argument   None
    """

    def test_q2( self ):

        exp1 = Bio.Martel.Expression.Literal( 'P' )
        exp2 = Bio.Martel.Any( 'gHiJk' )
        exp3 = Bio.Martel.Str( '.<>,;' )
        expression = exp1 | exp2 | exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( ".<>,;", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_r2()
    Usage
    Function   Overloaded or

    Returns    Results of the test
    Argument   None
    """

    def test_r2( self ):

        exp1 = Bio.Martel.Expression.Literal( 'P' )
        exp2 = Bio.Martel.Any( 'gHiJk' )
        exp3 = Bio.Martel.Str( '.<>,;' )
        expression = exp1 | exp2 | exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "H", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_s2()
    Usage
    Function   Overloaded or

    Returns    Results of the test
    Argument   None
    """

    def test_s2( self ):

        exp1 = Bio.Martel.Expression.Literal( 'P' )
        exp2 = Bio.Martel.Any( 'gHiJk' )
        exp3 = Bio.Martel.Str( '.<>,;' )
        expression = exp1 | exp2 | exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "P", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_t2()
    Usage
    Function   Overloaded or

    Returns    Results of the test
    Argument   None
    """

    def test_t2( self ):

        exp1 = Bio.Martel.Expression.Literal( 'P' )
        exp2 = Bio.Martel.Any( 'gHiJk' )
        exp3 = Bio.Martel.Str( '.<>,;' )
        expression = exp1 | exp2 | exp3
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "j", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


    """
    Title      test_v2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_v2( self ):

        exp1 = Bio.Martel.Any( 'GHIjk' )
        exp2 = Bio.Martel.Str( '[]{}~' )
        expression = exp1 | exp2
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "[]{}~", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_w2()
    Usage
    Function   Overloaded plus

    Returns    Results of the test
    Argument   None
    """

    def test_w2( self ):

        exp1 = Bio.Martel.Any( 'GHIjk' )
        exp2 = Bio.Martel.Str( '[]{}~' )
        expression = exp1 | exp2
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "[]{}", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


    """
    Title      test_x2()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_x2( self ):

        alt = Bio.Martel.Str( 'uvwx' )
        expression = Bio.Martel.Alt( alt )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "uvwx", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_y2()
    Usage
    Function   Alt

    Returns    Results of the test
    Argument   None
    """

    def test_y2( self ):

        alt = Bio.Martel.Str( 'uvwx' )
        expression = Bio.Martel.Alt( alt )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "Uvwx", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )

    """
    Title      test_z2()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_z2( self ):

        exp = Bio.Martel.Expression.Literal( '9' )
        expression = Bio.Martel.Seq( exp )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "9", tagtable )[ 0 ]
	return ( self.assert_condition( success == 1, "Failed" ) )

    """
    Title      test_a3()
    Usage
    Function   Seq

    Returns    Results of the test
    Argument   None
    """

    def test_a3( self ):

        exp = Bio.Martel.Expression.Literal( '9' )
        expression = Bio.Martel.Seq( exp )
        print expression
        tagtable, want_flg = Bio.Martel.Generate.generate( expression )
        success = tag( "8", tagtable )[ 0 ]
	return ( self.assert_condition( success == 0, "Failed" ) )


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

