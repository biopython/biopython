# SeqTestCase
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
from Bio import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Tools import Translate
import UnitTests.UnitTestCase
import UnitTests.UnitTestSuite
import UnitTests.UnitTestResults



class SeqTestCase( UnitTests.UnitTestCase.UnitTestCase ):

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
    Title      test_mod1_a()
    Usage
    Function   Test a sequence of length mod 1

    Returns    Results of the test
    Argument   None
    """

    def test_mod1_a( self ):

        s = "T"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_id[1]
        print 'ready to translate'
        protein = trans.translate_to_stop(dna)
        print protein.tostring()
	return ( self.assert_equals( '', protein.tostring() ) )

    """
    Title      test_mod2_b()
    Usage
    Function   Test a sequence of length mod 2

    Returns    Results of the test
    Argument   None
    """

    def test_mod2_b( self ):

        s = "TC"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_id[1]
        print 'ready to translate'
        protein = trans.translate_to_stop(dna)
        print protein.tostring()
	return ( self.assert_equals( '', protein.tostring() ) )

    """
    Title      test_mod1_c()
    Usage
    Function   Test a sequence of length mod 1

    Returns    Results of the test
    Argument   None
    """

    def test_mod1_c( self ):

        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_id[1]
        print 'ready to translate'
        protein = trans.translate_to_stop(dna)
        print protein.tostring()
	return ( self.assert_equals( 'ENSFSLDFL', protein.tostring() ) )

    """
    Title      test_last_id()
    Usage
    Function   Translate from table 15

    Returns    Results of the test
    Argument   None
    """

    def test_last_id( self ):

        s = "GAA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_id[ 15 ]
        print 'ready to translate'
        protein = trans.translate_to_stop(dna)
        print protein.tostring()
	return ( self.assert_equals( 'E', protein.tostring() ) )

    """
    Title      test_by name()
    Usage
    Function   Translate from table for Mitochondrial Mold

    Returns    Results of the test
    Argument   None
    """

    def test_by_name( self ):

        s = "ATA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_name[ 'SGC1' ]
#        trans = Translate.unambiguous_dna_by_name[ 'Vertebrate Mitochondrial' ]
        print 'ready to translate'
        protein = trans.translate_to_stop(dna)
        print protein.tostring()
	return ( self.assert_equals( 'M', protein.tostring() ) )

    """
    Title      test_by name2()
    Usage
    Function   Translate from table for Mitochondrial Mold

    Returns    Results of the test
    Argument   None
    """

    def test_by_name2( self ):

        s = "ATA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_name[ 'Vertebrate Mitochondria' ]
        print 'ready to translate'
        protein = trans.translate_to_stop(dna)
        print protein.tostring()
	return ( self.assert_equals( 'M', protein.tostring() ) )


    """
    Title      test_by name3()
    Usage
    Function   Translate from table for Echinoderm Mitochondrial

    Returns    Results of the test
    Argument   None
    """

    def test_by_name3( self ):

        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        trans = Translate.unambiguous_dna_by_name[ 'SGC8' ]
#        trans = Translate.unambiguous_dna_by_name[ 'Vertebrate Mitochondrial' ]
        protein = trans.translate_to_stop(dna)
        print protein.tostring()
	return ( self.assert_equals( 'ENSFSLDFLWNPSPSNDAWDSSY', protein.tostring() ) )



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

