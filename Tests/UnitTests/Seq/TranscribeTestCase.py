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
from Bio.Tools import Transcribe
import UnitTests.UnitTestCase
import UnitTests.UnitTestSuite
import UnitTests.UnitTestResults



class TranscribeTestCase( UnitTests.UnitTestCase.UnitTestCase ):

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
    Function   Transcribe

    Returns    Results of the test
    Argument   None
    """

    def test_a( self ):

        s = "ATA"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        print 'ready to transcribe'
        rna = Transcribe.unambiguous_transcriber.transcribe( dna )

	return ( self.assert_equals( 'AUA', rna.tostring() ) )

    """
    Title      test_b()
    Usage
    Function   Transcribe

    Returns    Results of the test
    Argument   None
    """

    def test_b( self ):

        s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
        dna = Seq.Seq(s, IUPAC.unambiguous_dna)
        print 'ready to transcribe'
        rna = Transcribe.unambiguous_transcriber.transcribe( dna )

	return ( self.assert_equals( 'GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU', rna.tostring() ) )

    """
    Title      test_c()
    Usage
    Function   Transcribe

    Returns    Results of the test
    Argument   None
    """

    def test_c( self ):

        s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
        rna = Seq.Seq(s, IUPAC.unambiguous_rna)
        print 'ready to transcribe'
        dna = Transcribe.unambiguous_transcriber.back_transcribe( rna )

	return ( self.assert_equals( 'GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT', dna.tostring() ) )


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

