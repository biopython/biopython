# Copyright 2002 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Based on ideas from
Immunocomputing: a survey
I.Antoniou, S.Gutnikov, V.Ivanov, Yu.Melnikov, A.Tarakanov
12. Forrest S., Perelson A. Aleen L. and Cherukuri R. 
Self-nonself disctimination in a computer. Proc. of IEEE symposium on reseqrch
in security and privacy. Oakland, USA, 1994, 202-212.


Immune system simulation.
Accepts an initial set of sequences to be protected.
Creates a set of randomly scrambled sequences and uses a lazy check to remove
those that trigger on members of the protected set.
The detector for a suspicious sequence checks for a close match to a scrambled sequence.
The detectors start out with equal weights.  When a detector finds a suspicious antigen,
its weight is incremented so its chances of being selected in the future increases.
Intended only for experimentation.
"""

import os
import sys
import string
from urllib import FancyURLopener
from urllib import urlencode

from Bio.RecordFile import RecordFile
from Bio.SGMLExtractor import SGMLExtractorHandle
from Bio.NetCatch import NetCatch
from Bio.NetCatch import ExtractUrls
from Bio.Seq import Seq
from Bio.Align.Generic import Alignment
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Alphabet import DNAAlphabet
from Bio.Alphabet import Gapped
from Bio.SGMLExtractor import SGMLExtractorHandle
from Bio.HotRand import HotRandom
from Bio.HotRand import hex_convert




def match_sequence( first, second, threshold ):
    len_first = len( first )
    len_second = len( second )
    if( len_first > len_second ):
        len_min = len_second
    else:
        len_min = len_first
    if( threshold > len_min ):
        threshold = len_min
    max_match = 0
    match_count = 0
    for j in range( 0, len_min ):
        if( first[ j ] == second[ j ] ):
            match_count = match_count + 1
            if( match_count > max_match ):
                max_match = match_count
        else:
            match_count = 0
    if( max_match >= threshold ):
        return 1
    else:
        return 0

class Lymphocyte:

    def __init__( self, residues ):
        self.residues = residues
        self.may_be_autoimmune = 1
        self.weight = 1


class Immune:
    """
    friendly should be an instance of Align.  It should contain the set of
    protected sequences.
    """

    def __init__( self, friendly, alphabet = [ 'a', 'c', 'g', 't' ], size = 20 ):
        self.hot_random = HotRandom()
        self.friendly = friendly
        self.alphabet = alphabet[:]
        self.lymphocyte_factory( size )

    def select_at_random( self, items ):
        selector = self.hot_random.hot_rand( len( items ) - 1 )
        return selector

    def scramble( self, seq ):
        seq = seq.lower()
        for index in range( 0, len( seq ) ):
            if( seq[ index ] not in self.alphabet ):
                seq[ index ] = seq[ select_at_random( self.alphabet ) ]
        return seq

    def found_antigen( self, detector, mystery_sequence, threshold = 7 ):
        return( match_sequence( detector, mystery_sequence, threshold ) )

    def lazy_auto_immune_check( self, seq ):
        auto_immune = 0
        for candidate in self.friendly.get_all_seqs():
            if( self.found_antigen( seq, candidate.seq.data ) ):
                auto_immune = 1
                break
        return auto_immune


    def compute_accum_weight( self ):
        accum_weight = 0
        for index in range( 0, len( self.lymphocytes ) ):
            lymphocyte = self.lymphocytes[ index ]
            accum_weight = accum_weight + lymphocyte.weight
            lymphocyte.accum_weight = accum_weight
            self.lymphocytes[ index ] = lymphocyte
        self.accum_weight = accum_weight
        return self.accum_weight

    def search_accum_weight( self, t):
        min = 0; max = len( self.lymphocytes ) - 1
        while 1:
            if max < min:
                return min
            m = (min + max) / 2
            if self.lymphocytes[ m ].accum_weight < t:
                min = m + 1
            elif self.lymphocytes[ m ].accum_weight > t:
                max = m - 1
            else:
                return m

    def pick_a_lymphocyte( self ):
        """
        Random selection biased by weight
        """
        weight = self.hot_random.hot_rand( self.accum_weight )
        index = self.search_accum_weight( weight )

    def random_test( self, mystery_sequence ):
        """
        A single test probably won't catch a corrupted sequence.
        Lots of tests are required
        """
        index = self.pick_a_lymphocyte()
        lymphocyte = self.lymphocytes[ index ]
        detector = lymphocyte.residues
        suspicious = self.found_antigen( detector, mystery_sequence )
        if suspicious:
            if( lymphocyte.may_be_autoimmune ):
                auto_immune = self.lazy_auto_immune_check( detector )
            if( auto_immune ):
                self.lymphocytes.remove( index )
                self.create_lymphocyte()
                suspicious = 0
            else:
                lymphocyte.may_be_autoimmune = auto_immune
                lymphocyte.weight = lymphocyte.weight + 1
                self.lymphocytes[ index ] = lymphocyte
                self.compute_accum_weight()
        return suspicious


    def create_lymphocyte( self ):
        lymphocyte = self.scramble( consensus.data )
        self.lymphocytes.append( Lymphocyte( lymphocyte ) )
        self.compute_accum_weight()



    def lymphocyte_factory( self, num_lymphocytes = 20 ):
        self.lymphocytes = []
        summary_info = SummaryInfo( self.friendly )
        consensus = summary_info.dumb_consensus()
        self.consensus = consensus
        for j in range( 0, num_lymphocytes ):
            lymphocyte = self.scramble( consensus.data )
            self.lymphocytes.append( Lymphocyte( lymphocyte ) )
        self.compute_accum_weight()





