# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.SeqFeature import Reference
from Bio.Seq import MutableSeq
"""Hold GenBank data in a straightforward format.

classes:
o Record - All of the information in a GenBank record.
o KabatReference - hold reference data for a record.
o Annotation - Hold the information in a Feature Table.
"""

class KabatReference( Reference ):

    def __init__( self ):
        Reference.__init__( self )

    def print_kabat_reference( self ):
        print self.authors
        print self.journal
        print self.pubmed_id

class Record:
    """Hold Kabat information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GenBank data.

    Attributes:
    kabatid - id of the sequence in the Kabat database
    creation_date
    date_last_mod - date of last modification
    definition - definition of sequence
    species
    nucleic_acid_refs - references to nucleic acid sequence
    amino_acid_refs - references to amino acid sequence
    annotation
    nucleotide_sequence_name
    amino_acid_sequence_name
    nucleotide_sequence
    amino_acid_sequence
    """
    def __init__(self):
        self.kabatid = ''
        self.creation_date = ''
        self.date_last_mod = ''
        self.definition = ''
        self.species = ''
        self.nucleotide_sequence_name = ''
        self.amino_acid_sequence_name = ''
        self.nucleotide_refs = []
        self.amino_acid_refs = []
        self.annotation = {}
        self.nucleotide_sequence = MutableSeq( "" )
        self.amino_acid_sequence = MutableSeq( "" )

    def print_kabat( self ):
        print 'Kabat id: %s' % self.kabatid
        print 'Creation date: %s' % self.creation_date
        print 'Last modification: %s' % self.date_last_mod
        print self.definition
        print 'Species; %s' % self.species
        print self.nucleotide_sequence_name
        print self.amino_acid_sequence_name
        for ref in self.amino_acid_refs:
            ref.print_kabat_reference()
        for ref in self.nucleotide_refs:
            ref.print_kabat_reference()
        for key in self.annotation.keys():
            val = self.annotation[ key ]
            print '%s: %s' % ( key, val )
        dna_seq = self.nucleotide_sequence.tostring()
        print_seq( dna_seq )
        amino_seq = self.amino_acid_sequence.tostring()
        print_seq( amino_seq )

def print_seq( seq ):
    print ""
    for j in range( 0, len( seq ), 80 ):
        print seq[ j: j + 80 ]


