# Copyright 2001 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Martel based parser to read Kabat formatted files.

This is a huge regular regular expression for Kabat, built using
the 'regular expressiona on steroids' capabilities of Martel.

http://immuno.bme.nwu.edu/

Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

"""
# standard library
#http://immuno.bme.nwu.edu/seqhunt.html
import string

# Martel
import Martel
from Martel import RecordReader

# - useful constants for dealing with the blank space in GenBank documents
# this is useful since blank space can be significant in GenBank flat files.
INDENT = 12
FEATURE_KEY_INDENT = 5
FEATURE_QUALIFIER_INDENT = 21

# --- first set up some helper constants and functions
blank_space = Martel.MaxRepeat(Martel.Str(" "), 1, 80 )
date = Martel.Group("date",
                    Martel.Re("[-\w]+"))
amino_3_letter_codes = [ 'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'PHE', 'GLY',\
                'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', \
                'ARG', 'SER','THR', 'VAL', 'TRP', 'TYR', '---' ]
# from http://www.lmb.uni-muenchen.de/groups/bioinformatics/02/ch_02.html

amino_1_letter_codes = 'ACDEFGHIKLMNPQRSTVWY-'

nucleotides = 'gcat-'

amino_alts = map( Martel.Str, amino_3_letter_codes )

codon = Martel.Group( "codon", Martel.MaxRepeat( Martel.Any( nucleotides ), 1, 3 ) )
amino_3_letter_code = Martel.Group( "amino_3_letter_code", \
    reduce( Martel.Alt, amino_alts ) )
amino_1_letter_code = Martel.Group( "amino_1_letter_code", \
    Martel.Any( amino_1_letter_codes ) )
residue = Martel.Group( "residue",
                           blank_space +
                           codon +
                           blank_space +
                           amino_3_letter_code +
                           blank_space +
                           amino_1_letter_code )


kabatid = Martel.Group("kabatid",
                    Martel.Rep1(Martel.Integer()))
pubmed_num = Martel.Group("pubmed_num",
                    Martel.Rep1(Martel.Integer()))
residue_num = Martel.Group("residue_num",
                    Martel.Rep1(Martel.Integer()))
kabat_num = Martel.Group("kabat_num",
                    Martel.Rep1(Martel.Integer()))
id_line = Martel.Group("id_line",
                          Martel.Str("KADBID") +
                          blank_space +
                          kabatid +
                          Martel.ToEol() )
creation_date_line = Martel.Group( "creation_date_line", \
                              Martel.Str( "CREATD" ) +
                              blank_space +
                              Martel.ToEol("creation_date") )
last_mod_date_line = Martel.Group( "last_mod_date_line", \
                              Martel.Str( "LMODIF" ) +
                              blank_space +
                              Martel.ToEol("last_mod_date") )
definition_line = Martel.Group( "definition_line", \
                              Martel.Str( "DEFINI" ) +
                              blank_space +
                              Martel.ToEol( "definition" ) )
species_line = Martel.Group( "species_line", \
                              Martel.Str( "SPECIE" ) +
                              blank_space +
                              Martel.ToEol( "species" ) )
amino_acid_sequence_name_line = Martel.Group( "amino_acid_sequence_name_line",
                                Martel.Str( "AANAME" ) +
                                blank_space +
                                Martel.ToEol( "amino_acid_sequence_name" ) )

nucleotide_sequence_name_line = Martel.Group( "nucleotide_sequence_name_line",
                                Martel.Str( "NNNAME" ) +
                                blank_space +
                                Martel.ToEol( "nucleotide_sequence_name" ) )
amino_acid_ref_author_line = Martel.Group( "amino_acid_ref_author_line",
                                Martel.Str( "AAREFA" ) +
                                Martel.ToEol( "amino_acid_ref_author" ) )
amino_acid_ref_journal_line = Martel.Group( "amino_acid_ref_journal_line",
                                Martel.Str( "AAREFJ" ) +
                                Martel.ToEol( "amino_acid_ref_journal" ) )
amino_acid_ref_pubmed_line = Martel.Group( "amino_acid_ref_pubmed_line",
                                Martel.Str( "AAPUBM" ) +
                                Martel.ToEol( "amino_acid_ref_pubmed" ) )
nucleotide_ref_author_line = Martel.Group( "nucleotide_ref_author_line",
                                Martel.Str( "NNREFA" ) +
                                Martel.ToEol( "nucleotide_ref_author" ) )
nucleotide_ref_journal_line = Martel.Group( "nucleotide_ref_journal_line",
                                Martel.Str( "NNREFJ" ) +
                                Martel.ToEol( "nucleotide_ref_journal" ) )
nucleotide_ref_pubmed_line = Martel.Group( "nucleotide_ref_pubmed_line",
                                Martel.Str( "NNPUBM" ) +
                                Martel.ToEol("nucleotide_ref_pubmed") )
annotation_key = Martel.Group( "annotation_key", Martel.RepN( Martel.Re( "\w" ), 4 ) )
annotation_val = Martel.Group( "annotation_val", Martel.ToEol() )
annotation_item = Martel.Group( "annotation_item",
                                annotation_key +
                                blank_space +
                                annotation_val )

annotation_line = Martel.Group( "annotation_line",
                                Martel.Str( "ANNOTA" ) +
                                blank_space + annotation_item )
residue_line = Martel.Group( "residue_line",
                             Martel.Str( "SEQTPA" ) +
                             blank_space +
                             residue_num +
                             blank_space +
                             kabat_num +
                             Martel.Opt( residue ) +
                             Martel.ToEol( ) )
kabat_record_end_line =      Martel.Group( "record_end",
                             Martel.Str( "RECEND" ) +
                             Martel.ToEol() +
                             Martel.ToEol() +
                             Martel.ToEol())
amino_acid_ref = Martel.Group( "amino_acid_ref",
                               amino_acid_ref_author_line +
                               amino_acid_ref_journal_line +
                               amino_acid_ref_pubmed_line )
nucleotide_ref = Martel.Group( "nucleotide_ref",
                               nucleotide_ref_author_line +
                               nucleotide_ref_journal_line +
                               nucleotide_ref_pubmed_line )
residue_lines = Martel.Group( "residue_lines", Martel.Rep( residue_line ) )
annotation = Martel.Group( "annotation", Martel.Rep( annotation_line ) )
refs = Martel.Group( "refs", Martel.Rep1( amino_acid_ref | nucleotide_ref ) )
kabat_record = id_line + creation_date_line + last_mod_date_line + \
               definition_line + species_line + amino_acid_sequence_name_line + \
               nucleotide_sequence_name_line + refs + annotation +  \
               residue_lines + kabat_record_end_line








