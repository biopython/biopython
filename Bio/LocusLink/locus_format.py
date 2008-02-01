


"""Martel based parser to read LocusLink flat files.

This is a huge regular expression for LocusLink,
built using the 'regular expressions on steroids' capabilities of
Martel.

A description of the format can be found in the 'ligand.doc' file
from the Ligand distribution, available from:

 http://www.ncbi.nih.gov/LocusLink


"""

# Martel
from Martel import Str
from Martel import Str1
from Martel import Alt
from Martel import Rep
from Martel import Group
from Martel import ToEol
from Martel import AnyEol
from Martel import Any
from Martel import Word
from Martel import Opt
from Martel import AssertNot

from Martel import RecordReader

import warnings
warnings.warn("Bio.LocusLink was deprecated, as NCBI's LocusLink was superceded by Entrez Gene. If you still need this module, please get in touch with the Biopython developers (biopython-dev@biopython.org) to avoid permanent removal of this module", DeprecationWarning)


# --- First set up some helper constants and functions
INDENT = 12

blank_spaces = Rep(Str1(" "))
point = Str1(".")

white_space = Rep( Any( " 	" ) )
locus_keys = [ \
        'LOCUSID', \
        'LOCUS_CONFIRMED', \
        'LOCUS_TYPE', \
        'ORGANISM', \
        'STATUS', \
        'NM', \
        'NP', \
        'CDD', \
        'PRODUCT', \
        'ASSEMBLY', \
        'CONTIG', \
        'EVID', \
        'XM', \
        'XP', \
        'ACCNUM', \
        'TYPE', \
        'PROT', \
        'OFFICIAL_SYMBOL', \
        'OFFICIAL_GENE_NAME', \
        'PREFERRED_PRODUCT', \
        'ALIAS_SYMBOL', \
        'SUMMARY', \
        'CHR', \
        'STS', \
        'COMP', \
        'ALIAS_PROT', \
        'UNIGENE', \
        'BUTTON', \
        'LINK', \
        'OMIM', \
        'MAP', \
        'MAPLINK', \
        'ECNUM', \
        'PROTOTYPE', \
        'DB_DESCR', \
        'DB_LINK', \
        'PMID', \
        'GRIF', \
        'SUBFUNC', \
        'GO', \
        'EXTANNOT'
        
        ]

accnum_block_keys = [ \
    'ACCNUM', \
    'TYPE', \
    'PROT' \
    ]
phenotype = Str1( 'PHENOTYPE' )
db = Str1( 'DB' )
accnum_block_key = Str( *accnum_block_keys )



valid_locus_key = Str( *locus_keys )
def define_locus_line( entry_tag ):
    
    return(   white_space + \
             Str1(entry_tag ) + \
             white_space + \
	     Str1( ":" ) + \
	     white_space + \
	     ToEol() )

def define_locus_group( entry_name, entry_tag ):   
    return Group( entry_name, \
                  define_locus_line( entry_tag ))

accnum_block = Group( 'accnum_block', \
                      define_locus_line( 'ACCNUM' ) + \
                      define_locus_line( 'TYPE' ) + \
                      Opt(define_locus_line( 'PROT' ) ) ) 

phenotype_block = Group( 'phenotype_block', \
                      define_locus_line( 'PHENOTYPE' ) + \
                      Opt( define_locus_line( 'PHENOTYPE_ID' ) ) )

db_block = Group( 'db_block', \
                      define_locus_line( 'DB_DESCR' ) + \
                      define_locus_line( 'DB_LINK' ) )

begin_record_line = Str1( '>>' ) + ToEol()
locus_line = Group( 'locus_line', \
                 white_space + AssertNot( accnum_block_key ) + AssertNot( phenotype ) + AssertNot( db ) + Word() + white_space + Str1( ':' ) + ToEol() )

locus_record = begin_record_line + Rep( Alt( locus_line, accnum_block, phenotype_block, db_block ) )
#locus_record = Rep( locus_line )
