# Copyright 1999 by Katharine Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import types
import string
from types import *
from Bio import File
from Bio import ParserSupport
from Bio import Gobase


tests = [ 'G405967.htm', 'G423945.htm', 'P2661406.htm', 'S4323016.htm' ]
#    'ps00488.txt', 'ps00546.txt' ]
record_parser = Gobase.RecordParser()

def print_list( list ):
    for item in list:
        print "testing %s" % item

def print_sequence_record( data ):
    print 'modlecule type %s ' % data.molecule_type
    print 'is plasmid %s ' % data.is_plasmid
    print 'shape %s' % data.shape
#    submission_date
#    update_date
    print 'Entrez Record %s' % data.entrez_record
    print 'Genbank accession %s' % data.genbank_accession

def print_gene_record( data ):
    print 'Gene class is %s' % data.gene_class
    print 'Plasmid encoded is %s' % data.is_plasmid
    print 'Partial gene %s' % data.is_partial_gene
    print 'Pseudo %s' % data.is_pseudo_gene
    print 'Transpliced gene %s' % data.is_transpliced_gene
    print 'Chloroplast origin %s' % data.chloroplast_origin
    print 'Contains intron %s' % data.contains_intron
    print 'ORF %s' % data.orf
    print 'Included in an intron %s' % data.included_in_intron
#    published_info
    print 'Genbank accession %s' % data.genbank_accession
    print 'Entrez record %s' % data.entrez_record
    print 'Product type %s' % data.product_type
    print 'Product class %s' % data.product_class



def print_protein_record( data ):
    print 'Product class is %s ' % data.product_class
    print 'Gene class is %s ' % data.gene_class
    print 'Partial protein %s ' % data.is_partial_protein
    print 'Plasmid %s' % data.is_plasmid
    print 'Function %s ' % data.function
    print 'Entrez record %s' % data.entrez_record


for test in tests:
    print "testing %s" % test
    datafile = os.path.join( 'Gobase', test )
    src_handle = open( datafile )
    data = record_parser.parse( src_handle )
    segments = string.split( str( data.__class__ ), '.' )
    klass = segments[ len( segments ) - 1 ]
    print klass
    if( klass == 'SequenceRecord' ):
        print_sequence_record( data )
    elif( klass == 'GeneRecord' ):
        print_gene_record( data )
    elif( klass == 'ProteinRecord' ):
        print_protein_record( data )
    print data.species_name
    print data.taxon_division
    print '\n'
