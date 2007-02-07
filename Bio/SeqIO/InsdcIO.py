# Copyright 2007 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.GenBank.Scanner import GenBankScanner, EmblScanner
from Bio.Alphabet import generic_protein

# NOTE
# ====
# The "brains" for parsing GenBank and EMBL files (and any
# other flat file variants from the INSDC in future) is in
# Bio.GenBank.Scanner (plus the _FeatureConsumer in Bio.GenBank)
#
# See also
# ========
# International Nucleotide Sequence Database Collaboration
# http://www.insdc.org/
# 
# GenBank
# http://www.ncbi.nlm.nih.gov/Genbank/
#
# EMBL Nucleotide Sequence Database
# http://www.ebi.ac.uk/embl/
#
# DDBJ (DNA Data Bank of Japan)
# http://www.ddbj.nig.ac.jp/

def GenBankIterator(handle) :
    """Breaks up a Genbank file into SeqRecord objects

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    #This calls a generator function:
    return GenBankScanner(debug=0).parse_records(handle)

def EmblIterator(handle) :
    """Breaks up an EMBL file into SeqRecord objects

    Every section from the LOCUS line to the terminating // becomes
    a single SeqRecord with associated annotation and features.
    
    Note that for genomes or chromosomes, there is typically only
    one record."""
    #This calls a generator function:
    return EmblScanner(debug=0).parse_records(handle)

def GenBankCdsFeatureIterator(handle, alphabet=generic_protein) :
    """Breaks up a Genbank file into SeqRecord objects for each CDS feature

    Every section from the LOCUS line to the terminating // can contain
    many CDS features.  These are returned as with the stated amino acid
    translation sequence (if given).
    """
    #This calls a generator function:
    return GenBankScanner(debug=0).parse_cds_features(handle, alphabet)
    
def EmblCdsFeatureIterator(handle, alphabet=generic_protein) :
    """Breaks up a EMBL file into SeqRecord objects for each CDS feature

    Every section from the LOCUS line to the terminating // can contain
    many CDS features.  These are returned as with the stated amino acid
    translation sequence (if given).
    """
    #This calls a generator function:
    return EmblScanner(debug=0).parse_cds_features(handle, alphabet)
