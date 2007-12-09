# Copyright 2006 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "swiss" (aka SwissProt/UniProt) file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.SwissProt module which offers more than just accessing
the sequences as SeqRecord objects."""

from Bio.SwissProt import SProt
    
#This is a normal function!
def SwissIterator(handle) :
    """Breaks up a Swiss-Prot/UniProt file into SeqRecord objects

    Every section from the ID line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    This parser is for the flat file "swiss" format as used by:
     * Swiss-Prot aka SwissProt
     * TrEMBL
     * UniProtKB aka UniProt Knowledgebase

    It does NOT read their new XML file format.
    http://www.expasy.org/sprot/

    For consistency with BioPerl and EMBOSS we call this the "swiss"
    format.
    """
    return SProt.Iterator(handle, SProt.SequenceParser())
