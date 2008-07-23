# Copyright 2008 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "ace" file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.Sequencing.Ace module which offers more than just accessing
the contig consensus sequences in an ACE file as SeqRecord objects."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide, generic_dna, generic_rna, Gapped
from Bio.Sequencing import Ace
    
#This is a generator function!
def AceIterator(handle) :
    """Returns SeqRecord objects from an ACE file.

    This uses the Bio.Sequencing.Ace module to do the hard work.  Note that
    by iterating over the file in a single pass, we are forced to ignore any
    WA, CT, RT or WR footer tags."""

    for ace_contig in Ace.parse(handle) :
        #Convert the ACE contig record into a SeqRecord...
        consensus_seq_str = ace_contig.sequence
        #Assume its DNA unless there is a U in it,
        if "U" in consensus_seq_str :
            if "T" in consensus_seq_str :
                #Very odd! Error?
                alpha = generic_ncleotide
            else :
                alpha = generic_rna
        else :
            alpha = generic_dna
            
        if "*" in consensus_seq_str :
            #For consistency with most other file formats, map
            #any * gaps into 0 gaps.
            assert "-" not in consensus_seq_str
            consensus_seq = Seq(consensus_seq_str.replace("*","-"),
                                Gapped(alpha, gap_char="-"))
        else :
            consensus_seq = Seq(consensus_seq_str, alpha)

        #TODO - Consensus base quality (BQ lines).  Note that any gaps
        #(* character) in the consensus does not get a quality entry.
        #This really needs Biopython support for per-letter-annotation.

        #TODO? - Base segments (BS lines) which indicates which read
        #phrap has chosen to be the consensus at a particular position.
        #Perhaps as SeqFeature objects?

        #TODO - Supporting reads (RD lines, plus perhaps QA and DS lines)
        #Perhaps as SeqFeature objects?
            
        seq_record = SeqRecord(consensus_seq,
                               id = ace_contig.name,
                               name = ace_contig.name)
        yield seq_record 
    #All done
