# Copyright 2008-2010 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "ace" file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.Sequencing.Ace module which offers more than just accessing
the contig consensus sequences in an ACE file as SeqRecord objects.
"""

from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide, generic_dna, generic_rna, Gapped
from Bio.Sequencing import Ace

__docformat__ = "restructuredtext en"


def AceIterator(handle):
    """Returns SeqRecord objects from an ACE file.

    This uses the Bio.Sequencing.Ace module to do the hard work.  Note that
    by iterating over the file in a single pass, we are forced to ignore any
    WA, CT, RT or WR footer tags.

    Ace files include the base quality for each position, which are taken
    to be PHRED style scores. Just as if you had read in a FASTQ or QUAL file
    using PHRED scores using Bio.SeqIO, these are stored in the SeqRecord's
    letter_annotations dictionary under the "phred_quality" key.

    >>> from Bio import SeqIO
    >>> with open("Ace/consed_sample.ace", "rU") as handle:
    ...     for record in SeqIO.parse(handle, "ace"):
    ...         print("%s %s... %i" % (record.id, record.seq[:10], len(record)))
    ...         print(max(record.letter_annotations["phred_quality"]))
    Contig1 agccccgggc... 1475
    90

    However, ACE files do not include a base quality for any gaps in the
    consensus sequence, and these are represented in Biopython with a quality
    of zero. Using zero is perhaps misleading as there may be very strong
    evidence to support the gap in the consensus. Previous versions of
    Biopython therefore used None instead, but this complicated usage, and
    prevented output of the gapped sequence as FASTQ format.

    >>> from Bio import SeqIO
    >>> with open("Ace/contig1.ace", "rU") as handle:
    ...     for record in SeqIO.parse(handle, "ace"):
    ...         print("%s ...%s..." % (record.id, record.seq[85:95]))
    ...         print(record.letter_annotations["phred_quality"][85:95])
    ...         print(max(record.letter_annotations["phred_quality"]))
    Contig1 ...AGAGG-ATGC...
    [57, 57, 54, 57, 57, 0, 57, 72, 72, 72]
    90
    Contig2 ...GAATTACTAT...
    [68, 68, 68, 68, 68, 68, 68, 68, 68, 68]
    90

    """
    for ace_contig in Ace.parse(handle):
        # Convert the ACE contig record into a SeqRecord...
        consensus_seq_str = ace_contig.sequence
        # Assume its DNA unless there is a U in it,
        if "U" in consensus_seq_str:
            if "T" in consensus_seq_str:
                # Very odd! Error?
                alpha = generic_nucleotide
            else:
                alpha = generic_rna
        else:
            alpha = generic_dna

        if "*" in consensus_seq_str:
            # For consistency with most other file formats, map
            # any * gaps into - gaps.
            assert "-" not in consensus_seq_str
            consensus_seq = Seq(consensus_seq_str.replace("*", "-"),
                                Gapped(alpha, gap_char="-"))
        else:
            consensus_seq = Seq(consensus_seq_str, alpha)

        # TODO? - Base segments (BS lines) which indicates which read
        # phrap has chosen to be the consensus at a particular position.
        # Perhaps as SeqFeature objects?

        # TODO - Supporting reads (RD lines, plus perhaps QA and DS lines)
        # Perhaps as SeqFeature objects?

        seq_record = SeqRecord(consensus_seq,
                               id=ace_contig.name,
                               name=ace_contig.name)

        # Consensus base quality (BQ lines).  Note that any gaps (originally
        # as * characters) in the consensus do not get a quality entry, so
        # we assign a quality of None (zero would be missleading as there may
        # be excelent support for having a gap here).
        quals = []
        i = 0
        for base in consensus_seq:
            if base == "-":
                quals.append(0)
            else:
                quals.append(ace_contig.quality[i])
                i += 1
        assert i == len(ace_contig.quality)
        seq_record.letter_annotations["phred_quality"] = quals

        yield seq_record
    # All done


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
