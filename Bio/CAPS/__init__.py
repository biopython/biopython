# Copyright 2005 by Jonathan Taylor.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Cleaved amplified polymorphic sequence (CAPS) markers.

A CAPS marker is a location a DifferentialCutsite as described below and a
set of primers that can be used to visualize this.  More information can
be found in the paper `Konieczny and Ausubel (1993)`_ (PMID 8106085).

.. _`Konieczny and Ausubel (1993)`: https://doi.org/10.1046/j.1365-313X.1993.04020403.x

"""

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq


class DifferentialCutsite:
    """Differential enzyme cutsite in an alignment.

    A differential cutsite is a location in an alignment where an enzyme cuts
    at least one sequence and also cannot cut at least one other sequence.

    Members:
     - start - Where it lives in the alignment.
     - enzyme - The enzyme that causes this.
     - cuts_in - A list of sequences (as indexes into the alignment) the
       enzyme cuts in.
     - blocked_in - A list of sequences (as indexes into the alignment) the
       enzyme is blocked in.

    """

    def __init__(self, **kwds):
        """Initialize a DifferentialCutsite.

        Each member (as listed in the class description) should be included as a
        keyword.
        """
        self.start = int(kwds["start"])
        self.enzyme = kwds["enzyme"]
        self.cuts_in = kwds["cuts_in"]
        self.blocked_in = kwds["blocked_in"]


class AlignmentHasDifferentLengthsError(Exception):
    """Exception where sequences in alignment have different lengths."""

    pass


class CAPSMap:
    """A map of an alignment showing all possible dcuts.

    Members:
     - alignment - The alignment that is mapped.
     - dcuts - A list of possible CAPS markers in the form of
       DifferentialCutsites.

    """

    def __init__(self, alignment, enzymes=None):
        """Initialize the CAPSMap.

        Required:
         - alignment - The alignment to be mapped.

        Optional:
         - enzymes - List of enzymes to be used to create the map.
           Defaults to an empty list.

        """
        if enzymes is None:
            enzymes = []
        if isinstance(alignment, MultipleSeqAlignment):
            self.sequences = [rec.seq for rec in alignment]
            self.size = len(self.sequences)
            self.length = len(self.sequences[0])
            for seq in self.sequences:
                if len(seq) != self.length:
                    raise AlignmentHasDifferentLengthsError
        else:  # Alignment object
            self.sequences = [Seq(s) for s in alignment]
            self.size, self.length = alignment.shape

        self.alignment = alignment
        self.enzymes = enzymes

        # look for dcuts
        self._digest()

    def _digest_with(self, enzyme):
        cuts = []  # list of lists, one per sequence
        all_seq_cuts = []

        # go through each sequence
        for seq in self.sequences:
            # grab all the cuts in the sequence
            seq_cuts = [cut - enzyme.fst5 for cut in enzyme.search(seq)]
            # maintain a list of all cuts in all sequences
            all_seq_cuts.extend(seq_cuts)
            cuts.append(seq_cuts)

        # we sort the all list and remove duplicates
        all_seq_cuts = sorted(set(all_seq_cuts))

        for cut in all_seq_cuts:
            # test for dcuts

            cuts_in = []
            blocked_in = []

            for i, seq in enumerate(self.sequences):
                if cut in cuts[i]:
                    cuts_in.append(i)
                else:
                    blocked_in.append(i)

            if cuts_in and blocked_in:
                self.dcuts.append(
                    DifferentialCutsite(
                        start=cut, enzyme=enzyme, cuts_in=cuts_in, blocked_in=blocked_in
                    )
                )

    def _digest(self):
        self.dcuts = []

        for enzyme in self.enzymes:
            self._digest_with(enzyme)
