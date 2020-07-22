# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Several routines used to extract information from FSSP sections.

filter: filters a passed summary section and alignment section according to a numeric
        attribute in the summary section. Returns new summary and alignment sections

For example, to filter in only  those records which have a zscore greater than
4.0 and lesser than 7.5:

new_sum, new_align = filter(sum, align, 'zscore', 4, 7.5)
"""

from Bio import FSSP
import copy
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class FSSPAlign(MultipleSeqAlignment):
    """Provision to do single Multi Sequence Alignment from FSSP files."""

    def _add_numbering_table(self, new_record):
        new_record.annotations["abs2pdb"] = {}
        new_record.annotations["pdb2abs"] = {}


class FSSPMultAlign(dict):
    """Provision to do multiple Multi Sequence Alignment from FSSP files."""

    def __init__(self):
        """Initialize the class."""
        self.abs_res = []
        self.pdb_res = []
        self.data = {}


def mult_align(sum_dict, align_dict):
    """Return multiple alignment instance (MultipleSeqAlignment)."""
    mult_align_dict = {}
    for j in align_dict.abs(1).pos_align_dict:
        mult_align_dict[j] = ""

    for i in range(1, len(align_dict) + 1):
        # loop on positions
        for j in align_dict.abs(i).pos_align_dict:
            # loop within a position
            mult_align_dict[j] += align_dict.abs(i).pos_align_dict[j].aa
    fssp_align = MultipleSeqAlignment([])
    for i in sorted(mult_align_dict):
        fssp_align.append(
            SeqRecord(Seq(mult_align_dict[i]), sum_dict[i].pdb2 + sum_dict[i].chain2)
        )
    return fssp_align


#
# Warning: this function really slows down when filtering large FSSP files.
# The reason is the use of copy.deepcopy() to copy align_dict into
# new_align_dict. I have to figure out something better.
# Took me ~160 seconds for the largest FSSP file (1reqA.fssp)
#


def filter(sum_dict, align_dict, filter_attribute, low_bound, high_bound):
    """Filter a passed summary section and alignment section.

    Filter according to a numeric attribute in the summary section.
    Return new summary and alignment sections.
    """
    new_sum_dict = FSSP.FSSPSumDict()
    new_align_dict = copy.deepcopy(align_dict)
    for prot_num in sum_dict:
        attr_value = getattr(sum_dict[prot_num], filter_attribute)
        if attr_value >= low_bound and attr_value <= high_bound:
            new_sum_dict[prot_num] = sum_dict[prot_num]
    prot_numbers = sorted(new_sum_dict)
    for pos_num in new_align_dict.abs_res_dict:
        new_align_dict.abs(pos_num).pos_align_dict = {}
        for prot_num in prot_numbers:
            new_align_dict.abs(pos_num).pos_align_dict[prot_num] = align_dict.abs(
                pos_num
            ).pos_align_dict[prot_num]
    return new_sum_dict, new_align_dict


def name_filter(sum_dict, align_dict, name_list):
    """Filter summary and alignment blocks for given names only.

    Accepts a list of names. Returns a new Summary block and Alignment block which
    contain the info only for those names passed.
    """
    new_sum_dict = FSSP.FSSPSumDict()
    new_align_dict = copy.deepcopy(align_dict)
    for cur_pdb_name in name_list:
        for prot_num in sum_dict:
            if sum_dict[prot_num].pdb2 + sum_dict[prot_num].chain2 == cur_pdb_name:
                new_sum_dict[prot_num] = sum_dict[prot_num]
    prot_numbers = sorted(new_sum_dict)
    for pos_num in new_align_dict.abs_res_dict:
        new_align_dict.abs(pos_num).pos_align_dict = {}
        for prot_num in prot_numbers:
            new_align_dict.abs(pos_num).pos_align_dict[prot_num] = align_dict.abs(
                pos_num
            ).pos_align_dict[prot_num]
    return new_sum_dict, new_align_dict
