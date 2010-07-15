from Bio import FSSP
import copy
from Bio.Align import Generic
from Bio import Alphabet
import time

class FSSPAlign(Generic.Alignment):
    def _add_numbering_table(self, new_record):
        new_record.annotations['abs2pdb'] = {}
        new_record.annotations['pdb2abs'] = {}

   
class FSSPMultAlign(dict):
    def __init__(self):
        self.abs_res = []
        self.pdb_res = []
        self.data = {}
def mult_align(sum_dict,align_dict):
   """Returns a biopython multiple alignment instance (Bio.Align.Generic)"""
   mult_align_dict = {}
   for j in align_dict.abs(1).pos_align_dict:
      mult_align_dict[j] = ''
   
   for i in range(1,len(align_dict)+1):
      # loop on positions
      for j in align_dict.abs(i).pos_align_dict:
         # loop within a position
         mult_align_dict[j] += align_dict.abs(i).pos_align_dict[j].aa
   seq_order = mult_align_dict.keys()
   seq_order.sort()
   fssp_align = Generic.Alignment(Alphabet.Gapped(
                                  Alphabet.IUPAC.extended_protein))
   for i in seq_order:
      fssp_align.add_sequence(sum_dict[i].pdb2+sum_dict[i].chain2,
                                 mult_align_dict[i])
#        fssp_align._add_numbering_table()
   return fssp_align


# Several routines used to extract information from FSSP sections
# filter:
# filters a passed summary section and alignment section according to a numeric
# attribute in the summary section. Returns new summary and alignment sections
# For example, to filter in only  those records which have a zscore greater than
# 4.0 and lesser than 7.5:
# new_sum, new_align = filter(sum, align, 'zscore', 4, 7.5)
#
# Warning: this function really slows down when filtering large FSSP files.
# The reason is the use of copy.deepcopy() to copy align_dict into
# new_align_dict. I have to figure out something better.
# Took me ~160 seconds for the largest FSSP file (1reqA.fssp)
#

def filter(sum_dict,align_dict,filter_attribute,low_bound, high_bound):
   """filters a passed summary section and alignment section according to a numeric
   attribute in the summary section. Returns new summary and alignment sections"""
   new_sum_dict = FSSP.FSSPSumDict()
   new_align_dict = copy.deepcopy(align_dict)
#   for i in align_dict:
#      new_align_dict[i]  = copy.copy(align_dict[i])
   # new_align_dict = copy.copy(align_dict)
   for prot_num in sum_dict:
      attr_value = getattr(sum_dict[prot_num],filter_attribute)
      if (attr_value >= low_bound and
          attr_value <= high_bound):
         new_sum_dict[prot_num] = sum_dict[prot_num]
   prot_numbers = new_sum_dict.keys()
   prot_numbers.sort()
   for pos_num in new_align_dict.abs_res_dict:
      new_align_dict.abs(pos_num).pos_align_dict = {}
      for prot_num in prot_numbers:
         new_align_dict.abs(pos_num).pos_align_dict[prot_num] = \
                   align_dict.abs(pos_num).pos_align_dict[prot_num]
   return new_sum_dict, new_align_dict

def name_filter(sum_dict, align_dict, name_list):
   """ Accepts a list of names. Returns a new Summary block and Alignment block which
       contain the info only for those names passed."""
   new_sum_dict = FSSP.FSSPSumDict()
   new_align_dict = copy.deepcopy(align_dict)
   for cur_pdb_name in name_list:
      for prot_num in sum_dict:
         if sum_dict[prot_num].pdb2+sum_dict[prot_num].chain2 == cur_pdb_name:
            new_sum_dict[prot_num] = sum_dict[prot_num]
   prot_numbers = new_sum_dict.keys()
   prot_numbers.sort()
   for pos_num in new_align_dict.abs_res_dict:
      new_align_dict.abs(pos_num).pos_align_dict = {}
      for prot_num in prot_numbers:
         new_align_dict.abs(pos_num).pos_align_dict[prot_num] = \
                   align_dict.abs(pos_num).pos_align_dict[prot_num]
   return new_sum_dict, new_align_dict

