from types import ListType, TupleType
# A superclass for reading [f]ixed-column type [f]lat-[f]ile records. (e.g.
class fff_rec:
   def __init__(self,inrec=''):
      self.data = inrec
   def __repr__(self):
      return str(self.data)
   __str__ = __repr__

   def __getslice__(self,low,high):
      return self.data[low:high]

   def __len__(self):
      return len(self.data)
   def __getitem__(self,index):
      if ((type(index) is TupleType or type(index) is ListType) and 
         len(index) == 2):
         return self.data[index[0]:index[1]]
      else:
         return self.data[index]

# Definition of the align section in a FSSP file
class align:
   abs_res_num = (0,4)
   pdb_res_num = (4,9)
   chain_id = 10
   res_name = 12
   ss1 = 15
   turn3 = 17
   turn4 = 18
   turn5 = (20,22)
   acc = (34,37)
   start_aa_list = 42

