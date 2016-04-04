# Copyright 2001 by Iddo Friedberg.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import FSSP
from Bio.FSSP import FSSPTools
import sys
import os
# import pickle

test_file = os.path.join('FSSP', '1cnv.fssp')
f = sys.stdout
f.write("\nRead in %s\n" % os.path.basename(test_file))
handle = open(test_file)
head_rec, sum_rec, align_rec = FSSP.read_fssp(handle)
handle.close()
f.write("...1cnv.fssp read\n")
for i in ["author", "compnd", "database", "header", "nalign",
          "pdbid", "seqlength", "source"]:
    f.write('head_rec.%s %s\n' % (i, str(getattr(head_rec, i))))
f.write("\nlen(sum_rec) = %d; head_rec.nalign = %d\n" %
        (len(sum_rec), head_rec.nalign))
f.write("The above two numbers should be the same\n")
f.write("\nCreate a multiple alignment instance using Bio.Align\n")
alignment = FSSPTools.mult_align(sum_rec, align_rec)
f.write("...Done\n")
# Percent ID filtering takes too long.. remove from test.

# f.write("\nFilter in percent ID's >= 15%\n")
# sum_ge_15, align_ge_15 = FSSPTools.filter(sum_rec, align_rec, 'pID', 15,100)

# f.write("\nnumber of records filtered in: %d\n" % len(sum_ge_15))
# k = sorted(sum_ge_15)
# f.write("\nRecords filtered in %s\n" % k)
# Pickling takes too long.. remove from test.
# f.write("\nLet's Pickle this\n")
# dump_file = os.path.join('FSSP', 'mydump.pik')
# pickle.dump((head_rec, sum_rec, align_rec),open(dump_file, 'w'))

f.write("\nFilter by name\n")
name_list = ['2hvm0', '1hvq0', '1nar0', '2ebn0']
f.write("\nname list %s\n" % str(name_list))
sum_newnames, align_newnames = FSSPTools.name_filter(sum_rec, align_rec,
                                                     name_list)
for key in sorted(sum_newnames):
    f.write("%s : %s\n" % (key, sum_newnames[key]))

new_dict = align_newnames['0P168'].pos_align_dict
for key in sorted(new_dict):
    f.write("%s : %s\n" % (key, new_dict[key]))
