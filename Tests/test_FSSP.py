from Bio import FSSP, Align
from Bio.FSSP import FSSPTools
import sys
import os
import cPickle
import time

print sys.argv
if len(sys.argv) == 1:
   test_file = 'FSSP/1cnv.fssp'
else:
   test_file = sys.argv[1]
f = sys.stdout
f.write("\nRead in %s\n" % os.path.basename(test_file))
head_rec, sum_rec, align_rec = FSSP.read_fssp(open(test_file))
f.write("...1cnv.fssp read\n")
for i in dir(head_rec):
   f.write('head_rec.%s %s\n' % (i, str(getattr(head_rec,i))))
f.write("\nlen(sum_rec) = %d; head_rec.nalign = %d\n" %
        (len(sum_rec), head_rec.nalign))
f.write("The above two numbers should be the same\n")
f.write("\nCreate a multiple alignment instance using Bio.Align\n")
alignment = FSSPTools.mult_align(sum_rec, align_rec)
f.write("...Done\n")
f.write("\nFilter in percent ID's >= 15%\n")
sum_ge_15, align_ge_15 = FSSPTools.filter(sum_rec, align_rec, 'pID', 15,100)

f.write("\nnumber of records filtered in: %d\n" % len(sum_ge_15))
f.write("\nRecords filtered in %s\n" % sum_ge_15.keys())
f.write("\nLet's Pickle this\n")
time_1 = time.time()
cPickle.dump((head_rec, sum_rec, align_rec),open('FSSP/mydump.pik','w'))
f.write("Took me %.2f seconds to pickle\n" % (time.time() - time_1))

f.write("\nFilter by name\n")
name_list = ['2hvm0', '1hvq0', '1nar0', '2ebn0']
f.write("\nname list %s\n" % str(name_list))
sum_newnames, align_newnames = FSSPTools.name_filter(sum_rec, align_rec,
                                                     name_list)
f.write("\n%s\n" % sum_newnames)
f.write("\n%s\n" % align_newnames['0P168'].pos_align_dict)
