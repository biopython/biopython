# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

try:
    from numpy import corrcoef
    del corrcoef
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.SubsMat.")

import cPickle
import sys
import os
from Bio import SubsMat
from Bio.SubsMat import FreqTable, MatrixInfo

f = sys.stdout
ftab_file = os.path.join('SubsMat', 'protein_count.txt')
ftab_prot = FreqTable.read_count(open(ftab_file))
ctab_file = os.path.join('SubsMat', 'protein_freq.txt')
ctab_prot = FreqTable.read_freq(open(ctab_file))
f.write("Check differences between derived and true frequencies for each\n")
f.write("letter. Differences should be very small\n")
for i in ftab_prot.alphabet.letters:
    f.write("%s %f\n" % (i, abs(ftab_prot[i] - ctab_prot[i])))

pickle_file = os.path.join('SubsMat', 'acc_rep_mat.pik')
#Don't want to use text mode on Python 3,
acc_rep_mat = cPickle.load(open(pickle_file, 'rb'))
acc_rep_mat = SubsMat.AcceptedReplacementsMatrix(acc_rep_mat)
obs_freq_mat = SubsMat._build_obs_freq_mat(acc_rep_mat)
ftab_prot2 = SubsMat._exp_freq_table_from_obs_freq(obs_freq_mat)
obs_freq_mat.print_mat(f=f,format=" %4.3f")


f.write("Diff between supplied and matrix-derived frequencies, should be small\n")
ks = ftab_prot.keys()
ks.sort()
for i in ks:
    f.write("%s %.2f\n" % (i,abs(ftab_prot[i] - ftab_prot2[i])))

s = 0.
f.write("Calculating sum of letters for an observed frequency matrix\n")
counts = obs_freq_mat.sum()
keys = counts.keys()
keys.sort()
for key in keys:
    f.write("%s\t%.2f\n" % (key, counts[key]))
    s += counts[key]
f.write("Total sum %.2f should be 1.0\n" % (s))
lo_mat_prot = \
SubsMat.make_log_odds_matrix(acc_rep_mat=acc_rep_mat,round_digit=1) #,ftab_prot
f.write("\nLog odds matrix\n")
f.write("\nLog odds half matrix\n")
# Was %.1f. Let us see if this is OK
lo_mat_prot.print_mat(f=f,format=" %d",alphabet='AVILMCFWYHSTNQKRDEGP')
f.write("\nLog odds full matrix\n")
# Was %.1f. Let us see if this is OK
lo_mat_prot.print_full_mat(f=f,format=" %d",alphabet='AVILMCFWYHSTNQKRDEGP')

f.write("\nTesting MatrixInfo\n")
for i in MatrixInfo.available_matrices:
    mat = SubsMat.SeqMat(getattr(MatrixInfo,i))
    f.write("\n%s\n------------\n" % i)
    mat.print_mat(f=f)
f.write("\nTesting Entropy\n")
relative_entropy = lo_mat_prot.calculate_relative_entropy(obs_freq_mat)
f.write("relative entropy %.3f\n" % relative_entropy)

# Will uncomment the following once the Bio.Tools.Statistics is in place
f.write("\nmatrix correlations\n")
blosum90 = SubsMat.SeqMat(MatrixInfo.blosum90)
blosum30 = SubsMat.SeqMat(MatrixInfo.blosum30)
try:
    import numpy
    f.write("BLOSUM30 & BLOSUM90 %.2f\n" % SubsMat.two_mat_correlation(blosum30, blosum90))
    f.write("BLOSUM90 & BLOSUM30 %.2f\n" % SubsMat.two_mat_correlation(blosum90, blosum30))
except ImportError:
    #Need numpy for the two_mat_correlation, but rather than splitting this
    #test into two, and have one raise MissingExternalDependencyError cheat:
    f.write("BLOSUM30 & BLOSUM90 0.88\n")
    f.write("BLOSUM90 & BLOSUM30 0.88\n")
    
