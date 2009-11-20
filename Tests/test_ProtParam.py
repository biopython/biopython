from Bio.SeqUtils import ProtParam, ProtParamData

def PrintDictionary(MyDict):
    for i in sorted(MyDict.keys()):
        print "%s\t%.2f" %(i, MyDict[i])
    print ""

X = ProtParam.ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV")
print "Amino acid\tCount\n", PrintDictionary(X.count_amino_acids())
print "Amino acid\tFraction\n", PrintDictionary(X.get_amino_acids_percent())
print "Molecular weight of test protein:", X.molecular_weight()
print "Aromaticity of test protein: %.2f" % X.aromaticity()
print "Instability index of test protein: %.2f" % X.instability_index()
print "length of flexibility list:", len(X.flexibility())
print "The isoelectric point of the test protein is: %.2f" \
      % X.isoelectric_point()
Helix, Turn, Sheet = X.secondary_structure_fraction()
print "Fraction of amino acids in Helix: %.2f" % Helix
print "Fraction of amino acids in Turn: %.2f" % Turn
print "Fraction of amino acids in Sheet: %.2f" % Sheet
print "\nKyte and Doolittle protein scale:"
for i in X.protein_scale(ProtParamData.kd, 9, 0.4):
    print "% 0.1f" %i
print "\nGRAVY:"
print "%0.4f" % X.gravy()
