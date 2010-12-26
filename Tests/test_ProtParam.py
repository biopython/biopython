from Bio.SeqUtils import ProtParam, ProtParamData

def PrintDictionary(MyDict):
    for i in sorted(MyDict):
        print "%s\t%.2f" %(i, MyDict[i])
    print ""

X = ProtParam.ProteinAnalysis("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV")

print "Amino acid\tCount"
PrintDictionary(X.count_amino_acids())

print "Amino acid\tFraction"
PrintDictionary(X.get_amino_acids_percent())

print "Molecular weight of test protein: %.2f" % X.molecular_weight()
print "Aromaticity of test protein: %.2f" % X.aromaticity()
print "Instability index of test protein: %.2f" % X.instability_index()
print "length of flexibility list: %i" % len(X.flexibility())
print "The isoelectric point of the test protein is: %.2f" \
      % X.isoelectric_point()
Helix, Turn, Sheet = X.secondary_structure_fraction()
print "Fraction of amino acids in Helix: %.2f" % Helix
print "Fraction of amino acids in Turn: %.2f" % Turn
print "Fraction of amino acids in Sheet: %.2f" % Sheet
print "\nKyte and Doolittle protein scale:"
expected = [-0.0783,+0.0358,+0.1258,+0.6950,+0.8775,+0.8350,+0.2925,+0.3383,
            -0.1733,-0.4142,-0.5292,-0.6108,-0.8308,-0.8100,-0.8208,-1.0283,
            -1.6300,-1.8233,-2.4267,-2.2292,-1.7817,-1.4742,-0.7467,-0.1608,
            +0.1108,+0.2142,+0.1792,-0.1217,-0.4808,-0.4333,-0.5167,-0.2833,
            +0.3758,+0.7225,+0.4958,+0.6033,+0.5625,+0.3108,-0.2408,-0.0575,
            -0.3717,-0.7800,-1.1242,-1.4083,-1.7550,-2.2642,-2.8575,-2.9175,
            -2.5358,-2.5325,-1.8142,-1.4667,-0.6058,-0.4483,+0.1300,+0.1225,
            +0.2825,+0.1650,+0.3317,-0.2000,+0.2683,+0.1233,+0.4092,+0.1392,
            +0.4192,+0.2758,-0.2350,-0.5750,-0.5983,-1.2067,-1.3867,-1.3583,
            -0.8708,-0.5383,-0.3675,+0.0667,+0.0825,-0.0150,+0.1817,+0.4692,
            +0.3017,+0.3800,+0.4825,+0.4675,+0.1575,-0.1783,-0.5175,-1.2017,
            -1.7033,-1.5500,-1.2375,-0.8500,-0.0583,+0.3125,+0.4242,+0.7133,
            +0.5633,+0.0483,-0.7167,-1.3158,-1.9217,-2.5033,-2.4117,-2.2483,
            -2.3758,-2.0633,-1.8900,-1.8667,-1.9292,-1.8625,-2.0050,-2.2708,
            -2.4050,-2.3508,-2.1758,-1.5533,-1.0350,-0.1983,-0.0233,+0.1800,
            +0.0317,-0.0917,-0.6375,-0.9650,-1.4500,-1.6008,-1.7558,-1.5450,
            -1.7900,-1.8133,-2.0125,-2.1383,-2.3142,-2.1525,-2.1425,-1.9733,
            -1.4742,-0.8083,-0.2100,+0.8067,+1.3092,+1.8367,+2.0283,+2.3558]
for i,e in zip(X.protein_scale(ProtParamData.kd, 9, 0.4), expected):
    assert abs(i-e)<0.01
print "ok"
print "\nGRAVY:"
print "%0.4f" % X.gravy()

