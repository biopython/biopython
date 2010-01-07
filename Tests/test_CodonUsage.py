from Bio.SeqUtils import CodonUsage
import os
import sys

# first make a CAI object
X = CodonUsage.CodonAdaptationIndex()
# now generate an index from a file
if os.path.exists("./CodonUsage/HighlyExpressedGenes.txt"):
    X.generate_index("./CodonUsage/HighlyExpressedGenes.txt")
elif os.path.exists("./Tests/CodonUsage/HighlyExpressedGenes.txt"):
    X.generate_index("./Tests/CodonUsage/HighlyExpressedGenes.txt")
else:
    print "Cannot find the file HighlyExpressedGene.txt\nMake sure you run the tests from within the Tests folder"
    sys.exit()
# alternatively you could use any predefined dictionary like this:
# from CaiIndices import SharpIndex # you can save your dictionary in this file.
# X.SetCaiIndex(SharpIndex)

print "The current index used:"
X.print_index()

print "-" * 60
print "codon adaptation index for test gene: %.2f" % X.cai_for_gene("ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGA")
