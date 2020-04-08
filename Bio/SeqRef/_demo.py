from Bio.SeqRef import PdbRef, GbRef

for r in (
    PdbRef("ABCD", "L"),
    PdbRef("AH5K"),
    GbRef("CE93289.1"),
    GbRef("CE3423234", 3),
    GbRef("CE8392922", "4"),
):
    print(repr(r))
