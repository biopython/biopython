# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Make sure the translation functions work.
# Start simple - unambiguous DNA to unambiguous protein

from __future__ import print_function

from Bio import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC

# First, test the transcription functions

s = "ATA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
rna = dna.transcribe()
assert str(rna) == "AUA"

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
rna = dna.transcribe()
assert str(rna) == 'GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU'

s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
rna = Seq.Seq(s, IUPAC.unambiguous_rna)
dna = rna.back_transcribe()
assert str(dna) == 'GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT'


# use the standard table

# Do some simple tests first
s = ""
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert str(protein) == ""

s = "TAA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert str(protein) == ""

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert str(protein) == 'ENSFSLDFL'

s = "GAA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(15, to_stop=True)
assert str(protein) == "E"

s = "ATA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate('Vertebrate Mitochondrial', to_stop=True)
assert str(protein) == "M"

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATAT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate('SGC8', to_stop=True)
assert str(protein) == 'ENSFSLDFLWNPSPSNDAWDSSY'

# use the standard table

s = "TCAAAAAGGTGCATCTAGATG"
print("Starting with %s" % s)
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert isinstance(protein.alphabet, IUPAC.IUPACProtein)

print("%i ungapped residues translated" % len(protein))

gapped_protein = dna.translate()
assert isinstance(gapped_protein.alphabet, Alphabet.HasStopCodon)
print(str(protein))

print("%i residues translated, including gaps" % len(gapped_protein))
print(str(gapped_protein))

# This has "AGG" as a stop codon
p2 = dna.translate(table=2, to_stop=True)
print("%i SGC1 has a stop codon" % len(p2))
print(str(p2))
p2 = dna.translate(table=2)
print("Actually, there are %i stops." % p2.count("*"))
print(str(p2))

# Make sure I can change the stop character
p2 = dna.translate(table=2, stop_symbol="+")
print("Yep, %i stops." % p2.count("+"))
print(str(p2))


# Some of the same things, with RNA
# (The code is the same, so I'm not doing all of the tests.)
rna = Seq.Seq(s.replace("T", "U"), IUPAC.unambiguous_rna)

protein_from_rna = rna.translate(to_stop=True)
assert protein.alphabet is protein_from_rna.alphabet
assert str(protein) == str(protein_from_rna)
print("RNA translation ... works.")

gapped_protein_from_rna = rna.translate()
assert len(gapped_protein) == len(gapped_protein_from_rna)
assert str(gapped_protein) == str(gapped_protein_from_rna)
print("RNA translation to stop ... works.")

# some tests for "by name"
# How about some forward ambiguity?
print("Forward ambiguous")
s = "RATGATTARAATYTA"
#     B  D  *  N  L
dna = Seq.Seq(s, IUPAC.ambiguous_dna)
protein = dna.translate('Vertebrate Mitochondrial')
print(str(protein))
stop_protein = dna.translate('SGC1', to_stop=True)
print(str(stop_protein))

# XXX (Backwards with ambiguity code is unfinished!)
