# Make sure the translation functions work.
# Start simple - unambiguous DNA to unambiguous protein

from Bio import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC

# First, test the transcription functions

s = "ATA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
rna = dna.transcribe()
assert rna.tostring()=="AUA"

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
rna = dna.transcribe()
assert rna.tostring()=='GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU'

s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
rna = Seq.Seq(s, IUPAC.unambiguous_rna)
dna = rna.back_transcribe()
assert dna.tostring()=='GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT'


# use the standard table

# Do some simple tests first
s = "T"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert  protein.tostring()==""

s = "TC"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert protein.tostring()==""

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert protein.tostring()=='ENSFSLDFL'

s = "GAA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(15, to_stop=True)
assert protein.tostring()=="E"

s = "ATA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate('Vertebrate Mitochondrial', to_stop=True)
assert protein.tostring()=="M"

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate('SGC8', to_stop=True)
assert protein.tostring()=='ENSFSLDFLWNPSPSNDAWDSSY'

# use the standard table

s = "TCAAAAAGGTGCATCTAGATG"
print "Starting with", s
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = dna.translate(to_stop=True)
assert isinstance(protein.alphabet, IUPAC.IUPACProtein)

print len(protein), "ungapped residues translated"

gapped_protein = dna.translate()
assert isinstance(gapped_protein.alphabet, Alphabet.HasStopCodon)
print protein.tostring()

print len(gapped_protein), "residues translated, including gaps"
print gapped_protein.tostring()

# This has "AGG" as a stop codon
p2 = dna.translate(table=2, to_stop=True)
print len(p2), "SGC1 has a stop codon"
print p2.tostring()
p2 = dna.translate(table=2)
print "Actually, there are", p2.count("*"), "stops."
print p2.tostring()

# Make sure I can change the stop character
p2 = dna.translate(table=2, stop_symbol="+")
print "Yep,", p2.count("+"), "stops."
print p2.tostring()


# Some of the same things, with RNA
# (The code is the same, so I'm not doing all of the tests.)
rna = Seq.Seq(s.replace("T", "U"), IUPAC.unambiguous_rna)

print "RNA translation ...",
protein_from_rna = rna.translate(to_stop=True)
assert protein.alphabet is protein_from_rna.alphabet
assert protein.tostring() == protein_from_rna.tostring()
print "works."

print "RNA translation to stop ...",
gapped_protein_from_rna = rna.translate()
assert len(gapped_protein) == len(gapped_protein_from_rna)
assert gapped_protein.tostring() == gapped_protein_from_rna.tostring()
print "works."

# some tests for "by name"
# How about some forward ambiguity?
print "Forward ambiguous"
s = "RATGATTARAATYTA"
#     B  D  *  N  L
dna = Seq.Seq(s, IUPAC.ambiguous_dna)
protein = dna.translate('Vertebrate Mitochondrial')
print protein.tostring()
stop_protein = dna.translate('SGC1', to_stop=True)
print stop_protein.tostring()

# XXX (Backwards with ambiguity code is unfinished!)

