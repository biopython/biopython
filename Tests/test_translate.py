# Make sure the translation functions work.
# Start simple - unambiguous DNA to unambiguous protein

from Bio import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio import Transcribe, Translate

# First, test the transcription functions

s = "ATA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
rna = Transcribe.unambiguous_transcriber.transcribe( dna )
assert rna.tostring()=="AUA"

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
rna = Transcribe.unambiguous_transcriber.transcribe( dna )
assert rna.tostring()=='GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU'

s = "GAAAAUUCAUUUUCUUUGGACUUUCUCUGAAAUCCGAGUCCUAGGAAAGAUGCGUGAGAUUCUUCAUAUU"
rna = Seq.Seq(s, IUPAC.unambiguous_rna)
dna = Transcribe.unambiguous_transcriber.back_transcribe( rna )
assert dna.tostring()=='GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT'


# use the standard table
trans = Translate.unambiguous_dna_by_id[1]

# Do some simple tests first
s = "T"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = trans.translate_to_stop(dna)
assert  protein.tostring()==""

s = "TC"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = trans.translate_to_stop(dna)
assert protein.tostring()==""

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
protein = trans.translate_to_stop(dna)
assert protein.tostring()=='ENSFSLDFL'

s = "GAA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
trans = Translate.unambiguous_dna_by_id[ 15 ]
protein = trans.translate_to_stop(dna)
assert protein.tostring()=="E"

s = "ATA"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
trans = Translate.unambiguous_dna_by_name[ 'Vertebrate Mitochondrial' ]
protein = trans.translate_to_stop(dna)
assert protein.tostring()=="M"

s = "GAAAATTCATTTTCTTTGGACTTTCTCTGAAATCCGAGTCCTAGGAAAGATGCGTGAGATTCTTCATATT"
dna = Seq.Seq(s, IUPAC.unambiguous_dna)
trans = Translate.unambiguous_dna_by_name[ 'SGC8' ]
protein = trans.translate_to_stop(dna)
assert protein.tostring()=='ENSFSLDFLWNPSPSNDAWDSSY'

# use the standard table
trans = Translate.unambiguous_dna_by_id[1]

s = "TCAAAAAGGTGCATCTAGATG"
print "Starting with", s
dna = Seq.Seq(s, IUPAC.unambiguous_dna)


protein = trans.translate_to_stop(dna)
assert isinstance(protein.alphabet, IUPAC.IUPACProtein)

print len(protein), "ungapped residues translated"


gapped_protein = trans.translate(dna)
assert isinstance(gapped_protein.alphabet, Alphabet.HasStopCodon)
print protein.tostring()

print len(gapped_protein), "residues translated, including gaps"
print gapped_protein.tostring()

# This has "AGG" as a stop codon
p2 = Translate.unambiguous_dna_by_id[2].translate_to_stop(dna)
print len(p2), "SGC1 has a stop codon"
print p2.tostring()
p2 = Translate.unambiguous_dna_by_id[2].translate(dna)
print "Actually, there are", p2.data.count("*"), "stops."
print p2.tostring()

# Make sure I can change the stop character
p2 = Translate.unambiguous_dna_by_id[2].translate(dna, "+")
print "Yep,", p2.data.count("+"), "stops."
print p2.tostring()


# back translation is not unique!
back_dna = trans.back_translate(protein)
print back_dna.tostring()
assert len(back_dna) == len(protein) * 3
assert isinstance(back_dna.alphabet, IUPAC.IUPACUnambiguousDNA)

# but forward again better give the same results
# (Note: the alphabets will differ - translate returns a gap encoding)
double_back_protein = trans.translate(back_dna)
assert double_back_protein.data == protein.data

# Try the same trick with stops
back_dna2 = trans.back_translate(gapped_protein)
assert len(back_dna2) == 3*len(gapped_protein)
double_back_protein2 = trans.translate(back_dna2)
assert gapped_protein.data == double_back_protein2.data
print repr(gapped_protein.data), "==", repr(double_back_protein2.data)

# Some of the same things, with RNA
# (The code is the same, so I'm not doing all of the tests.)
rna = Seq.Seq(s.replace("T", "U"), IUPAC.unambiguous_rna)
rna_trans = Translate.unambiguous_rna_by_id[1]

print "RNA translation ...",
protein_from_rna = rna_trans.translate_to_stop(rna)
assert protein.alphabet is protein_from_rna.alphabet
assert protein.data == protein_from_rna.data
print "works."

print "RNA translation to stop ...",
gapped_protein_from_rna = rna_trans.translate(rna)
assert len(gapped_protein) == len(gapped_protein_from_rna)
assert gapped_protein.data == gapped_protein_from_rna.data
print "works."

back_rna = rna_trans.back_translate(protein_from_rna)
assert back_dna.data.replace("T", "U") == back_rna.data

# some tests for "by name"
trans = Translate.unambiguous_dna_by_name[ 'Vertebrate Mitochondrial' ]
trans = Translate.unambiguous_dna_by_name[ 'SGC1' ]

# How about some forward ambiguity?
print "Forward ambiguous"
s = "RATGATTARAATYTA"
#     B  D  *  N  L
dna = Seq.Seq(s, IUPAC.ambiguous_dna)
trans = Translate.ambiguous_dna_by_id[1]
protein = trans.translate(dna)
print protein.tostring()
stop_protein = trans.translate_to_stop(dna)
print stop_protein.tostring()

# XXX (Backwards with ambiguity code is unfinished!)

