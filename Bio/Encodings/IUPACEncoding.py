# Set up the IUPAC alphabet properties


from Bio.PropertyManager import default_manager
from Bio import Alphabet
from Bio.Alphabet import IUPAC

from Bio.Tools import Transcribe, Translate

set_prop = default_manager.class_property


#  transcriber objects

set_prop[Alphabet.DNAAlphabet]["transcriber"] = \
             Transcribe.generic_transcriber

set_prop[IUPAC.IUPACAmbiguousDNA]["transcriber"] = \
             Transcribe.ambiguous_transcriber

set_prop[IUPAC.IUPACUnambiguousDNA]["transcriber"] = \
             Transcribe.unambiguous_transcriber


set_prop[Alphabet.RNAAlphabet]["transcriber"] = \
             Transcribe.generic_transcriber

set_prop[IUPAC.IUPACAmbiguousRNA]["transcriber"] = \
             Transcribe.ambiguous_transcriber

set_prop[IUPAC.IUPACUnambiguousRNA]["transcriber"] = \
             Transcribe.unambiguous_transcriber


# translator objects
for name, obj in Translate.unambiguous_dna_by_name.items():
    property = "translator.name." + name
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    set_prop[obj.table.protein_alphabet.__class__][property] = obj

for name, obj in Translate.unambiguous_rna_by_name.items():
    property = "translator.name." + name
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "rna_translator.name." + name
    set_prop[obj.table.protein_alphabet.__class__][property] = obj


for id, obj in Translate.unambiguous_dna_by_id.items():
    property = "translator.id.%d" % id
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    set_prop[obj.table.protein_alphabet.__class__][property] = obj
    if id == 1:
        set_prop[obj.table.nucleotide_alphabet.__class__]["translator"] = obj
        set_prop[obj.table.protein_alphabet.__class__]["translator"] = obj


for id, obj in Translate.unambiguous_rna_by_id.items():
    property = "translator.id.%d" % id
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "rna_translator.id.%d" % id
    set_prop[obj.table.protein_alphabet.__class__][property] = obj
    if id == 1:
        set_prop[obj.table.nucleotide_alphabet.__class__]["translator"] = obj
        set_prop[obj.table.protein_alphabet.__class__]["rna_translator"] = obj
