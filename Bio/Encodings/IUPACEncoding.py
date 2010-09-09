"""Properties once used for transcription and translation (DEPRECATED).

This module is deprecated, and is expected to be removed in the next release.
If you use this module, please contact the Biopython developers via the
mailing lists.
"""
#NOTE - Adding a deprecation warning would affect Bio.Alphabet.IUPAC

# Set up the IUPAC alphabet properties


from Bio.PropertyManager import default_manager
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data import IUPACData

from Bio import Transcribe, Translate

set_prop = default_manager.class_property

# weight tables
set_prop[IUPAC.IUPACUnambiguousDNA]["weight_table"] = \
             IUPACData.unambiguous_dna_weights
set_prop[IUPAC.IUPACAmbiguousDNA]["weight_table"] = \
             IUPACData.avg_ambiguous_dna_weights
set_prop[IUPAC.IUPACUnambiguousRNA]["weight_table"] = \
             IUPACData.unambiguous_rna_weights
set_prop[IUPAC.IUPACAmbiguousRNA]["weight_table"] = \
             IUPACData.avg_ambiguous_rna_weights
set_prop[IUPAC.IUPACProtein]["weight_table"] = \
             IUPACData.protein_weights
set_prop[IUPAC.ExtendedIUPACProtein]["weight_table"] = \
             IUPACData.avg_extended_protein_weights

set_prop[IUPAC.IUPACUnambiguousDNA]["weight_range_table"] = \
             IUPACData.unambiguous_dna_weight_ranges
set_prop[IUPAC.IUPACAmbiguousDNA]["weight_range_table"] = \
             IUPACData.ambiguous_dna_weight_ranges
set_prop[IUPAC.IUPACUnambiguousRNA]["weight_range_table"] = \
             IUPACData.unambiguous_rna_weight_ranges
set_prop[IUPAC.IUPACAmbiguousRNA]["weight_range_table"] = \
             IUPACData.ambiguous_rna_weight_ranges
set_prop[IUPAC.IUPACProtein]["weight_range_table"] = \
             IUPACData.protein_weight_ranges
set_prop[IUPAC.ExtendedIUPACProtein]["weight_range_table"] = \
             IUPACData.extended_protein_weight_ranges



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
for name, obj in Translate.unambiguous_dna_by_name.iteritems():
    property = "translator.name." + name
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    set_prop[obj.table.protein_alphabet.__class__][property] = obj

for name, obj in Translate.unambiguous_rna_by_name.iteritems():
    property = "translator.name." + name
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "rna_translator.name." + name
    set_prop[obj.table.protein_alphabet.__class__][property] = obj


for id, obj in Translate.unambiguous_dna_by_id.iteritems():
    property = "translator.id.%d" % id
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    set_prop[obj.table.protein_alphabet.__class__][property] = obj
    if id == 1:
        set_prop[obj.table.nucleotide_alphabet.__class__]["translator"] = obj
        set_prop[obj.table.protein_alphabet.__class__]["translator"] = obj


for id, obj in Translate.unambiguous_rna_by_id.iteritems():
    property = "translator.id.%d" % id
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "rna_translator.id.%d" % id
    set_prop[obj.table.protein_alphabet.__class__][property] = obj
    if id == 1:
        set_prop[obj.table.nucleotide_alphabet.__class__]["translator"] = obj
        set_prop[obj.table.protein_alphabet.__class__]["rna_translator"] = obj

# ambiguous translator objects
for name, obj in Translate.ambiguous_dna_by_name.iteritems():
    property = "translator.name." + name
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "ambiguous_translator.name." + name
    set_prop[obj.table.protein_alphabet.__class__][property] = obj

for name, obj in Translate.ambiguous_rna_by_name.iteritems():
    property = "translator.name." + name
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "ambiguous_rna_translator.name." + name
    set_prop[obj.table.protein_alphabet.__class__][property] = obj


for id, obj in Translate.ambiguous_dna_by_id.iteritems():
    property = "translator.id.%d" % id
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "ambiguous_translator.id.%d" % id
    set_prop[obj.table.protein_alphabet.__class__][property] = obj
    if id == 1:
        set_prop[obj.table.nucleotide_alphabet.__class__]["translator"] = obj
        set_prop[obj.table.protein_alphabet.__class__]["ambiguous_translator"] = obj


for id, obj in Translate.ambiguous_rna_by_id.iteritems():
    property = "translator.id.%d" % id
    set_prop[obj.table.nucleotide_alphabet.__class__][property] = obj
    property = "ambiguous_rna_translator.id.%d" % id
    set_prop[obj.table.protein_alphabet.__class__][property] = obj
    if id == 1:
        set_prop[obj.table.nucleotide_alphabet.__class__]["translator"] = obj
        set_prop[obj.table.protein_alphabet.__class__]["ambiguous_rna_translator"] = obj
