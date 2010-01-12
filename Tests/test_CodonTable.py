# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Data.CodonTable import *

#Check the extension of stop codons to include well defined ambiguous ones
assert list_ambiguous_codons(['TGA', 'TAA'],IUPACData.ambiguous_dna_values) == ['TGA', 'TAA', 'TRA']
assert list_ambiguous_codons(['TAG', 'TGA'],IUPACData.ambiguous_dna_values) == ['TAG', 'TGA']
assert list_ambiguous_codons(['TAG', 'TAA'],IUPACData.ambiguous_dna_values) == ['TAG', 'TAA', 'TAR']
assert list_ambiguous_codons(['UAG', 'UAA'],IUPACData.ambiguous_rna_values) == ['UAG', 'UAA', 'UAR']
assert list_ambiguous_codons(['TGA', 'TAA', 'TAG'],IUPACData.ambiguous_dna_values) == ['TGA', 'TAA', 'TAG', 'TAR', 'TRA']


#Basic sanity test,
for n in ambiguous_generic_by_id.keys():
    assert ambiguous_rna_by_id[n].forward_table["GUU"] == "V"
    assert ambiguous_rna_by_id[n].forward_table["GUN"] == "V"
    if n != 23 :
        assert ambiguous_rna_by_id[n].forward_table["UUN"] == "X" #F or L

    assert ambiguous_dna_by_id[n].forward_table["GTT"] == "V"
    if n != 23 :
        assert ambiguous_dna_by_id[n].forward_table["TTN"] == "X" #F or L
    assert ambiguous_dna_by_id[n].forward_table["GTN"] == "V"

    if n != 23 :
        assert ambiguous_generic_by_id[n].forward_table.get("TTN") == "X"
    assert ambiguous_generic_by_id[n].forward_table["ACN"] == "T"
    assert ambiguous_generic_by_id[n].forward_table["GUU"] == "V"
    assert ambiguous_generic_by_id[n].forward_table["GUN"] == "V"
    if n != 23 :
        assert ambiguous_generic_by_id[n].forward_table["UUN"] == "X" #F or L
    assert ambiguous_generic_by_id[n].forward_table["GTT"] == "V"
    if n != 23 :
        assert ambiguous_generic_by_id[n].forward_table["TTN"] == "X" #F or L
    assert ambiguous_generic_by_id[n].forward_table["GTN"] == "V"
    #And finally something evil, an RNA-DNA mixture:
    if n != 23 :
        assert ambiguous_generic_by_id[n].forward_table["UTN"] == "X" #F or L
    assert ambiguous_generic_by_id[n].forward_table["UTU"] == "F"

    #R = A or G, so URR = UAA or UGA / TRA = TAA or TGA = stop codons
    if "UAA" in unambiguous_rna_by_id[n].stop_codons \
    and "UGA" in unambiguous_rna_by_id[n].stop_codons:
        try:
            print ambiguous_dna_by_id[n].forward_table["TRA"]
            assert False, "Should be a stop only"
        except KeyError:
            pass
        try:
            print ambiguous_rna_by_id[n].forward_table["URA"]
            assert False, "Should be a stop only"
        except KeyError:
            pass
        try:
            print ambiguous_generic_by_id[n].forward_table["URA"]
            assert False, "Should be a stop only"
        except KeyError:
            pass
        assert "URA" in ambiguous_generic_by_id[n].stop_codons
        assert "URA" in ambiguous_rna_by_id[n].stop_codons
        assert "TRA" in ambiguous_generic_by_id[n].stop_codons
        assert "TRA" in ambiguous_dna_by_id[n].stop_codons

    if "UAG" in unambiguous_rna_by_id[n].stop_codons \
    and "UAA" in unambiguous_rna_by_id[n].stop_codons \
    and "UGA" in unambiguous_rna_by_id[n].stop_codons:
        try:
            print ambiguous_dna_by_id[n].forward_table["TAR"]
            assert False, "Should be a stop only"
        except KeyError:
            pass
        try:
            print ambiguous_rna_by_id[n].forward_table["UAR"]
            assert False, "Should be a stop only"
        except KeyError:
            pass
        try:
            print ambiguous_generic_by_id[n].forward_table["UAR"]
            assert False, "Should be a stop only"
        except KeyError:
            pass
        try:
            print ambiguous_generic_by_id[n].forward_table["URR"]
            assert False, "Should be a stop OR an amino"
        except TranslationError:
            pass
        assert "UAR" in ambiguous_generic_by_id[n].stop_codons
        assert "UAR" in ambiguous_rna_by_id[n].stop_codons
        assert "TAR" in ambiguous_generic_by_id[n].stop_codons
        assert "TAR" in ambiguous_dna_by_id[n].stop_codons
        assert "URA" in ambiguous_generic_by_id[n].stop_codons
        assert "URA" in ambiguous_rna_by_id[n].stop_codons
        assert "TRA" in ambiguous_generic_by_id[n].stop_codons
        assert "TRA" in ambiguous_dna_by_id[n].stop_codons

    if "UUG" in unambiguous_rna_by_id[n].start_codons \
    and "CUG" in unambiguous_rna_by_id[n].start_codons \
    and "AUG" in unambiguous_rna_by_id[n].start_codons \
    and "UUG" not in unambiguous_rna_by_id[n].start_codons:
        assert "NUG" not in ambiguous_dna_by_id[n].start_codons
        assert "RUG" not in ambiguous_dna_by_id[n].start_codons
        assert "WUG" not in ambiguous_dna_by_id[n].start_codons
        assert "KUG" not in ambiguous_dna_by_id[n].start_codons
        assert "SUG" not in ambiguous_dna_by_id[n].start_codons
        assert "DUG" not in ambiguous_dna_by_id[n].start_codons
del n

#Table 2 Vertebrate Mitochondrial has
#TAA and TAG -> TAR, plus AGA and AGG -> AGR
assert "AGR" in ambiguous_dna_by_id[2].stop_codons
assert "TAR" in ambiguous_dna_by_id[2].stop_codons
assert "AGR" in ambiguous_rna_by_id[2].stop_codons
assert "UAR" in ambiguous_rna_by_id[2].stop_codons
assert "AGR" in ambiguous_generic_by_id[2].stop_codons
assert "UAR" in ambiguous_generic_by_id[2].stop_codons
assert "TAR" in ambiguous_generic_by_id[2].stop_codons
assert ambiguous_generic_by_id[1].stop_codons == ambiguous_generic_by_name["Standard"].stop_codons
assert ambiguous_generic_by_id[4].stop_codons == ambiguous_generic_by_name["SGC3"].stop_codons
assert ambiguous_generic_by_id[15].stop_codons == ambiguous_generic_by_name['Blepharisma Macronuclear'].stop_codons

print "Done"
