# Information about the IUPAC alphabets

protein_letters = "ACDEFGHIKLMNPQRSTVWY"
extended_protein_letters = "ACDEFGHIKLMNPQRSTVWYBXZ"
ambiguous_dna_letters = "GATCRYWSMKHBVDN"
unambiguous_dna_letters = "GATC"
ambiguous_rna_letters = "GAUCRYWSMKHBVDN"
unambiguous_rna_letters = "GAUC"

#   B == 5-bromouridine
#   D == 5,6-dihydrouridine
#   S == thiouridine
#   W == wyosine
extended_dna_letters = "GATCBDSW"

# are there extended forms?
#extended_rna_letters = "GAUCBDSW"

ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }
ambiguous_rna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "U": "U",
    "M": "AC",
    "R": "AG",
    "W": "AU",
    "S": "CG",
    "Y": "CU",
    "K": "GU",
    "V": "ACG",
    "H": "ACU",
    "D": "AGU",
    "B": "CGU",
    "X": "GAUC",
    "N": "GAUC",
    }

ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }


def _make_ranges(dict):
    d = {}
    for key, value in dict.items():
        d[key] = (value, value)
    return d

# From bioperl's SeqStats.pm
unambiguous_dna_weights = {
    "A": 347.,
    "C": 323.,
    "G": 363.,
    "T": 322.,
    }
unambiguous_dna_weight_ranges = _make_ranges(unambiguous_dna_weights)

unambiguous_rna_weights = {
    "A": unambiguous_dna_weights["A"] + 16.,  # 16 for the oxygen
    "C": unambiguous_dna_weights["C"] + 16.,
    "G": unambiguous_dna_weights["G"] + 16.,
    "U": 340.,
}
unambiguous_rna_weight_ranges = _make_ranges(unambiguous_rna_weights)

def _make_ambiguous_ranges(dict, weight_table):
    range_d = {}
    avg_d = {}
    for letter, values in dict.items():
        weights = map(weight_table.get, values)
        range_d[letter] = (min(weights), max(weights))
        total_w = 0.0
        for w in weights:
            total_w = total_w + w
        avg_d[letter] = total_w / len(weights)
    return range_d, avg_d

ambiguous_dna_weight_ranges, avg_ambiguous_dna_weights = \
               _make_ambiguous_ranges(ambiguous_dna_values,
                                      unambiguous_dna_weights)

ambiguous_rna_weight_ranges, avg_ambiguous_rna_weights = \
               _make_ambiguous_ranges(ambiguous_rna_values,
                                      unambiguous_rna_weights)

protein_weights = {
    "A": 89.,
    "C": 121.,
    "D": 133.,
    "E": 147.,
    "F": 165.,
    "G": 75.,
    "H": 155.,
    "I": 131.,
    "K": 146.,
    "L": 131.,
    "M": 149.,
    "N": 132.,
    "P": 115.,
    "Q": 146.,
    "R": 174.,
    "S": 105.,
    "T": 119.,
    "V": 117.,
    "W": 204.,
    "Y": 181.,
}
extended_protein_values = {
    "A": "A",
    "B": "ND",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "V": "V",
    "W": "W",
    "X": "ACDEFGHIKLMNPQRSTVWY",
    "Y": "Y",
    "Z": "QE",
}
    
protein_weight_ranges = _make_ranges(protein_weights)

extended_protein_weight_ranges, avg_extended_protein_weights = \
               _make_ambiguous_ranges(extended_protein_values,
                                      protein_weights)



