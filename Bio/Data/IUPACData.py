# Information about the IUPAC alphabets

protein_letters = "ACDEFGHIKLMNPQRSTVWY"
extended_protein_letters = "ACDEFGHIKLMNPQRSTVWYBXZJUO"
#   B = "Asx";  aspartic acid or asparagine (D or N)
#   X = "Xxx";  unknown or 'other' amino acid
#   Z = "Glx";  glutamic acid or glutamine (E or Q)
#   J = "Xle";  leucine or isoleucine (L or I, used in mass-spec)
#   U = "Sec";  selenocysteine
#   O = "Pyl";  pyrrolysine
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

ambiguous_rna_complement = {
    "A": "U",
    "C": "G",
    "G": "C",
    "U": "A",
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


def _make_ranges(mydict):
    d = {}
    for key, value in mydict.iteritems():
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

def _make_ambiguous_ranges(mydict, weight_table):
    range_d = {}
    avg_d = {}
    for letter, values in mydict.iteritems():
        #Following line is a quick hack to skip undefined weights for U and O
        if len(values)==1 and values[0] not in weight_table : continue
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
    "A": 89.09,
    "C": 121.16,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.18,
    "K": 146.19,
    "L": 131.18,
    "M": 149.21,
    "N": 132.12,
    #"O": 0.0, # Needs to be recorded!
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    #"U": 168.05, # To be confirmed
    "V": 117.15,
    "W": 204.23,
    "Y": 181.19
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
    "J": "IL",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "O": "O",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "U": "U",
    "V": "V",
    "W": "W",
    "X": "ACDEFGHIKLMNPQRSTVWY",
    #TODO - Include U and O in the possible values of X?
    #This could alter the extended_protein_weight_ranges ...
    "Y": "Y",
    "Z": "QE",
}
    
protein_weight_ranges = _make_ranges(protein_weights)

extended_protein_weight_ranges, avg_extended_protein_weights = \
               _make_ambiguous_ranges(extended_protein_values,
                                      protein_weights)



