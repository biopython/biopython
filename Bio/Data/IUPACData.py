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


# For Center of Mass Calculation.
# Taken from http://www.chem.qmul.ac.uk/iupac/AtWt/ & PyMol
atom_weights = {
    'H'  :   1.00794,
    'He' :   4.002602,
    'Li' :   6.941,
    'Be' :   9.012182,
    'B'  :  10.811,
    'C'  :  12.0107,
    'N'  :  14.0067,
    'O'  :  15.9994,
    'F'  :  18.9984032,
    'Ne' :  20.1797,
    'Na' :  22.989770,
    'Mg' :  24.3050,
    'Al' :  26.981538,
    'Si' :  28.0855,
    'P'  :  30.973761,
    'S'  :  32.065,
    'Cl' :  35.453,
    'Ar' :  39.948,
    'K'  :  39.0983,
    'Ca' :  40.078,
    'Sc' :  44.955910,
    'Ti' :  47.867,
    'V'  :  50.9415,
    'Cr' :  51.9961,
    'Mn' :  54.938049,
    'Fe' :  55.845,
    'Co' :  58.933200,
    'Ni' :  58.6934,
    'Cu' :  63.546,
    'Zn' :  65.39,
    'Ga' :  69.723,
    'Ge' :  72.64,
    'As' :  74.92160,
    'Se' :  78.96,
    'Br' :  79.904,   
    'Kr' :  83.80,
    'Rb' :  85.4678,
    'Sr' :  87.62,
    'Y'  :  88.90585,
    'Zr' :  91.224,
    'Nb' :  92.90638,
    'Mo' :  95.94,
    'Tc' :  98.0,
    'Ru' : 101.07,
    'Rh' : 102.90550,
    'Pd' : 106.42,
    'Ag' : 107.8682,
    'Cd' : 112.411,
    'In' : 114.818,
    'Sn' : 118.710,
    'Sb' : 121.760,
    'Te' : 127.60,
    'I'  : 126.90447,
    'Xe' : 131.293,
    'Cs' : 132.90545,
    'Ba' : 137.327,
    'La' : 138.9055,
    'Ce' : 140.116,
    'Pr' : 140.90765,
    'Nd' : 144.24,
    'Pm' : 145.0,
    'Sm' : 150.36,
    'Eu' : 151.964,
    'Gd' : 157.25,
    'Tb' : 158.92534,
    'Dy' : 162.50,
    'Ho' : 164.93032,
    'Er' : 167.259,
    'Tm' : 168.93421,
    'Yb' : 173.04,
    'Lu' : 174.967,
    'Hf' : 178.49,
    'Ta' : 180.9479,
    'W'  : 183.84,
    'Re' : 186.207,
    'Os' : 190.23,
    'Ir' : 192.217,
    'Pt' : 195.078,
    'Au' : 196.96655,
    'Hg' : 200.59,
    'Tl' : 204.3833,
    'Pb' : 207.2,
    'Bi' : 208.98038,
    'Po' : 208.98,
    'At' : 209.99,
    'Rn' : 222.02,
    'Fr' : 223.02,
    'Ra' : 226.03,
    'Ac' : 227.03,
    'Th' : 232.0381,
    'Pa' : 231.03588,
    'U'  : 238.02891,
    'Np' : 237.05,
    'Pu' : 244.06,
    'Am' : 243.06,
    'Cm' : 247.07,
    'Bk' : 247.07,
    'Cf' : 251.08,
    'Es' : 252.08,
    'Fm' : 257.10,
    'Md' : 258.10,
    'No' : 259.10,
    'Lr' : 262.11,
    'Rf' : 261.11,
    'Db' : 262.11,
    'Sg' : 266.12,
    'Bh' : 264.12,
    'Hs' : 269.13,
    'Mt' : 268.14,
}

# pKa of amino acids side chains
# taken from EMBOSS (/usr/share/EMBOSS/data/Epk.dat)
protein_pka_side_chain = {
    "A": None,
    "C": 8.5,
    "D": 3.9,
    "E": 4.1,
    "F": None,
    "G": None,
    "H": 6.5,
    "I": None,
    "K": 10.8,
    "L": None,
    "M": None,
    "N": None,
    "P": None,
    "Q": None,
    "R": 12.5,
    "S": None,
    "T": None,
    "V": None,
    "W": None,
    "Y": 10.1
    }

# Relative hydrophobicity of amino acids taken
# from /usr/share/EMBOSS/data/Eaa_hydropathy.dat

# Kyte-Doolittle
protein_hydropathy_KD = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3
    }

# Sweet & Eisenberg
protein_hydropathy_OHM = {
    "A": -0.40,
    "C": 0.17,
    "D": -1.31,
    "E": -3.5,
    "F": -1.22,
    "G": 1.92,
    "H": -0.64,
    "I": 1.25,
    "K": -0.67,
    "L": 1.22,
    "M": 1.02,
    "N": -0.92,
    "P": -0.49,
    "Q": -0.91,
    "R": -0.59,
    "S": -0.55,
    "T": -0.28,
    "V": 0.91,
    "W": 0.50,
    "Y": 1.67
    }

# Consensus
protein_hydropathy_consensus = {
    "A": 0.62,
    "C": 0.29,
    "D": -0.90,
    "E": -0.74,
    "F": 1.19,
    "G": 0.48,
    "H": -0.40,
    "I": 1.38,
    "K": -1.50,
    "L": 1.06,
    "M": 0.64,
    "N": -0.78,
    "P": 0.12,
    "Q": -0.85,
    "R": -2.53,
    "S": -0.18,
    "T": -0.05,
    "V": 1.08,
    "W": 0.81,
    "Y": 0.26
    }
