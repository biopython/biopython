import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis, IsoelectricPoint

records = {
    'czlowiek': 'MGKKRAPQKGKTVTKPQEIIVDESKLNWKPVDIPDTLDDFGGFYGLEEIDGVDVKVVDGKVTFVTKKDSKVLKDSNKEKVGDDQESVENESGSDSESELLEFKNLDDIKEGELSAASYSSSDEDEQGNIESSKLTDPSEDVDEDVDEDIPIVEKLISNFSQK',
    'pies': 'METIDSKQNINRESLLEERRKKLAKWKQKKAQFDAQKEHQTSRNDIVTNSLEGKQTTEKFTERQERVKEELRKRKNEFRKSDEPVSVKPSKKKSKRSKVKKKISFDFSDDDDSEIGVSFR',
    'kot': 'WKQKKAQFDAQKEHQTSRNDIVTNSLEGKQTTEKFTERQECCCCCCCRVMGKKRAPQKGKTVTKPQEIIVDESKLNWKPVDIPDTLDDFGGFYGLEEIDGVDVKVVDGKVTFVTKKDSKVLKDSNKEKVGDDQESVENESGSDSESELLEFKNLDDIKEGELSAASYSSSDEDEQGNIESSKLTDPSEDVDEDVDEDIPIVEKLISNFSQK'
}


class Analyzer:
    """
    A class that takes in records from SeqAnalysis.database and performs
    multiple operations on it and returns pd.DataFrame with values
    """

    def makedict(self, records) -> dict:
        """
        Makes a dictionary from records, where ids are keys and sequences are values

        :return: returns a created dictionary
        """
        dict = {key: value.seq[0] for key, value in records}
        return dict

    def median_stdev_avg(self) -> tuple:
        """
        Calculates standard deviation of lenths of given sequences

        :return: changes self.stdev values to standard value calculate from sequences
        """
        list_of_seq = self.dict.values()
        list_of_lenths = [len(_) for _ in list_of_seq]
        stdev = float(np.std(list_of_lenths))
        avg = float(np.average(list_of_lenths))
        median = float(np.median(list_of_lenths))
        return stdev, avg, median

    def seq_lenth(self) -> list:
        """
        Makes a list of sequences lenths

        :return: returns a list of sequences lenths
        """
        return [len(_[1]) for _ in self.tuples]

    def id(self) -> list:
        """
        Makes a list of protein IDs

        :return: ruturns a list od protein ids
        """
        return [_[0] for _ in self.tuples]

    def molecular_weight(self) -> tuple:
        """
        Main funcion calculating proteins values (10 values).

        :return: Funcion returns a tuple containning:
        1. molecular_weight - calculated malecular weight from protein sequences in a list
        2. pi - calculated isoelectric point for protein sequences in a list
        3. amino_acids - counted standard aa in protein sequences, in a list of dicts
        4. aromacity - calculated aromacity for protein sequences according to Lobry 1994 method, in a list
        5. amino_acids_percent - calculated percantage based aa content in protein sequences, in a list
        6. instability - calculated instability indexes for protein sequences according to Guruprasad et al 1990, in a list
        7. flexibility - calculated flexibility of protein sequences according to Vihinen, 1994, in a list
        8. secondary_structure_fraction - calculated fraction of helix, turn and sheet for protein sequences, in a list
        of tuples with 3 elements (helix, turn, sheet) fractions accordingly
        9. seq_gravy - calculated gravy of protein sequences according to Kyte and Doolittle, in a list
        10.mol_ext_coefficient - calculated molar extinction coefficients for protein sequences, in a list
        """
        l_objects = [ProteinAnalysis(_[1]) for _ in self.tuples]
        molecular_weight = [_.molecular_weight() for _ in l_objects]
        pi = [_.isoelectric_point() for _ in l_objects]
        amino_acids = [_.count_amino_acids() for _ in l_objects]
        aromacity = [_.aromaticity() for _ in l_objects]
        amino_acids_percent = [_.get_amino_acids_percent() for _ in l_objects]
        instability = [_.instability_index() for _ in l_objects]
        flexibility = [sum(__) / len(__) for __ in [_.flexibility() for _ in l_objects]]
        secondary_structure_fraction = [_.secondary_structure_fraction() for _ in l_objects]
        seq_gravy = [_.gravy() for _ in l_objects]
        mol_ext_coefficient = [_.molar_extinction_coefficient()[0] for _ in l_objects]
        return molecular_weight, pi, amino_acids, aromacity, amino_acids_percent, instability, flexibility, \
               secondary_structure_fraction, seq_gravy, mol_ext_coefficient

    def n_and_c_term(self) -> tuple:
        """
        Counts occurrences of given aa on 'C' and 'N' terminus of protein

        :return: returns a list of two dicts, where list[0] is a dictionary with information about occurrences
        of aa on 'N' terminus, and list[1] about 'C' terminus respectively
        """
        list_c_term = [_[-1] for _ in self.dict.values()]
        dict_c_term = {i: list_c_term.count(i) / len(set(list_c_term)) for i in set(list_c_term)}

        list_n_term = [_[0] for _ in self.dict.values()]
        dict_n_term = {i: list_n_term.count(i) / len(set(list_n_term)) for i in set(list_n_term)}
        return dict_n_term, dict_c_term

    def __init__(self, records):
        self.dict = records  # Analyzer.makedict(self, records)  # creates a dict with seq ids and seqeances itself
        self.tuples = [(id, seq) for id, seq in self.dict.items()]
        self.stdev, self.average, self.median = Analyzer.median_stdev_avg(self)
        self.ids = Analyzer.id(self)
        self.lenths = Analyzer.seq_lenth(self)
        self.objects = Analyzer.molecular_weight(self)

        self.molecular_weight, \
        self.pi, \
        self.amino_acids, \
        self.aromaticity, \
        self.amino_acids_percent, \
        self.instability, \
        self.flexibility, \
        self.secondary_structure_fraction, \
        self.gravy, \
        self.mol_ext_coeff \
        = Analyzer.molecular_weight(self)

        self.n_term, self.c_term = Analyzer.n_and_c_term(self)

    def results(self):
        dict_1 = {
            'median_id': self.median,
            'mean_id': self.average,
            'sd_id': self.stdev,
            'id': self.ids,
            'sequence_length': self.lenths,
            'molecular_weight': self.molecular_weight,
            'theoretical_pl': self.pi,
            'instability_index': self.instability,
            'flexibility': self.flexibility,
            'secondary_structure_fraction': self.secondary_structure_fraction,
            'extinction_coefficient': self.mol_ext_coeff,
            'gravy': self.gravy,
            'c_term_freq': self.c_term,
            'n_term_freq': self.n_term,
            'amino_acid_comp': self.amino_acids,
            'amino_acid_percent': self.amino_acids_percent,
            'aromacity': self.aromaticity
        }

        return dict_1


Analyzer = Analyzer(records)

print(Analyzer.results())



""" dict_2
{'system': {'median_id': '0',
              'mean_id': '0',
               'sd_id': '0',
                'id': ['system'],
                'sequence_length': [340.75],
                'molecular_weight': [100.77],
                 'theoretical_pl': [0.2],
                'extinction_coefficient': [61.04],
                 'instability_index': [0.93],
                  'gravy': [56.1], 
                  'est_half_life': [692.71], 
                  'c_term_freq': {'A': 0.105, 'G': 0.895}, 
                  'n_term_freq': {'A': 0.94, 'G': 0.06000000000000005}, 
                  'amino_acid_comp': [{'T': 420, 'L': 473, 'A': 209, 'S': 147, 'G': 94, 'M': 402, 'I': 302, 'F': 261, 'V': 257, 'E': 303, 'H': 403, 'N': 372, 'Q': 78, 'K': 387, 'P': 99, 'W': 123, 'R': 352}], 
                  'atomic_comp': [{'H': 123, 'N': 387, 'O': 52}]}, """
