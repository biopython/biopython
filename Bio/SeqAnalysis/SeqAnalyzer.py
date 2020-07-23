# Copyright 2019-2020 by moozeq, KacDom, erpeg and s0lczi. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


import json

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class SeqAnalyzer:
    """
    A class that takes in records from SeqAnalysis.database, performs multiple operations on them and returns
    2 dictionaries with result values(rounded to three decimal places). It is important to point out that in cases where
    aa sequences have positions containing characters: 'B', 'Z', 'X', 'J', which characters mean that aa in this
    position aren't certain - the class rejects such sequences, doesn't perform any calculation on them and doesn't
    pass them on to SeqVisualizer. The same goes to protein sequences shorter than 10 aa.
    For further information about said dictionaries please see 'results()' function documentation.
    The operations that the class performs upon received data from SeqAnalysis.database is it:

    1. Calculates molecular weight of protein sequences
    2. Calculates isoelectric point for protein sequences
    3. Counts standard aa in protein sequences
    4. Calculates aromaticity for protein sequences according to Lobry 1994 method
    5. Calculates amino acid percentage content based aa content in protein sequences
    6. Calculates instability indexes for protein sequences according to Guruprasad et al 1990
    7. Calculates flexibility of protein sequences according to Vihinen, 1994
    8. Calculates percentage of helix, turn and sheet in protein sequences
    9. Calculates GRAVY of protein sequences according to Kyte and Doolittle
    10. Calculates molar extinction coefficients for protein sequences
    11. Calculates median length of sequences
    12. Calculates average length of sequences
    13. Calculates standard deviation od sequence length
    14. Calculates 'N' and 'C' terminus amino acid occurrences frequencies
    """

    def __init__(self, data: dict, analysis_root_dir: str):
        """
        This method is called when an object is created from the class and it allows the class to initialize and
        calculate via static methods usage(median_standard_deviation_average_length, main_calculations
        terminal_aa_counter) attributes of a SeqAnalyzer. Calculated values are rounded to round_values decimal places

        :param data: takes in a dict with protein IDs as keys and SeqRecord objects as values
        :param analysis_root_dir: root directory for analysis
        """
        round_value = 3
        id_seq_dictionary = {
            id: aa_seq
            for id, aa_seq
            in {id: str(aa_seq.seq) for id, aa_seq in data.items()}.items()
            if not (set('BZXJ').intersection(aa_seq)) and len(aa_seq) >= 10
        }
        self.round = self.terminus_aa = [(seq[0], seq[-1]) for seq in id_seq_dictionary.values()]
        self.lengths = [len(seq) for seq in id_seq_dictionary.values()]
        self.ids = [ids for ids in id_seq_dictionary.keys()]
        self.standard_deviation, self.average_length, self.median = \
            SeqAnalyzer.median_standard_deviation_average_length(id_seq_dictionary, round_value)
        self.molecular_weight, self.pi, self.amino_acids, self.aromaticity, self.amino_acids_percent, \
            self.instability, self.flexibility, self.secondary_structure_fraction, self.gravy, \
            self.mol_ext_coefficient = SeqAnalyzer.main_calculations(id_seq_dictionary, round_value)
        self.aa_n_terminus, self.aa_c_terminus = SeqAnalyzer.terminal_aa_counter(id_seq_dictionary, round_value)
        with open(f'{analysis_root_dir}/results.json', 'w') as res_file:
            json.dump(self.results(), res_file, indent=4)

    @staticmethod
    def median_standard_deviation_average_length(id_seq_dictionary: dict, round_value: int) -> tuple:
        """
        Calculates standard deviation of lengths of given sequences, median of given sequences lengths and
        average sequence length calculated from given protein sequences. Calculated values (rounded to round_value
        decimal places) are later assigned to apropriate self. values :
        self.standard_deviation - standard deviation
        self.average_length - average length
        self.median - median length

        :parameter: takes in a dict with protein ID as keys and protein sequences as values

        :return: Calculated median, standard deviation of lengths and average length values for protein sequences in
        a form of a tuple

        :Example:
        >>> SeqAnalyzer.median_standard_deviation_average_length({'protein1': 'WKQTNSLEGKQ', 'protein2': 'WKQQTNSLEGKQ'}, 3)
        (0.5, 11.5, 11.5)
        """
        list_of_lengths = [len(seq) for seq in id_seq_dictionary.values()]
        standard_deviation = round(float(np.std(list_of_lengths)), round_value)
        average_length = round(float(np.average(list_of_lengths)), round_value)
        median_length = round(float(np.median(list_of_lengths)), round_value)

        return standard_deviation, average_length, median_length

    @staticmethod
    def main_calculations(id_seq_dictionary: dict, round_value: int) -> tuple:
        """
        Main function calculating proteins values (10 values):
        1. Calculated molecular weight of protein sequences (molecular_weight )
        2. Calculated isoelectric point for protein sequences (pi)
        3. Counted standard aa in protein sequences, in a list of dicts (amino_acids)
        4. Calculated aromaticity for protein sequences according to Lobry 1994 method, in a list (aromaticity)
        5. Calculated amino acid percentage content based aa content in protein sequences, in a list (amino_acids_percent)
        6. Calculated instability indexes for protein sequences according to Guruprasad et al 1990, in a list(instability)
        7. Calculated flexibility of protein sequences according to Vihinen, 1994, in a list (flexibility)
        8. Calculated fraction of helix, turn and sheet for protein sequences, in a dict with (helix, turn, sheet) as
        keys and ther percentage in protein structure as values
        9. Calculated gravy of protein sequences according to Kyte and Doolittle, in a list (seq_gravy )
        10. Calculated molar extinction coefficients for protein sequences, in a dict (mol_ext_coefficient)
        Calculated values are rounded to round_value decimal places

        :return: Output parameters are in list of:
        1. molecular_weight - of ints and/or floats
        2. pi - of ints and/or floats
        3. amino_acids - dictionaries
        4. aromaticity - of ints and/or floats
        5. amino_acids_percent - of ints and/or floats
        6. instability - of ints and/or floats
        7. flexibility - of ints and/or floats
        8. secondary_structure_fraction - of dicts with three keys/values with ints and/or floats as values
        9. seq_gravy - of ints and/or floats
        10. mol_ext_coefficient - of dicts with two keys/values with ints and/or floats as values

        :Example:   #need to finish
        >>> SeqAnalyzer.main_calculations({'protein1': 'WKQTNSTDLHEYTLEGKQ', 'protein2': 'TDLHEYTWKQQTNSLEGKQ'}, 3)
        ([2178.314, 2306.443], [5.445, 5.414], [{'A': 0, 'C': 0, 'D': 1, 'E': 2, 'F': 0, 'G': 1, 'H': 1, 'I': 0, 'K': 2, 'L': 2, 'M': 0, 'N': 1, 'P': 0, 'Q': 2, 'R': 0, 'S': 1, 'T': 3, 'V': 0, 'W': 1, 'Y': 1}, {'A': 0, 'C': 0, 'D': 1, 'E': 2, 'F': 0, 'G': 1, 'H': 1, 'I': 0, 'K': 2, 'L': 2, 'M': 0, 'N': 1, 'P': 0, 'Q': 3, 'R': 0, 'S': 1, 'T': 3, 'V': 0, 'W': 1, 'Y': 1}], [0.111, 0.105], [{'A': 0.0, 'C': 0.0, 'D': 5.556, 'E': 11.111, 'F': 0.0, 'G': 5.556, 'H': 5.556, 'I': 0.0, 'K': 11.111, 'L': 11.111, 'M': 0.0, 'N': 5.556, 'P': 0.0, 'Q': 11.111, 'R': 0.0, 'S': 5.556, 'T': 16.667, 'V': 0.0, 'W': 5.556, 'Y': 5.556}, {'A': 0.0, 'C': 0.0, 'D': 5.263, 'E': 10.526, 'F': 0.0, 'G': 5.263, 'H': 5.263, 'I': 0.0, 'K': 10.526, 'L': 10.526, 'M': 0.0, 'N': 5.263, 'P': 0.0, 'Q': 15.789, 'R': 0.0, 'S': 5.263, 'T': 15.789, 'V': 0.0, 'W': 5.263, 'Y': 5.263}], [17.928, 19.737], [1.005, 1.009], [{'Helix': 22.222, 'Turn': 16.667, 'Sheet': 22.222}, {'Helix': 21.053, 'Turn': 15.789, 'Sheet': 21.053}], [-1.661, -1.758], [{'Cysteines reduced': 6990, 'Cysteines residues': 6990}, {'Cysteines reduced': 6990, 'Cysteines residues': 6990}])

        """
        list_of_proteinanalysis_objects = [ProteinAnalysis(seq) for seq in id_seq_dictionary.values()]
        molecular_weight = [round(seq.molecular_weight(), round_value) for seq in list_of_proteinanalysis_objects]
        pi = [round(seq.isoelectric_point(), round_value) for seq in list_of_proteinanalysis_objects]
        amino_acids = [seq.count_amino_acids() for seq in list_of_proteinanalysis_objects]
        aromaticity = [round(seq.aromaticity(), round_value) for seq in list_of_proteinanalysis_objects]
        amino_acids_percent = [dict(seq.get_amino_acids_percent()) for seq in list_of_proteinanalysis_objects]
        for dict_aa_percent in amino_acids_percent:  # rounding values in aa percentage occurrences dictionaries
            for key, value in dict_aa_percent.items():
                dict_aa_percent[key] = round(value * 100, round_value)

        instability = [round(seq.instability_index(), round_value) for seq in list_of_proteinanalysis_objects]
        flexibility = [round(sum(seq) / len(seq), round_value) for seq
                       in [aa.flexibility() for aa in list_of_proteinanalysis_objects]]
        secondary_structure_fraction = [
            {'Helix': round(v1 * 100, round_value), 'Turn': round(v2 * 100, round_value),
             'Sheet': round(v3 * 100, round_value)}
            for v1, v2, v3
            in (seq.secondary_structure_fraction() for seq in list_of_proteinanalysis_objects)
        ]
        gravy = [round(seq.gravy(), round_value) for seq in list_of_proteinanalysis_objects]
        mol_ext_coefficient = [
            {'Cysteines reduced': v1, 'Cysteines residues': v2} for v1, v2 in
            (seq.molar_extinction_coefficient() for seq in list_of_proteinanalysis_objects)
        ]
        return molecular_weight, pi, amino_acids, aromaticity, amino_acids_percent, instability, flexibility, secondary_structure_fraction, gravy, mol_ext_coefficient

    @staticmethod
    def terminal_aa_counter(id_seq_dictionary: dict, round_value: int) -> tuple:
        """
        Counts occurrences of given aa on 'C' and 'N' terminus of protein and makes a dictionary for each terminus
        ('C' and 'N') with those aa occurrences. Calculated values are rounded to round_value decimal places

        :parameter: takes in a dict with protein ID as keys and protein sequences as values

        :return: returns a tuple with two dicts, where aa_n_terminius is a dictionary with information about occurrences
        of aa on 'N' terminus, and aa_c_terminus about 'C' terminus respectively

        :Example:
        >>> SeqAnalyzer.terminal_aa_counter({'protein1': 'WKQTNSTDLHEYTLEGKQ', 'protein2': 'TDLHEYTWKQQTNSLEGKQ'}, 3)\
        == ({'Q': 2.0}, {'T': 0.5, 'W': 0.5})
        True
        """
        list_c_terminus = [seq[-1] for seq in id_seq_dictionary.values()]
        aa_n_terminus = {
            i: round(list_c_terminus.count(i) / len(set(list_c_terminus)), round_value) for i in set(list_c_terminus)
        }

        list_n_terminus = [seq[0] for seq in id_seq_dictionary.values()]
        aa_c_terminus = {
            i: round(list_n_terminus.count(i) / len(set(list_n_terminus)), round_value) for i in set(list_n_terminus)
        }
        return aa_n_terminus, aa_c_terminus

    def results(self) -> tuple:
        """
        Function created with the aim to export calculated data about protein sequences for further visualization. The
        function returns said data in a form of two dictionaries. One dict (dict_all_sequences) contains data about
        protein sequences as a whole and second dict (dict_each_sequence) contains a data for each sequence separately
        where in (dict_each_sequence) ids of sequences are keys and values are dicts of data for given sequence

        :return: returns two dictionaries: dict_all_sequences and dict_each_sequence
        """
        dict_all_sequences = {
            'median_length': self.median,
            'average_length': self.average_length,
            'standard_deviation_length': self.standard_deviation,
            'id': self.ids,
            'sequence_length': self.lengths,
            'molecular_weight': self.molecular_weight,
            'theoretical_pi': self.pi,
            'instability_index': self.instability,
            'flexibility': self.flexibility,
            'gravy': self.gravy,
            'aromaticity': self.aromaticity,
            'c_term_freq': self.aa_c_terminus,
            'n_term_freq': self.aa_n_terminus,
            'secondary_structure_fraction': self.secondary_structure_fraction,
            'extinction_coefficient': self.mol_ext_coefficient,
            'amino_acid_comp': self.amino_acids,
            'amino_acid_percent': self.amino_acids_percent,
        }

        dict_each_sequence = {ids: {} for ids in self.ids}  # creating empty dict of dicts

        for index, key in enumerate(dict_each_sequence):
            dict_each_sequence[key]['id'] = [self.ids[index]]
            dict_each_sequence[key]['sequence_length'] = [self.lengths[index]]
            dict_each_sequence[key]['molecular_weight'] = [self.molecular_weight[index]]
            dict_each_sequence[key]['theoretical_pi'] = [self.pi[index]]
            dict_each_sequence[key]['instability_index'] = [self.instability[index]]
            dict_each_sequence[key]['flexibility'] = [self.flexibility[index]]
            dict_each_sequence[key]['gravy'] = [self.gravy[index]]
            dict_each_sequence[key]['aromaticity'] = [self.aromaticity[index]]
            dict_each_sequence[key]['secondary_structure_fraction'] = [self.secondary_structure_fraction[index]]
            dict_each_sequence[key]['extinction_coefficient'] = [self.mol_ext_coefficient[index]]
            dict_each_sequence[key]['amino_acid_comp'] = [self.amino_acids[index]]
            dict_each_sequence[key]['amino_acid_percent'] = [self.amino_acids_percent[index]]
            dict_each_sequence[key]['c_term_freq'] = {self.terminus_aa[index][1]: 1}
            dict_each_sequence[key]['n_term_freq'] = {self.terminus_aa[index][0]: 1}

        return dict_all_sequences, dict_each_sequence
