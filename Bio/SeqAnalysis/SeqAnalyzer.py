import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis, IsoelectricPoint

records = {
    'czlowiek': 'MGKKRAPQPIVEKLISNFSQK',
    'pies': 'MSVKPSKKKSKRSKVKKKISFDFSDDDDSEIGVSFR',
}


class SeqAnalyzer:
    """
    A class that takes in records from SeqAnalysis.database, performs
    multiple operations on it and returns 2 dictionaries with result values. For further information about
    said dictionaries please see 'results()' function documentation.
    The operations that the class performs upon received data from SeqAnalysis.database is it:
    1. Calculates molecular weight of protein sequences
    2. Calculates isoelectric point for protein sequences
    3. Counts standard aa in protein sequences, in a list of dicts (amino_acids)
    4. Calculates aromaticity for protein sequences according to Lobry 1994 method
    5. Calculates amino acid percentage content based aa content in protein sequences
    6. Calculates instability indexes for protein sequences according to Guruprasad et al 1990
    7. Calculates flexibility of protein sequences according to Vihinen, 1994
    8. Calculates fraction of helix, turn and sheet for protein sequence
    9. Calculates GRAVY of protein sequences according to Kyte and Doolittle
    10. Calculates molar extinction coefficients for protein sequences
    11. Calculates median length of sequences
    12. Calculates average length of sequences
    13. Calculates standard deviation od sequence length
    14. Calculates 'N' and 'C' terminus amino acid occurrences frequencies
    """

    def __init__(self, data: dict):
        """
        This method is called when an object is created from the class and it allows the class to initialize and
        calculate via staticmethods usage(median_standard_deviation_average_length,  main_calculations
        terminal_aa_counter) attributes of a SeqAnalyzer

        :param data: takes in a dict with protein IDs as keys and SeqRecord objects as values
        """
        id_seq_dictionary = records  # {key: value.seq[0] for key, value in data}
        self.lengths = [len(seq) for seq in id_seq_dictionary.values()]
        self.ids = [ids for ids in id_seq_dictionary.keys()]
        self.standard_deviation, self.average_length, self.median = SeqAnalyzer.median_standard_deviation_average_length(
            id_seq_dictionary)
        self.molecular_weight, \
        self.pi, \
        self.amino_acids, \
        self.aromaticity, \
        self.amino_acids_percent, \
        self.instability, \
        self.flexibility, \
        self.secondary_structure_fraction, \
        self.gravy, \
        self.mol_ext_coefficient \
            = SeqAnalyzer.main_calculations(id_seq_dictionary)
        self.aa_n_terminus, self.aa_c_terminus = SeqAnalyzer.terminal_aa_counter(id_seq_dictionary)

    @staticmethod
    def median_standard_deviation_average_length(id_seq_dictionary: dict) -> tuple:
        """
        Calculates standard deviation of lengths of given sequences, median of given sequences lengths and
        average sequence length calculated from given protein sequences. Calculated values are later assigned to
        apropriate self. values :
        self.standard_deviation - standard deviation
        self.average_length - average length
        self.median - median length

        :parameter: takes in a dict with protein ID as keys and protein sequences as values

        :return: Calculated median, standard deviation of lengths and average length values for protein sequences in
        a form of a tuple

        :Example:
        >>> SeqAnalyzer.median_standard_deviation_average_length({'protein1': 'WKQTNSLEGKQ', 'protein2': 'WKQQTNSLEGKQ'})
        0.5, 11.5, 11.5

        """
        list_of_lengths = [len(seq) for seq in id_seq_dictionary.values()]
        standard_deviation = round(float(np.std(list_of_lengths)), 3)
        average_length = round(float(np.average(list_of_lengths)), 3)
        median_length = round(float(np.median(list_of_lengths)), 3)
        return standard_deviation, average_length, median_length

    @staticmethod
    def main_calculations(id_seq_dictionary: dict) -> tuple:
        """
        Main function calculating proteins values (10 values):
        1. Calculated molecular weight of protein sequences (molecular_weight )
        2. Calculated isoelectric point for protein sequences (pi)
        3. Counted standard aa in protein sequences, in a list of dicts (amino_acids)
        4. Calculated aromaticity for protein sequences according to Lobry 1994 method, in a list (aromaticity)
        5. Calculated amino acid percentage content based aa content in protein sequences, in a list (amino_acids_percent )
        6. Calculated instability indexes for protein sequences according to Guruprasad et al 1990, in a list (instability)
        7. Calculated flexibility of protein sequences according to Vihinen, 1994, in a list (flexibility)
        8. Calculated fraction of helix, turn and sheet for protein sequences, in a list
        of tuples with 3 elements (helix, turn, sheet) fractions accordingly (secondary_structure_fraction)
        9. Calculated gravy of protein sequences according to Kyte and Doolittle, in a list (seq_gravy )
        10. Calculated molar extinction coefficients for protein sequences, in a list (mol_ext_coefficient)

        :return: Output parameters are in list of:
        1. molecular_weight - of ints and/or floats
        2. pi - of ints and/or floats
        3. amino_acids - dictionaries
        4. aromaticity - of ints and/or floats
        5. amino_acids_percent - of ints and/or floats
        6. instability - of ints and/or floats
        7. flexibility - of ints and/or floats
        8. secondary_structure_fraction - of tuples (of ints and/or floats) with 3 elements
        9. seq_gravy - of ints and/or floats
        10. mol_ext_coefficient - of two elemental tuples of ints and/or floats

        :Example:   #need to finish
        >>> SeqAnalyzer.main_calculations({'protein1': 'WKQTNSLEGKQ', 'protein2': 'WKQQTNSLEGKQ'})
        [1318.4349, 1446.5641000000003],
        [8.59112548828125, 8.59112548828125],
        """
        list_of_proteinanalysis_objects = [ProteinAnalysis(seq) for seq in id_seq_dictionary.values()]
        molecular_weight = [round(seq.molecular_weight(), 3) for seq in list_of_proteinanalysis_objects]
        pi = [round(seq.isoelectric_point(), 3) for seq in list_of_proteinanalysis_objects]
        amino_acids = [seq.count_amino_acids() for seq in list_of_proteinanalysis_objects]
        aromaticity = [round(seq.aromaticity(), 3) for seq in list_of_proteinanalysis_objects]
        amino_acids_percent = [dict(seq.get_amino_acids_percent()) for seq in list_of_proteinanalysis_objects]
        for dict_aa_percent in amino_acids_percent:  # rounding values in aa percentage occurrences dictionaries
            for key, value in dict_aa_percent.items():
                dict_aa_percent[key] = round(value, 3)

        instability = [round(seq.instability_index(), 3) for seq in list_of_proteinanalysis_objects]
        flexibility = [round(sum(seq) / len(seq), 3) for seq in
                       [aa.flexibility() for aa in list_of_proteinanalysis_objects]]
        secondary_structure_fraction = [[round(value, 3) for value in [v1, v2, v3]] for v1, v2, v3 in
                                        (seq.secondary_structure_fraction() for seq in
                                         list_of_proteinanalysis_objects)]
        gravy = [round(seq.gravy(), 3) for seq in list_of_proteinanalysis_objects]
        mol_ext_coefficient = [[round(value, 3) for value in [v1, v2]] for v1, v2 in
                               (seq.molar_extinction_coefficient() for seq in list_of_proteinanalysis_objects)]
        return molecular_weight, pi, amino_acids, aromaticity, amino_acids_percent, instability, flexibility, \
               secondary_structure_fraction, gravy, mol_ext_coefficient

    @staticmethod
    def terminal_aa_counter(id_seq_dictionary: dict) -> tuple:
        """
        Counts occurrences of given aa on 'C' and 'N' terminus of protein and makes a dictionary for each terminus
        ('C' and 'N') with those aa occurrences

        :parameter: takes in a dict with protein ID as keys and protein sequences as values

        :return: returns a tuple with two dicts, where aa_n_terminius is a dictionary with information about occurrences
        of aa on 'N' terminus, and aa_c_terminus about 'C' terminus respectively

        :Example:
        >>> SeqAnalyzer.terminal_aa_counter({'protein1': 'WKQTNSLEGKQ', 'protein2': 'WKQQTNSLEGKQ'})
        'n_term_freq': {'Q': 2.0}, 'c_term_freq': {'W': 2.0}

        """
        list_c_terminus = [seq[-1] for seq in id_seq_dictionary.values()]
        aa_n_terminus = {i: round(list_c_terminus.count(i) / len(set(list_c_terminus)), 3) for i in set(list_c_terminus)}

        list_n_terminus = [seq[0] for seq in id_seq_dictionary.values()]
        aa_c_terminus = {i: round(list_n_terminus.count(i) / len(set(list_n_terminus)), 3) for i in set(list_n_terminus)}
        return aa_n_terminus, aa_c_terminus

    def results(self):
        """
        Function created with the aim to export calculated data about protein sequences for further visualization. The
        function returns said data in a form of two dictionaries. One dict (dict_all_sequences) contains data about
        protein sequences as a whole and second dict (dict_each_sequence) contains a data for each sequence separately
        where in (dict_each_sequence) ids of sequences are keys and values are dicts of data for given sequence

        :return: returns two dictionaries: dict_all_sequences and dict_each_sequence
        """
        dict_all_sequences = \
            {
                'median_id': self.median,
                'average_length_id': self.average_length,
                'standard_deviation_id': self.standard_deviation,
                'id': self.ids,
                'sequence_length': self.lengths,
                'molecular_weight': self.molecular_weight,
                'theoretical_pl': self.pi,
                'instability_index': self.instability,
                'flexibility': self.flexibility,
                'secondary_structure_fraction': self.secondary_structure_fraction,
                'extinction_coefficient': self.mol_ext_coefficient,
                'gravy': self.gravy,
                'c_term_freq': self.aa_c_terminus,
                'n_term_freq': self.aa_n_terminus,
                'amino_acid_comp': self.amino_acids,
                'amino_acid_percent': self.amino_acids_percent,
                'aromaticity': self.aromaticity
            }

        dict_each_sequence = {ids: {} for ids in self.ids}  # creating empty dict of dicts

        for index, key in enumerate(dict_each_sequence):
            dict_each_sequence[key]['sequence_length'] = self.lengths[index]
            dict_each_sequence[key]['molecular_weight'] = self.molecular_weight[index]
            dict_each_sequence[key]['theoretical_pl'] = self.pi[index]
            dict_each_sequence[key]['instability_index'] = self.instability[index]
            dict_each_sequence[key]['flexibility'] = self.instability[index]
            dict_each_sequence[key]['secondary_structure_fraction'] = self.secondary_structure_fraction[index]
            dict_each_sequence[key]['extinction_coefficient'] = self.mol_ext_coefficient[index]
            dict_each_sequence[key]['amino_acid_comp'] = self.amino_acids[index]
            dict_each_sequence[key]['amino_acid_percent'] = self.amino_acids_percent[index]
            dict_each_sequence[key]['aromaticity'] = self.aromaticity[index]
            dict_each_sequence[key]['gravy'] = self.gravy[index]

        return dict_all_sequences, dict_each_sequence


wynik = SeqAnalyzer(records)
print(wynik.results())
