import pandas as pd
import numpy as np

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
        slownik = {_[0]: _[1][0] for _ in records}
        return slownik

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

    def __init__(self, records):
        self.dict = analyzer.makedict(self, records)  # creates a dict with seq ids and seqeances itself
        self.tuples = [(id, seq) for id, seq in self.dict.values()]
        self.stdev, self.average, self.median = analyzer.median_stdev_avg(self)

    def results(self):
        data = pd.DataFrame(data={
            'id': [id[0] for id in self.tuples].append(['median', 'standard_deviation', 'average_lenth']),
            'sequence_length': [len(seq[1]) for seq in self.tuples].append([self.median, self.stdev, self.average]),
            'molecular_weight': [],
            'theoretical_pl': [],
            'extinction_coefficient': [],
            'instability_index': [],
            'gravy': [],
            'est_half_life': [],
            'c_term_freq': [{'A': 0.2, 'G': 0.8}, {'A': 0.3, 'L': 0.7}],
            'n_term_freq': [{'A': 0.8, 'G': 0.2}, {'A': 0.7, 'L': 0.3}],
            'amino_acid_comp': [{'A': 12, 'G': 24}, {'L': 124, 'A': 52}],
            'atomic_comp': [{'C': 42, 'H': 335, 'N': 55}, {'C': 55, 'H': 523, 'N': 13, 'O': 53}]
        })

        df = pd.DataFrame(data)
        return df


""" dict_1
{'median_id': '0',
 'mean_id': '0',
 'sd_id': '0',
 'id': ['remain', 'correct', 'margin', 'spring', 'toilet', 'certain'], 
 'sequence_length': [123.91, 142.41, 322.0, 896.65, 1145.16, 379.1],
 'molecular_weight': [121.91, 489.38, 356.37, 328.72, 184.74, 338.53],
 'theoretical_pl': [1.83, 0.79, 3.34, 3.96, 0.66, 3.08],
 'extinction_coefficient': [71.9, 97.17, 78.98, 36.03, 89.86, 34.71], 
 'instability_index': [1.99, 3.99, 2.68, 3.94, 4.54, 2.29],
 'gravy': [76.67, 3.81, 85.85, 50.89, 16.16, 86.14], 
 'est_half_life': [342.98, 1057.0, 296.45, 958.88, 1196.21, 112.16], 
 'c_term_freq': {'A': 0.136, 'G': 0.864},
 'n_term_freq': {'A': 0.903, 'G': 0.09699999999999998}, 
 'amino_acid_comp': [{'C': 271, 'K': 152, 'Y': 313, 'Q': 129, 'H': 400, 'W': 443, 'L': 74, 'F': 467, 'A': 146, 'G': 102, 'M': 233, 'N': 357, 'S': 291, 'R': 376, 'D': 322, 'P': 124, 'V': 83, 'I': 266, 'T': 263}, 
                     {'V': 265, 'T': 350, 'D': 120, 'R': 160, 'E': 63, 'Y': 129}, {'F': 276, 'G': 371, 'E': 329, 'N': 208, 'W': 64, 'Y': 263, 'P': 361, 'I': 327, 'A': 180, 'T': 359, 'D': 203, 'C': 318},
                     {'W': 47, 'T': 438, 'M': 211, 'I': 364, 'N': 476, 'C': 273, 'F': 438, 'L': 450, 'S': 404, 'A': 408, 'Y': 195, 'V': 292, 'K': 341, 'R': 477, 'G': 278, 'Q': 64, 'H': 242, 'D': 491}, 
                     {'H': 216, 'A': 237, 'L': 367, 'T': 173, 'W': 14, 'N': 206, 'S': 121, 'V': 306, 'P': 18, 'I': 478, 'Q': 137, 'C': 374, 'E': 131},
                     {'D': 492, 'I': 90, 'V': 372, 'L': 439, 'T': 129, 'C': 116, 'W': 195, 'K': 150, 'E': 9}], 
 'atomic_comp': [{'C': 430, 'O': 312, 'N': 51}, {'N': 192, 'C': 462, 'O': 438}, {'H': 220, 'C': 388, 'N': 244}, {'H': 113, 'N': 358, 'O': 279}, {'O': 80, 'C': 202, 'H': 275}, {'H': 88, 'O': 494, 'C': 151}]} """

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

