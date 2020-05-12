import pandas as pd
from Bio import SeqAnalysis
import numpy as np

#seq_db = SeqAnalysis.database(['YP_025292.1', '1JOY'])
records = seq_db.get()


class analyzer(records):
    """
    A class that takes in records from SeqAnalysis.database and performs
    multiple operations on it and returns pd.DataFrame with values
    """

    def makedict(self) -> dict:
        """
        Makes a dictionary from records, where ids are keys and sequences are values
        :return:
        returns a created dictionary
        """
        slownik = {}
        for _ in records:
            slownik[_[0]] = _[1][0]
        return slownik

    def median_stdev_avg(self):
        """
        Calculates standard deviation of lenths of given sequences
        :return:
        changes self.stdev values to standard value calculate from sequences
        """
        list_of_seq = self.dict.values()
        list_of_lenths = []
        for _ in list_of_seq:
            list_of_lenths.append(len(_))
        stdev = np.std(list_of_lenths)
        avg = np.average(list_of_lenths)
        median = np.median(list_of_lenths)
        return stdev, avg, median

    def __init__(self):
        self.dict = analyzer.makedict(self)  # creates a dict with seq ids and seqeances itself
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


seq_analyzer = SeqAnalysis.analyzer(records)
results = seq_analyzer.results()
