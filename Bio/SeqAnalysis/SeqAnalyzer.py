import pandas as pd
from Bio import SeqAnalysis
import numpy as np

seq_db = SeqAnalysis.database(['YP_025292.1', '1JOY'])
records = seq_db.get()








class analyzer(records):
    """
    A class that takes in records from SeqAnalysis.database and performs
    multiple operations on it and returns pd.DataFrame with values
    """
    def standard_deviaton(self):
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
        self.stdev = stdev


    def __init__(self):
        self.dict = analyzer.makedict(records) #creates a dict with seq ids and seqeances itself
        self.stdev = 0
        self.median = 0
        self.average = 0

    def results(self):
        data = {
            'id': [],
            'sequence_length': [],
            'molecular_weight': [],
            'theoretical_pl': [1.4, 24.5],
            'amino_acid_comp': [{'A': 12, 'G': 24}, {'L': 124, 'A': 52}],
            'atomic_comp': [{'C': 42, 'H': 335, 'N': 55}, {'C': 55, 'H': 523, 'N': 13, 'O': 53}],
            'extinction_coefficient': [1],
            'instability_index': [],
            'gravy': [],
            'est_half_life': [],
            'c_term_freq': [{'A': 0.2, 'G': 0.8}, {'A': 0.3, 'L': 0.7}],
            'n_term_freq': [{'A': 0.8, 'G': 0.2}, {'A': 0.7, 'L': 0.3}],
            'standard_deviation': self.stdev,
            'median':
        }
        df = pd.DataFrame(data)
        return df

seq_analyzer = SeqAnalysis.analyzer(records)
results = seq_analyzer.results()
