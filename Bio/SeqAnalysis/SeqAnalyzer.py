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

    def makedict(records: list[tuple]) -> dict :
        """
        :param records: takes in a list of tuples, where tuple[0] = seq id, and tuple[1] = SeqRecord
        :return:
        returns a dict out of records, where id is the key and its sequence its value
        """
        dict = {}
        for _ in records:
            dict[_[0]] = (_[1].seq)[0]
        return dict


    def standard_deviaton(self):
        """
        Calculates standard deviation of lenths of given sequences
        :return:
        changes self.stdev values to standard value calculate from sequences
        """
        list_of_seq = dict.values()
        list_of_lenths = []
        for _ in list_of_seq:
            list_of_lenths.append(len(_))
        stdev = np.std(list_of_lenths)
        self.stdev = stdev


    def __init__(self):
        self.dict = analyzer.makedict(records) #creates a dict with seq ids and seqeances itself
        self.stdev = 0





#creates a list of analyzer() objects
objects = analyzer.make_obj(dict)



seq_analyzer = SeqAnalysis.analyzer(records)

