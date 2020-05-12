import pandas as pd
#from Bio import SeqAnalysis

seq_db = SeqAnalysis.database(['YP_025292.1', '1JOY'])
records = seq_db.get()

class analyzer(records):

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
        pass



    def __init__(self):
        self.dict = analyzer.makedict(records) #creates a dict with seq ids and seqeances itself
        self.stdev = 





#creates a list of analyzer() objects
objects = analyzer.make_obj(dict)



seq_analyzer = SeqAnalysis.analyzer(records)

