from Bio.SeqAnalysis.downloader import read_list_of_ids, SequenceDatabase

filename = "C:/Users/solak/Desktop/Docs/Pracka/list.txt"
ids = read_list_of_ids(filename)

s = SequenceDatabase(ids, "directory")
