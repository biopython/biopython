"""Example of connecting with exPASy and parsing SwissProt records."""
# standard library
import string

# biopython
from Bio.WWW import ExPASy
from Bio.SwissProt import SProt
from Bio import File

# 'O23729', 'O23730', 'O23731', Chalcone synthases from Orchid

ids = ['O23729', 'O23730', 'O23731']

all_results = ''
for id in ids:
    results = ExPASy.get_sprot_raw(id)
    all_results = all_results + results.read()

s_parser = SProt.RecordParser()
s_iterator = SProt.Iterator(File.StringHandle(all_results), s_parser)

while 1:
    cur_record = s_iterator.next()

    if cur_record is None:
        break
    
    print "description:", cur_record.description
    for ref in cur_record.references:
        print "authors:", ref.authors
        print "title:", ref.title

    print "classification:", cur_record.organism_classification
    print

