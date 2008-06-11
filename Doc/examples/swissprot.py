"""Example of connecting with exPASy and parsing SwissProt records."""

# biopython
from Bio import ExPASy, SwissProt

# 'O23729', 'O23730', 'O23731', Chalcone synthases from Orchid

ids = ['O23729', 'O23730', 'O23731']

for id in ids:
    handle = ExPASy.get_sprot_raw(id)
    record = SwissProt.read(handle)
    print "description:", record.description
    for ref in record.references:
        print "authors:", ref.authors
        print "title:", ref.title

    print "classification:", record.organism_classification
    print

