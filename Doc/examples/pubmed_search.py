"""Example script showing how to interact with PubMed."""
# standard library
import string

# biopython
from Bio import PubMed
from Bio import Medline

# do the search and get the ids
search_term = 'orchid'
orchid_ids = PubMed.search_for(search_term)

print orchid_ids

# access Medline through a dictionary interface that returns PubMed Records
rec_parser = Medline.RecordParser()
medline_dict = PubMed.Dictionary(parser = rec_parser)

for id in orchid_ids[0:5]:
    cur_record = medline_dict[id]
    print 'title:', string.rstrip(cur_record.title)
    print 'authors:', cur_record.authors
    print 'source:', string.strip(cur_record.source)
    print
     
