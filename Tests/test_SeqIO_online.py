"""Testing online code for fetching sequences, and parsing them

Uses Bio.SeqIO to parse files downloaded with Bio.GenBank, Bio.WWW.NCBI, 
Bio.ExPASy etc.

Goals:
    Make sure that all retrieval is working as expected.
    May catch some format changes early too.
"""
import requires_internet

#We want to test these:
from Bio import GenBank
from Bio import ExPASy

#In order to check any sequences returned
from Bio import SeqIO
from StringIO import StringIO
from Bio.SeqUtils.CheckSum import seguid

#This lets us set the email address to be sent to NCBI Entrez:
from Bio import Entrez
Entrez.email = "biopython-dev@biopython.org"

def checksum_summary(record) :
    if len(record.seq) < 25 :
        short = record.seq.tostring()
    else :
        short = record.seq.tostring()[:19] \
              + "..." + record.seq.tostring()[-3:]
    return "%s [%s] len %i" \
           % (short, seguid(record.seq), len(record.seq))

#####################################################################

print "Checking Bio.ExPASy.get_sprot_raw()"
id_list = ["O23729"]
for identifier in id_list :
    print "- Fetching %s" % identifier
    handle = ExPASy.get_sprot_raw(identifier)
    records = list(SeqIO.parse(handle, "swiss"))
    assert len(records)==1
    record = records[0]
    print "  Got " + checksum_summary(record)
    assert record.id == identifier
del id_list, handle, identifier, records, record

#####################################################################

print "Checking Bio.GenBank.download_many()"
id_list = ['6273290', '6273289']
handle = GenBank.download_many(id_list)
for identifier, record in zip(id_list, SeqIO.parse(handle, "genbank")) :
    print "- Fetched GI %s as genbank" % identifier
    print "  Got " + checksum_summary(record)
    assert record.annotations["gi"] == identifier
handle.close()
del id_list, handle, identifier, record

#####################################################################

print "Checking Bio.GenBank.NCBIDictionary()"
for database, format, entry in [("genome","fasta","X52960"),
                                ("genome","genbank","X52960"),
                                ("nucleotide", "fasta", "6273291"),
                                ("nucleotide", "genbank", "6273291"),
                                ("protein", "fasta", "16130152"),
                                ("protein", "genbank", "16130152")] :
    print "- Fetching %s from %s as %s" % (entry, database, format)
    ncbi_dict = GenBank.NCBIDictionary(database, format)
    data = ncbi_dict[entry]
    assert isinstance(data, str)
    handle = StringIO(data) #turn into a handle
    records = list(SeqIO.parse(handle, format))
    handle.close()
    assert len(records) == 1
    print "  Got " + checksum_summary(records[0])
del database, format, entry, handle, data, records

#####################################################################
