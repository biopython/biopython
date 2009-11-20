# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for fetching sequences, and parsing them

Uses Bio.SeqIO to parse files downloaded with Bio.GenBank, Bio.WWW.NCBI, 
Bio.ExPASy etc.

Goals:
    Make sure that all retrieval is working as expected.
    May catch some format changes early too.
"""
import sys
import requires_internet
requires_internet.check()

#We want to test these:
from Bio import Entrez
from Bio import ExPASy

#In order to check any sequences returned
from Bio import SeqIO
from StringIO import StringIO
from Bio.SeqUtils.CheckSum import seguid

#This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython-dev@biopython.org"

def checksum_summary(record):
    if len(record.seq) < 25:
        short = record.seq.tostring()
    else:
        short = record.seq.tostring()[:19] \
              + "..." + record.seq.tostring()[-3:]
    return "%s [%s] len %i" \
           % (short, seguid(record.seq), len(record.seq))

#####################################################################

print "Checking Bio.ExPASy.get_sprot_raw()"
id_list = ["O23729"]
for identifier in id_list:
    print "- Fetching %s" % identifier
    handle = ExPASy.get_sprot_raw(identifier)
    records = list(SeqIO.parse(handle, "swiss"))
    assert len(records)==1
    record = records[0]
    print "  Got " + checksum_summary(record)
    assert record.id == identifier
del id_list, handle, identifier, records, record

#####################################################################

print "Checking Bio.Entrez.efetch()"
for database, format, entry in [("genome","fasta","X52960"),
                                ("genome","gb","X52960"),
                                ("nucleotide", "fasta", "6273291"),
                                ("nucleotide", "gb", "6273291"),
                                ("protein", "fasta", "16130152"),
                                ("protein", "gb", "16130152")]:
    print "- Fetching %s from %s as %s" % (entry, database, format)
    handle = Entrez.efetch(db=database,
                           id=entry,
                           rettype=format)
    record = SeqIO.read(handle, format) # checks there is exactly one record
    handle.close()
    print "  Got " + checksum_summary(record)
    assert (entry in record.name) or (entry in record.id) \
        or ("gi" in record.annotations and record.annotations["gi"]==entry), \
           "%s got %s, %s" % (entry, record.name, record.id)
del database, format, entry, handle, record

#####################################################################
