# Copyright 2007-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os

from Bio.SeqUtils import GC, quick_FASTA_reader
from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid
from Bio.SeqUtils.lcc import lcc_simp, lcc_mult
from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import single_letter_alphabet
from Bio import SeqIO


######################
# quick_FASTA_reader #
######################

dna_fasta_filename = "Fasta/f002"

tuple_records = quick_FASTA_reader(dna_fasta_filename)
assert len(tuple_records)==3
seq_records = list(SeqIO.parse(open(dna_fasta_filename),"fasta"))
assert len(seq_records)==3
for tuple_record, seq_record in zip(tuple_records, seq_records):
    assert tuple_record == (seq_record.description, seq_record.seq.tostring())
    print "%s has GC%% of %0.1f" % (seq_record.name, GC(seq_record.seq))

##############
# CodonUsage #
##############

print
print "Codon Adaption Index (CAI)"
CAI = CodonAdaptationIndex()
# Note - this needs a whole number of codons, and a DNA seq AS A STRING.
print "Example CAI %0.5f using E. coli (default)" \
      % CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG")

#We need a FASTA file of CDS sequences to count the codon usage...
dna_fasta_filename = "fasta.tmp"
dna_genbank_filename = "GenBank/NC_005816.gb"
record = SeqIO.read(open(dna_genbank_filename), "genbank")
records = []
for feature in record.features:
    if feature.type == "CDS" \
    and not feature.sub_features:
        start = feature.location.start.position
        end = feature.location.end.position
        table = int(feature.qualifiers["transl_table"][0])
        if feature.strand == -1:
            seq = record.seq[start:end].reverse_complement()
        else:
            seq = record.seq[start:end]
        #Double check we have the CDS sequence expected
        #TODO - Use any cds_start option if/when added to deal with the met
        a = "M" + str(seq[3:].translate(table))
        b = feature.qualifiers["translation"][0]+"*"
        assert a == b, "%r vs %r" % (a,b)
        records.append(SeqRecord(seq, id=feature.qualifiers["protein_id"][0],
                                 description=feature.qualifiers["product"][0]))
del start, end, table, seq
if os.path.isfile(dna_fasta_filename):
    os.remove(dna_fasta_filename)
handle = open(dna_fasta_filename, "w")
SeqIO.write(records, handle, "fasta")
handle.close()

CAI = CodonAdaptationIndex()
# Note - this needs a FASTA file which containing non-ambiguous DNA coding
# sequences - which should each be a whole number of codons.
CAI.generate_index(dna_fasta_filename)
print "Example CAI %0.5f using %s" \
      % (CAI.cai_for_gene("ATGCGTATCGATCGCGATACGATTAGGCGGATG"),
         record.annotations["source"])

os.remove(dna_fasta_filename)
del record, records
del dna_genbank_filename
del dna_fasta_filename

print

###################
# crc64 collision #
###################

# Example of crc64 collision from Sebastian Bassi using the
# immunoglobulin lambda light chain variable region from Homo sapiens
# Both sequences share the same CRC64 checksum: 44CAAD88706CC153
str_light_chain_one = "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGK" \
                + "APKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADY" \
                + "YCSSYAGSSTLVFGGGTKLTVL"
str_light_chain_two = "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGK" \
                + "APKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADY" \
                + "YCCSYAGSSTWVFGGGTKLTVL"

#Explicit testing of crc64 collision:
assert str_light_chain_one != str_light_chain_two
assert crc32(str_light_chain_one) != crc32(str_light_chain_two)
assert crc64(str_light_chain_one) == crc64(str_light_chain_two)
assert gcg(str_light_chain_one) != gcg(str_light_chain_two)
assert seguid(str_light_chain_one) != seguid(str_light_chain_two)

###########################
# main checksum/LCC tests #
###########################

#Print some output, which the test harness will check
examples = [str_light_chain_one, str_light_chain_two,
            "ATGCGTATCGATCGCGATACGATTAGGCGGAT"]

def u_crc32(seq):
    #NOTE - On Python 2 crc32 could return a signed int, but on Python 3 it is
    #always unsigned
    #Docs suggest should use crc32(x) & 0xffffffff for consistency.
    return crc32(seq) & 0xffffffff 

for i, seq_str in enumerate(examples):
    print "Example %i, length %i, %s..." % (i+1, len(seq_str), seq_str[:10])

    #Avoid cross platforms with printing floats by doing conversion explicitly
    def simple_LCC(s):
        return "%0.2f" % lcc_simp(s)

    def windowed_LCC(s):
        return ", ".join(["%0.2f" % v for v in lcc_mult(s,20)])

    for checksum in [u_crc32, crc64, gcg, seguid, simple_LCC, windowed_LCC]:
        #First using a string:
        value = checksum(seq_str)
        print " %s = %s" % (checksum.__name__, value)
        #Secondly check it works with a Seq object
        assert value == checksum(Seq(seq_str, single_letter_alphabet))
        #Finally check it works with a MutableSeq object
        assert value == checksum(MutableSeq(seq_str, single_letter_alphabet))

