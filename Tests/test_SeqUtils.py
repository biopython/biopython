# Copyright 2007-2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.SeqUtils import GC, quick_FASTA_reader
from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid
from Bio.SeqUtils.lcc import lcc_simp, lcc_mult
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import single_letter_alphabet
from Bio import SeqIO

######################
# quick_FASTA_reader #
######################

tuple_records = quick_FASTA_reader("Fasta/f002")
assert len(tuple_records)==3
seq_records = list(SeqIO.parse(open("Fasta/f002"),"fasta"))
assert len(seq_records)==3
for tuple_record, seq_record in zip(tuple_records, seq_records) :
    assert tuple_record == (seq_record.description, seq_record.seq.tostring())
    print "%s has GC%% of %0.1f" % (seq_record.name, GC(seq_record.seq))

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

for i, seq_str in enumerate(examples) :
    print "Example %i, length %i, %s..." % (i+1, len(seq_str), seq_str[:10])

    #Avoid cross platforms with printing floats by doing conversion explicitly
    def simple_LCC(s) :
        return "%0.2f" % lcc_simp(s)

    def windowed_LCC(s) :
        return ", ".join(["%0.2f" % v for v in lcc_mult(s,20)])

    for checksum in [crc32, crc64, gcg, seguid, simple_LCC, windowed_LCC] :
        #First using a string:
        value = checksum(seq_str)
        print " %s = %s" % (checksum.__name__, value)
        #Secondly check it works with a Seq object
        assert value == checksum(Seq(seq_str, single_letter_alphabet))
        #Finally check it works with a MutableSeq object
        assert value == checksum(MutableSeq(seq_str, single_letter_alphabet))
