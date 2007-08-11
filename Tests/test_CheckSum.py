# Copyright 2007 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.SeqUtils.CheckSum import crc32, crc64, gcg, seguid
from Bio.Seq import Seq
from Bio.Alphabet import single_letter_alphabet

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
assert str_light_chain_one <> str_light_chain_two
assert crc32(str_light_chain_one) <> crc32(str_light_chain_two)
assert crc64(str_light_chain_one) == crc64(str_light_chain_two)
assert gcg(str_light_chain_one) <> gcg(str_light_chain_two)
assert seguid(str_light_chain_one) <> seguid(str_light_chain_two)

##############
# Main tests #
##############

#Print some output, which the test harness will check
examples = [str_light_chain_one, str_light_chain_two,
            "ATGCGTATCGATCGCGATACGATTAGGCGGAT"]

for i, seq_str in enumerate(examples) :
    print "Example %i, length %i, %s..." % (i+1, len(seq_str), seq_str[:10])
    for checksum in [crc32, crc64, gcg, seguid] :
    
        #First using a string,
        seq_checksum = checksum(seq_str)
        print " %s = %s" % (checksum.__name__, seq_checksum)
        
        #Next using a Seq object
        seq_obj = Seq(seq_str, single_letter_alphabet)
        try :
            assert seq_checksum == checksum(seq_obj)
        except Exception, e:
   	        print " %s failed on Seq object, %s" % (checksum.__name__, str(e))
   	    

