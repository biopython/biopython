# Copyright 2008 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "phd" file format.

PHD files are output by PHRED and used by PHRAP and CONSED.

You are expected to use this module via the Bio.SeqIO functions.
See also the underlying Bio.Sequencing.Phd module."""

from Bio.SeqRecord import SeqRecord
from Bio.Sequencing import Phd
    
#This is a generator function!
def PhdIterator(handle) :
    """Returns SeqRecord objects from a PHD file.

    This uses the Bio.Sequencing.Phy module to do the hard work.
    """

    phd_records = Phd.parse(handle)
    for phd_record in phd_records:
        #Convert the PHY record into a SeqRecord...
        seq_record = SeqRecord(phd_record.seq,
                               id = phd_record.file_name,
                               name = phd_record.file_name)
        #Just re-use the comments dictionary as the SeqRecord's annotations
        seq_record.annotations = phd_record.comments
        yield seq_record 
    #All done

if __name__ == "__main__" :
    print "Quick self test"
    handle = open("../../Tests/Phd/Phd1")
    for record in PhdIterator(handle) :
        print record
    handle.close()
    print "Done"
        
    
