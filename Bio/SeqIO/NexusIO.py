# Copyright 2006 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.SeqIO support for the "nexus" file format.

You were expected to use this module via the Bio.SeqIO functions.
This module has now been replaced by Bio.AlignIO.NexusIO, and is
deprecated."""

import warnings
warnings.warn("Bio.SeqIO.NexusIO is deprecated.  You can continue to read" \
              + " 'nexus' files with Bio.SeqIO, but this is now" \
              + " handled via Bio.AlignIO internally.",
              DeprecationWarning)

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus

#You can get a couple of example files here:
#http://www.molecularevolution.org/resources/fileformats/
    
#This is a generator function!
def NexusIterator(handle) :
    """Returns SeqRecord objects from a Nexus file.

    Thus uses the Bio.Nexus module to do the hard work."""
    n = Nexus.Nexus(handle)
    for id in n.original_taxon_order :
        if id in n.matrix :
            seq = n.matrix[id]
        else :
            #Missing the sequence?
            seq = Seq("", n.alphabet)
        #ToDo - Can we extract any annotation too?
        yield SeqRecord(seq, id=id, name=id, description="")
    #All done
    return
