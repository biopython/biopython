# Copyright 2006 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "swiss" (aka SwissProt/UniProt) file format.

You are expected to use this module via the Bio.SeqIO functions.
See also the Bio.SwissProt module which offers more than just accessing
the sequences as SeqRecord objects."""

from Bio.SwissProt import SProt
import cStringIO
    
#This is a generator function!
def SwissIterator(handle) :
    """Breaks up a Swiss-Prot/UniProt file into SeqRecord objects.

    Every section from the ID line to the terminating // becomes
    a single SeqRecord with associated annotation and features.

    This parser is for the flat file "swiss" format as used by:
     * Swiss-Prot aka SwissProt
     * TrEMBL
     * UniProtKB aka UniProt Knowledgebase

    It does NOT read their new XML file format.
    http://www.expasy.org/sprot/

    For consistency with BioPerl and EMBOSS we call this the "swiss"
    format.
    """
    parser = SProt.SequenceParser()
    lines = []
    for line in handle:
        lines.append(line)
        if line[:2]=='//':
            handle = cStringIO.StringIO("".join(lines))
            record = parser.parse(handle)
            lines = []
            yield record
    #If there are more lines, it could only be a partial record.
    #Should we try and parse them anyway?
 

if __name__ == "__main__" :
    print "Quick self test..."

    example_filename = "../../Tests/SwissProt/sp008"

    import os
    if not os.path.isfile(example_filename):
        print "Missing test file %s" % example_filename
    else :
        #Try parsing it!
        handle = open(example_filename)
        records = SwissIterator(handle)
        for record in records:
            print record.name
            print record.id
            print record.annotations['keywords']
            print repr(record.annotations['organism'])
            print record.seq.tostring()[:20] + "..."
        handle.close()
