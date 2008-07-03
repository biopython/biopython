# Copyright 2008 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Bio.AlignIO support for the "nexus" file format.

You are expected to use this module via the Bio.AlignIO functions
(or the Bio.SeqIO functions).

See also the Bio.Nexus module (which this code calls internally),
as this offers more than just accessing the alignment or its
sequences as SeqRecord objects.
"""

from Bio.Nexus import Nexus
from Bio.Align.Generic import Alignment
from Bio.SeqRecord import SeqRecord

#You can get a couple of example files here:
#http://www.molecularevolution.org/resources/fileformats/
    
#This is a generator function!
def NexusIterator(handle, seq_count=None) :
    """Returns SeqRecord objects from a Nexus file.

    Thus uses the Bio.Nexus module to do the hard work.

    NOTE - We only expect ONE alignment matrix per Nexus file,
    meaning this iterator will only yield one Alignment."""
    n = Nexus.Nexus(handle)
    if not n.matrix :
        #No alignment found
        raise StopIteration
    alignment = Alignment(n.alphabet)

    #Bio.Nexus deals with duplicated names by adding a '.copy' suffix.
    #The original names and the modified names are kept in these two lists:
    assert len(n.unaltered_taxlabels) == len(n.taxlabels)
    
    if seq_count :
        assert seq_count == len(n.unaltered_taxlabels)
        
    for old_name, new_name in zip (n.unaltered_taxlabels, n.taxlabels) :
        assert new_name.startswith(old_name)
        seq = n.matrix[new_name] #already a Seq object with the alphabet set
        #ToDo - Can we extract any annotation too?
        #ToDo - Avoid abusing the private _records list
        alignment._records.append(SeqRecord(seq,
                                            id=new_name,
                                            name=old_name,
                                            description=""))
    #All done
    yield alignment

if __name__ == "__main__" :
    from StringIO import StringIO
    print "Quick self test"
    print
    print "Repeated names without a TAXA block"
    handle = StringIO("""#NEXUS
    [TITLE: NoName]

    begin data;
    dimensions ntax=4 nchar=50;
    format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";

    matrix
    CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---- 
    ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG 
    CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
    CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
    ;
    end; 
    """)
    for a in NexusIterator(handle) :
        print a
        for r in a :
            print repr(r.seq), r.name, r.id
    print "Done"

    print
    print "Repeated names with a TAXA block"
    handle = StringIO("""#NEXUS
    [TITLE: NoName]

    begin taxa
    CYS1_DICDI
    ALEU_HORVU
    CATH_HUMAN
    CYS1_DICDI;
    end;

    begin data;
    dimensions ntax=4 nchar=50;
    format interleave datatype=protein   gap=- symbols="FSTNKEYVQMCLAWPHDRIG";

    matrix
    CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---- 
    ALEU_HORVU          MAHARVLLLA LAVLATAAVA VASSSSFADS NPIRPVTDRA ASTLESAVLG 
    CATH_HUMAN          ------MWAT LPLLCAGAWL LGV------- -PVCGAAELS VNSLEK----
    CYS1_DICDI          -----MKVIL LFVLAVFTVF VSS------- --------RG IPPEEQ---X
    ;
    end; 
    """)
    for a in NexusIterator(handle) :
        print a
        for r in a :
            print repr(r.seq), r.name, r.id
    print "Done"
    print
    print "Reading an empty file"
    assert 0 == len(list(NexusIterator(StringIO())))
    print "Done"
