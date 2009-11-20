# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio import UniGene

#Start of the UniGene file for Monodelphis domestica downloaded from:
#ftp://ftp.ncbi.nih.gov/repository/UniGene/Monodelphis_domestica
handle = open("UniGene/Mdm_partial.data")

ugparser = UniGene.Iterator(handle, UniGene.RecordParser())
for record in ugparser:
    assert isinstance(record.ID, str)
    assert isinstance(record.title, str)
    assert isinstance(record.species, str)
    assert isinstance(record.express, list)
    assert isinstance(record.sequence, list)

    print record.ID
    print "Title: '%s'" % record.title
    print "Expressed:", record.express
    print "Chromosome:", record.chromosome
    if record.sequence:
        print "Sequences:"
        for s in record.sequence:
            assert isinstance(s, UniGene.UnigeneSequenceRecord)
            print s
    else:
        print "No sequences"
            
    assert record.species == "Mdm"
    #Should be no PROTSIM lines in this file!
    assert isinstance(record.protsim, list)
    assert len(record.protsim) == 0

    print
print "Done"
handle.close()
