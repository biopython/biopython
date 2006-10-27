from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Nexus import Nexus

#You can get a couple of example files here:
#http://www.molecularevolution.org/resources/fileformats/
    
#This is a generator function!
def NexusIterator(handle) :
    """Returns SeqRecord objects from a Nexus file

    Thus uses the Bio.Nexus module to do the hard work."""

    #Quick hack to get Bio.Nexus to cope with a StringIO from StringIO
    if not hasattr(handle,"name") :
        handle.name = None
        
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
