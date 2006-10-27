from Bio.Alphabet import generic_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#This is a generator function!
def ClustalIterator(handle, alphabet = generic_alphabet) :
    """Reads a Clustalw file returning a SeqRecord object iterator

    The entire file is loaded at once, but the SeqRecord objects
    are only created "on request".
    
    For more information on the file format, please see:
    http://www.bioperl.org/wiki/ClustalW_multiple_alignment_format
    """
    line = handle.readline()
    if not line: return
    if not line[:7] == 'CLUSTAL':
        raise SyntaxError("Did not find CLUSTAL header")
        #import sys
        #print >> sys.stderr, 'Warning file does not start with CLUSTAL header'

    seqs = {}
    ids = []
    while True:
        line = handle.readline()
        if not line: break
        if line[0] == ' ': continue
        fields = line.rstrip().split()
        if not len(fields): continue
        name, seq = fields
        if not name in ids: ids.append(name)
        seqs.setdefault(name, '')
        seqs[name] += seq.upper()

    for id in ids :
        yield SeqRecord(Seq(seqs[id], alphabet), id=id)
    
