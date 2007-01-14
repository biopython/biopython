# Copyright 2006, 2007 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Clustal parsing based on earlier code by Thomas Sicheritz-Ponten,
# copyright 2001.  See Bio/SeqIO/generic.py

#For reading alignments:
from Bio.Alphabet import generic_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#For writing alignments:
from Bio.SeqIO.Interfaces import SequenceWriter
from Bio.Clustalw import ClustalAlignment

#This is a generator function!
def ClustalIterator(handle, alphabet = generic_alphabet) :
    """Reads a Clustalw file returning a SeqRecord object iterator

    The entire file is loaded at once, but the SeqRecord objects
    are only created "on request".

    For more information on the file format, please see:
    http://www.bioperl.org/wiki/ClustalW_multiple_alignment_format

    You might like to look at Bio.Clustalw which has an interface
    to the command line tool clustalw, and can also clustal alignment
    files into Bio.Clustalw.ClustalAlignment objects.
         
    We call this the "clustal" format which is consist with EMBOSS.
    Sadly BioPerl calls it the "clustalw" format, so we can't match
    them both.
    """
    line = handle.readline()
    if not line: return
    if not line[:7] == 'CLUSTAL':
        raise SyntaxError("Did not find CLUSTAL header")

    seqs = {}
    ids = []
    while True:
        line = handle.readline()
        if not line: break
        if line[0] == ' ': continue
        fields = line.rstrip().split()
        if not len(fields): continue
        
        #We expect there to be two fields, but on older files
        #there may be a third entry containing a letter count.
        if len(fields) < 2 or len(fields) > 3:
            raise SyntaxError("Could not parse line:\n%s" % line)

        name, seq = fields[0], fields[1]
        if not name in ids: ids.append(name)
        seqs.setdefault(name, '')
        seqs[name] += seq.upper()

        if len(fields) == 3 :
            #This MAY be an old style file with a letter count...
            try :
                oddity = int(fields[2])
            except ValueError :
                raise SyntaxError("Could not parse line, odd third field:\n%s" % line)
            #Check this equals the number of letters (excluding gaps)
            #so far for this sequence?
            if len(seqs[name].replace("-","")) <> oddity :
                raise SyntaxError("Could not parse line, odd third field:\n%s" % line)

    for id in ids :
        yield SeqRecord(Seq(seqs[id], alphabet), id=id)
    
class ClustalWriter(SequenceWriter):
    """Write Clustal sequence alignments"""
    def __init__(self, handle, truncate=10):
        """Creates the writer object

        Use the method write_file() to actually record your sequence records."""
        self.handle = handle
        self.truncate = truncate
    
    def write_file(self, records) :
        """Use this to write an entire file containing the given records.

        records - a SeqRecord iterator, or list of SeqRecords

        This code uses Bio.Clustalw.ClustalAlignment to do the hard
        work.  If you are working with alignment objects then using
        Bio.Clustalw.ClustalAlignment directly would be best.
        """
        # ToDo - decide if using Bio.Clustalw.ClustalAlignment is
        # actually the best way to handle this.
        #
        # Copying that thirty lines of code (with slight tweaks)
        # would be much simpler, and would probably run quicker and
        # use less memory as it doesn't build a ClustalAlignment
        # object.
        #
        # The downside is code duplication.
        alignment_length = None
        alignment = ClustalAlignment()
        for record in records :
            if alignment_length is None :
                alignment_length = len(record.seq)
            elif alignment_length <> len(record.seq) :
                raise ValueError, "Sequences of different lengths"
            
            #ToDo, check alphabet for this sequence matches that
            #specified for the alignment.  Not sure how the
            #alphabet.contains() method is intended to be used,
            #but it doesn't make sense to me right now.

            #Doing this works, but ClustalAlignment will use
            #the record.descrption when outputing the records.
            #alignment._records.append(record)
            alignment.add_sequence(record.id, record.seq.tostring())

        self.handle.write(str(alignment))
        self.handle.close()
