# Copyright 2008-2010 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Bio.AlignIO support for the "nexus" file format.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

See also the Bio.Nexus module (which this code calls internally),
as this offers more than just accessing the alignment or its
sequences as SeqRecord objects.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord  
from Bio.Nexus import Nexus
from Bio.Align import MultipleSeqAlignment
from Interfaces import AlignmentWriter
from Bio import Alphabet

#You can get a couple of example files here:
#http://www.molecularevolution.org/resources/fileformats/
    
#This is a generator function!
def NexusIterator(handle, seq_count=None):
    """Returns SeqRecord objects from a Nexus file.

    Thus uses the Bio.Nexus module to do the hard work.

    You are expected to call this function via Bio.SeqIO or Bio.AlignIO
    (and not use it directly).

    NOTE - We only expect ONE alignment matrix per Nexus file,
    meaning this iterator will only yield one MultipleSeqAlignment.
    """
    n = Nexus.Nexus(handle)
    if not n.matrix:
        #No alignment found
        raise StopIteration

    #Bio.Nexus deals with duplicated names by adding a '.copy' suffix.
    #The original names and the modified names are kept in these two lists:
    assert len(n.unaltered_taxlabels) == len(n.taxlabels)
    
    if seq_count and seq_count != len(n.unaltered_taxlabels):
        raise ValueError("Found %i sequences, but seq_count=%i" \
               % (len(n.unaltered_taxlabels), seq_count))

    #ToDo - Can we extract any annotation too?
    records = (SeqRecord(n.matrix[new_name], id=new_name, \
                         name=old_name, description="") \
               for old_name, new_name \
               in zip (n.unaltered_taxlabels, n.taxlabels))
    #All done
    yield MultipleSeqAlignment(records, n.alphabet)

class NexusWriter(AlignmentWriter):
    """Nexus alignment writer.

    Note that Nexus files are only expected to hold ONE alignment
    matrix.

    You are expected to call this class via the Bio.AlignIO.write() or
    Bio.SeqIO.write() functions.
    """
    def write_file(self, alignments):
        """Use this to write an entire file containing the given alignments.

        alignments - A list or iterator returning MultipleSeqAlignment objects.
                     This should hold ONE and only one alignment.
        """
        align_iter = iter(alignments) #Could have been a list
        try:
            first_alignment = align_iter.next()
        except StopIteration:
            first_alignment = None
        if first_alignment is None:
            #Nothing to write!
            return 0
        
        #Check there is only one alignment...
        try:
            second_alignment = align_iter.next()
        except StopIteration:
            second_alignment = None
        if second_alignment is not None:
            raise ValueError("We can only write one Alignment to a Nexus file.")

        #Good.  Actually write the single alignment,
        self.write_alignment(first_alignment)
        return 1 #we only support writing one alignment!

    def write_alignment(self, alignment):
        #Creates an empty Nexus object, adds the sequences,
        #and then gets Nexus to prepare the output.
        if len(alignment) == 0:
            raise ValueError("Must have at least one sequence")
        if alignment.get_alignment_length() == 0:
            raise ValueError("Non-empty sequences are required")
        minimal_record = "#NEXUS\nbegin data; dimensions ntax=0 nchar=0; " \
                         + "format datatype=%s; end;"  \
                         % self._classify_alphabet_for_nexus(alignment._alphabet)
        n = Nexus.Nexus(minimal_record)
        n.alphabet = alignment._alphabet
        for record in alignment:
            n.add_sequence(record.id, record.seq.tostring())
        n.write_nexus_data(self.handle)
    
    def _classify_alphabet_for_nexus(self, alphabet):
        """Returns 'protein', 'dna', 'rna' based on the alphabet (PRIVATE).

        Raises an exception if this is not possible."""
        #Get the base alphabet (underneath any Gapped or StopCodon encoding)
        a = Alphabet._get_base_alphabet(alphabet)

        if not isinstance(a, Alphabet.Alphabet):
            raise TypeError("Invalid alphabet")
        elif isinstance(a, Alphabet.ProteinAlphabet):
            return "protein"
        elif isinstance(a, Alphabet.DNAAlphabet):
            return "dna"
        elif isinstance(a, Alphabet.RNAAlphabet):
            return "rna"
        else:
            #Must be something like NucleotideAlphabet or
            #just the generic Alphabet (default for fasta files)
            raise ValueError("Need a DNA, RNA or Protein alphabet")

if __name__ == "__main__":
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
    for a in NexusIterator(handle):
        print a
        for r in a:
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
    for a in NexusIterator(handle):
        print a
        for r in a:
            print repr(r.seq), r.name, r.id
    print "Done"
    print
    print "Reading an empty file"
    assert 0 == len(list(NexusIterator(StringIO())))
    print "Done"
    print
    print "Writing..."
    
    handle = StringIO()
    NexusWriter(handle).write_file([a])
    handle.seek(0)
    print handle.read()

    handle = StringIO()
    try:
        NexusWriter(handle).write_file([a,a])
        assert False, "Should have rejected more than one alignment!"
    except ValueError:
        pass
