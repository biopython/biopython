# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from StringIO import StringIO
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Generic import Alignment
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

test_write_read_alignment_formats = sorted(AlignIO._FormatToWriter)
test_write_read_align_with_seq_count = test_write_read_alignment_formats \
                                     + ["fasta", "tab"]

# test_files is a list of tuples containing:
# - string:  file format
# - integer: number of sequences per alignment
# - integer: number of alignments
# - string:  relative filename
#
# Most of the input files are also used by test_SeqIO.py,
# and by other additional tests as noted below.
test_files = [
#Following examples are also used in test_Clustalw.py
    ("clustal", 2, 1, 'Clustalw/cw02.aln'),
    ("clustal", 7, 1, 'Clustalw/opuntia.aln'),
    ("clustal", 5, 1, 'Clustalw/hedgehog.aln'),
    ("clustal", 2, 1, 'Clustalw/odd_consensus.aln'),
    ("clustal",20, 1, 'Clustalw/protein.aln'), #Used in the tutorial
    ("clustal",20, 1, 'Clustalw/promals3d.aln'), #Nonstandard header
#Following examples are also used in test_GFF.py
    ("fasta", 3, 1, 'GFF/multi.fna'), #Trivial nucleotide alignment
#Following example is also used in test_Nexus.py
    ("nexus", 9, 1, 'Nexus/test_Nexus_input.nex'),
    ("stockholm", 2, 1, 'Stockholm/simple.sth'),
    ("stockholm", 6, 1, 'Stockholm/funny.sth'),
    ("phylip", 6, 1, 'Phylip/reference_dna.phy'),
    ("phylip", 6, 1, 'Phylip/reference_dna2.phy'),
    ("phylip",10, 1, 'Phylip/hennigian.phy'),
    ("phylip",10, 1, 'Phylip/horses.phy'),
    ("phylip",10, 1, 'Phylip/random.phy'),
    ("phylip", 3, 1, 'Phylip/interlaced.phy'),
    ("phylip", 4, 1, 'Phylip/interlaced2.phy'),
    ("phylip-relaxed", 12, 1, 'ExtendedPhylip/primates.phyx'),
    ("phylip-sequential", 3, 1, 'Phylip/sequential.phy'),
    ("phylip-sequential", 4, 1, 'Phylip/sequential2.phy'),
    ("emboss", 4, 1, 'Emboss/alignret.txt'),
    ("emboss", 2, 5, 'Emboss/needle.txt'),
    ("emboss", 2, 1, 'Emboss/needle_asis.txt'),
    ("emboss", 2, 1, 'Emboss/water.txt'),
    ("emboss", 2, 1, 'Emboss/water2.txt'),
    ("emboss", 2, 1, 'Emboss/matcher_simple.txt'),
    ("emboss", 2, 5, 'Emboss/matcher_pair.txt'),
    ("fasta-m10", 2, 4, 'Fasta/output001.m10'),
    ("fasta-m10", 2, 6, 'Fasta/output002.m10'),
    ("fasta-m10", 2, 3, 'Fasta/output003.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output004.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output005.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output006.m10'),
    ("fasta-m10", 2, 9, 'Fasta/output007.m10'),
    ("fasta-m10", 2, 12,'Fasta/output008.m10'),
    ("ig", 16, 1, 'IntelliGenetics/VIF_mase-pro.txt'),
    ("pir", 2, 1,  'NBRF/clustalw.pir'),
    ]

def str_summary(text, max_len=40):
    if len(text) <= max_len:
        return text
    else:
        return text[:max_len-4] + "..." + text[-3:]

def alignment_summary(alignment, index="  ", vertical_threshold=5):
    """Returns a concise summary of an Alignment object as a string."""
    answer = []
    alignment_len = alignment.get_alignment_length()
    rec_count = len(alignment)
    if rec_count < vertical_threshold:
        #Show each sequence row horizontally
        for record in alignment:
            answer.append("%s%s %s" \
            % (index,str_summary(record.seq.tostring()),record.id))
    else:
        #Show each sequence row vertically
        for i in range(min(5,alignment_len)):
            answer.append(index + str_summary(alignment[:,i]) \
                                + " alignment column %i" % i)
        if alignment_len > 5:
            i = alignment_len - 1
            answer.append(index + str_summary("|" * rec_count) \
                                + " ...")
            answer.append(index + str_summary(alignment[:,i]) \
                                + " alignment column %i" % i)
    return "\n".join(answer)


def check_simple_write_read(alignments, indent=" "):
    #print indent+"Checking we can write and then read back these alignments"
    for format in test_write_read_align_with_seq_count:
        records_per_alignment = len(alignments[0])
        for a in alignments:
            if records_per_alignment != len(a):
                records_per_alignment = None
        #Can we expect this format to work?
        if not records_per_alignment \
        and format not in test_write_read_alignment_formats:
            continue

        print indent+"Checking can write/read as '%s' format" % format

        #Going to write to a handle...
        handle = StringIO()

        try:
            c = AlignIO.write(alignments, handle=handle, format=format)
            assert c == len(alignments)
        except ValueError, e:
            #This is often expected to happen, for example when we try and
            #write sequences of different lengths to an alignment file.
            print indent+"Failed: %s" % str(e)
            #Carry on to the next format:
            continue

        #First, try with the seq_count
        if records_per_alignment:
            handle.flush()
            handle.seek(0)
            try:
                alignments2 = list(AlignIO.parse(handle=handle, format=format, \
                                                 seq_count=records_per_alignment))
            except ValueError, e:
                #This is BAD.  We can't read our own output.
                #I want to see the output when called from the test harness,
                #run_tests.py (which can be funny about new lines on Windows)
                handle.seek(0)
                raise ValueError("%s\n\n%s\n\n%s" \
                                  % (str(e), repr(handle.read()), repr(alignments2)))
            simple_alignment_comparison(alignments, alignments2, format)

        if format in test_write_read_alignment_formats:
            #Don't need the seq_count
            handle.flush()
            handle.seek(0)
            try:
                alignments2 = list(AlignIO.parse(handle=handle, format=format))
            except ValueError, e:
                #This is BAD.  We can't read our own output.
                #I want to see the output when called from the test harness,
                #run_tests.py (which can be funny about new lines on Windows)
                handle.seek(0)
                raise ValueError("%s\n\n%s\n\n%s" \
                                  % (str(e), repr(handle.read()), repr(alignments2)))
            simple_alignment_comparison(alignments, alignments2, format)

        if len(alignments)>1:
            #Try writing just one Alignment (not a list)
            handle = StringIO()
            SeqIO.write(alignments[0], handle, format)
            assert handle.getvalue() == alignments[0].format(format)

def simple_alignment_comparison(alignments, alignments2, format):
    assert len(alignments) == len(alignments2)
    for a1, a2 in zip(alignments, alignments2):
        assert a1.get_alignment_length() == a2.get_alignment_length()
        assert len(a1) == len(a2)
        for r1, r2 in zip(a1,a2):
            #Check the bare minimum (ID and sequence) as
            #many formats can't store more than that.

            #Check the sequence
            assert r1.seq.tostring() == r2.seq.tostring()

            #Beware of different quirks and limitations in the
            #valid character sets and the identifier lengths!
            if format in ["phylip", "phylip-sequential"]:
                assert r1.id.replace("[","").replace("]","")[:10] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format=="phylip-relaxed":
                assert r1.id.replace(" ", "").replace(':', '|') == r2.id, \
                        "'%s' vs '%s'" % (r1.id, r2.id)
            elif format=="clustal":
                assert r1.id.replace(" ","_")[:30] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format=="stockholm":
                assert r1.id.replace(" ","_") == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format=="fasta":
                assert r1.id.split()[0] == r2.id
            else:
                assert r1.id == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
    return True

#Check Phylip files reject duplicate identifiers.
def check_phylip_reject_duplicate():
    """
    Ensure that attempting to write sequences with duplicate IDs after
    truncation fails for Phylip format.
    """
    handle = StringIO()
    sequences = [SeqRecord(Seq('AAAA'), id='longsequencename1'),
                 SeqRecord(Seq('AAAA'), id='longsequencename2'),
                 SeqRecord(Seq('AAAA'), id='other_sequence'),]
    alignment = MultipleSeqAlignment(sequences)
    try:
        # This should raise a ValueError
        AlignIO.write(alignment, handle, 'phylip')
        assert False, "Duplicate IDs after truncation are not allowed."
    except ValueError, e:
        # Expected - check the error
        assert "Repeated name 'longsequen'" in str(e)

check_phylip_reject_duplicate()


#Check parsers can cope with an empty file
for t_format in AlignIO._FormatToIterator:
     handle = StringIO()
     alignments = list(AlignIO.parse(handle, t_format))
     assert len(alignments) == 0

#Check writers can cope with no alignments
for t_format in list(AlignIO._FormatToWriter)+list(SeqIO._FormatToWriter):
     handle = StringIO()
     assert 0 == AlignIO.write([], handle, t_format), \
            "Writing no alignments to %s format should work!" \
            % t_format

#Check writers reject non-alignments
list_of_records = list(AlignIO.read(open("Clustalw/opuntia.aln"),"clustal"))
for t_format in list(AlignIO._FormatToWriter)+list(SeqIO._FormatToWriter):
    handle = StringIO()
    try:
        AlignIO.write([list_of_records], handle, t_format)
        assert False, "Writing non-alignment to %s format should fail!" \
            % t_format
    except (TypeError, AttributeError, ValueError):
        pass
    del handle
del list_of_records, t_format

#Main tests...
for (t_format, t_per, t_count, t_filename) in test_files:
    print "Testing reading %s format file %s with %i alignments" \
          % (t_format, t_filename, t_count)
    assert os.path.isfile(t_filename), t_filename

    #Try as an iterator using handle
    alignments  = list(AlignIO.parse(handle=open(t_filename,"r"), format=t_format))
    assert len(alignments)  == t_count, \
         "Found %i alignments but expected %i" % (len(alignments), t_count)
    for alignment in alignments:
        assert len(alignment) == t_per, \
            "Expected %i records per alignment, got %i" \
            % (t_per, len(alignment))

    #Try using the iterator with a for loop and a filename not handle
    alignments2 = []
    for record in AlignIO.parse(t_filename, format=t_format):
        alignments2.append(record)
    assert len(alignments2) == t_count

    #Try using the iterator with the next() method
    alignments3 = []
    seq_iterator = AlignIO.parse(handle=open(t_filename,"r"), format=t_format)
    while True:
        try:
            record = seq_iterator.next()
        except StopIteration:
            break
        assert record is not None, "Should raise StopIteration not return None"
        alignments3.append(record)

    #Try a mixture of next() and list (a torture test!)
    seq_iterator = AlignIO.parse(handle=open(t_filename,"r"), format=t_format)
    try:
        record = seq_iterator.next()
    except StopIteration:
        record = None
    if record is not None:
        alignments4 = [record]
        alignments4.extend(list(seq_iterator))
    else:
        alignments4 = []
    assert len(alignments4) == t_count

    #Try a mixture of next() and for loop (a torture test!)
    seq_iterator = AlignIO.parse(handle=open(t_filename,"r"), format=t_format)
    try:
        record = seq_iterator.next()
    except StopIteration:
        record = None
    if record is not None:
        alignments5 = [record]
        for record in seq_iterator:
            alignments5.append(record)
    else:
        alignments5 = []
    assert len(alignments5) == t_count

    # Check Bio.AlignIO.read(...)
    if t_count == 1:
        alignment = AlignIO.read(handle=open(t_filename), format=t_format)
        assert isinstance(alignment, Alignment)
    else:
        try:
            alignment = AlignIO.read(open(t_filename), t_format)
            assert False, "Bio.AlignIO.read(...) should have failed"
        except ValueError:
            #Expected to fail
            pass

    #Print the alignment
    for i,alignment in enumerate(alignments):
        if i < 3 or i+1 == t_count:
            print " Alignment %i, with %i sequences of length %i" \
                  % (i,
                     len(alignment),
                     alignment.get_alignment_length())
            print alignment_summary(alignment)
        elif i==3:
            print " ..."

    #Check AlignInfo.SummaryInfo likes the alignment
    summary = AlignInfo.SummaryInfo(alignment)
    dumb_consensus = summary.dumb_consensus()
    #gap_consensus = summary.gap_consensus()
    if t_format != "nexus":
        #Hack for bug 2535
        pssm = summary.pos_specific_score_matrix()
        rep_dict = summary.replacement_dictionary()
        try:
            info_content = summary.information_content()
        except ValueError, e:
            if str(e) != "Error in alphabet: not Nucleotide or Protein, supply expected frequencies":
                raise e
            pass

    if t_count==1 and t_format not in ["nexus","emboss","fasta-m10"]:
        #print " Trying to read a triple concatenation of the input file"
        data = open(t_filename,"r").read()
        handle = StringIO()
        handle.write(data + "\n\n" + data + "\n\n" + data)
        handle.seek(0)
        assert 3 == len(list(AlignIO.parse(handle=handle, format=t_format, seq_count=t_per)))

    #Some alignment file formats have magic characters which mean
    #use the letter in this position in the first sequence.
    #They should all have been converted by the parser, but if
    #not reversing the record order might expose an error.  Maybe.
    alignments.reverse()
    check_simple_write_read(alignments)

print "Finished tested reading files"
