# Copyright 2007 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from types import *
from Bio.SeqIO import *
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from StringIO import StringIO

#List of non-alignment file formats we can read AND write:
test_write_read_non_alignment_formats = ["fasta"]
#Longer list including alignment only file formats we can read AND write:
test_write_read_alignment_formats = test_write_read_non_alignment_formats[:]
test_write_read_alignment_formats.extend(["stockholm", "phylip"])

# test_files is a list of tuples containing:
# - string:  file format
# - boolean: alignment (requires all seqs be same length)
# - string:  relative filename
# - integer: number of sequences

test_files = [ \
#Following examples are also used in test_Clustalw.py
    ("clustal",True,  'Clustalw/cw02.aln', 2),
    ("clustal",True,  'Clustalw/opuntia.aln', 7),
#Following nucleic examples are also used in test_Fasta2.py
    ("fasta",  False, 'Nucleic/lupine.nu', 1),
    ("fasta",  False, 'Nucleic/elderberry.nu', 1),
    ("fasta",  False, 'Nucleic/phlox.nu', 1),
    ("fasta",  False, 'Nucleic/centaurea.nu', 1),
    ("fasta",  False, 'Nucleic/wisteria.nu', 1),
    ("fasta",  False, 'Nucleic/sweetpea.nu', 1),
    ("fasta",  False, 'Nucleic/lavender.nu', 1),
#Following protein examples are also used in test_Fasta2.py
    ("fasta",  False, 'Amino/aster.pro', 1),
    ("fasta",  False, 'Amino/loveliesbleeding.pro', 1),
    ("fasta",  False, 'Amino/rose.pro', 1),
    ("fasta",  False, 'Amino/rosemary.pro', 1),
#Following examples are also used in test_Fasta.py
    ("fasta",  False, 'Fasta/f001', 1), #Protein
    ("fasta",  False, 'Fasta/f002', 3), #DNA
    #("fasta", False, 'Fasta/f003', 2), #Protein with comments
    ("fasta",  False, 'Fasta/fa01', 2), #Protein with gaps
#Following examples are also used in test_GFF.py
    ("fasta",  False, 'GFF/NC_001802.fna', 1), #upper case
    ("fasta",  False, 'GFF/NC_001802lc.fna', 1), #lower case
    ("fasta",  True,  'GFF/multi.fna', 3), #Trivial nucleotide alignment
#Following examples are also used in test_Nexus.py (but not all of them?)
    ("nexus",  True,  'Nexus/f1.nex', 7),
    ("nexus",  True,  'Nexus/f2.nex', 9),
    ("nexus",  True,  'Nexus/test_Nexus_input.nex', 9),
#Following examples are also used in test_SwissProt.py
    ("swiss",  False, 'SwissProt/sp001', 1),
    ("swiss",  False, 'SwissProt/sp002', 1),
    ("swiss",  False, 'SwissProt/sp003', 1),
    ("swiss",  False, 'SwissProt/sp004', 1),
    ("swiss",  False, 'SwissProt/sp005', 1),
    ("swiss",  False, 'SwissProt/sp006', 1),
    ("swiss",  False, 'SwissProt/sp007', 1),
    ("swiss",  False, 'SwissProt/sp008', 1),
    ("swiss",  False, 'SwissProt/sp009', 1),
    ("swiss",  False, 'SwissProt/sp010', 1),
    ("swiss",  False, 'SwissProt/sp011', 1),
    ("swiss",  False, 'SwissProt/sp012', 1),
    ("swiss",  False, 'SwissProt/sp013', 1),
    ("swiss",  False, 'SwissProt/sp014', 1),
    ("swiss",  False, 'SwissProt/sp015', 1),
#Following examples are also used in test_GenBank.py
    ("genbank",False, 'GenBank/noref.gb', 1),
    ("genbank",False, 'GenBank/cor6_6.gb', 6),
    ("genbank",False, 'GenBank/iro.gb', 1),
    ("genbank",False, 'GenBank/pri1.gb', 1),
    ("genbank",False, 'GenBank/arab1.gb', 1),
    ("genbank",False, 'GenBank/protein_refseq.gb', 1),
    ("genbank",False, 'GenBank/extra_keywords.gb', 1),
    ("genbank",False, 'GenBank/one_of.gb', 1),
    ("genbank",False, 'GenBank/NT_019265.gb', 1),
    ("genbank",False, 'GenBank/origin_line.gb', 1),
    ("genbank",False, 'GenBank/blank_seq.gb', 1),
    ("genbank",False, 'GenBank/dbsource_wrap.gb', 1),
#Following files are also used in test_GFF.py
    ("genbank",False, 'GFF/NC_001422.gbk', 1),
#Following files are currently only used here:
    ("embl",      False, 'EMBL/TRBG361.embl', 1),
    ("stockholm", True,  'Stockholm/simple.sth', 2),
    ("stockholm", True,  'Stockholm/funny.sth', 5),
#Following PHYLIP files are currently only used here (test_SeqIO)
#and are mostly from Joseph Felsenstein's PHYLIP v3.6 documentation:
    ("phylip",    True,  'Phylip/reference_dna.phy', 6),
    ("phylip",    True,  'Phylip/reference_dna2.phy', 6),
    ("phylip",    True,  'Phylip/hennigian.phy', 10),
    ("phylip",    True,  'Phylip/horses.phy', 10),
    ("phylip",    True,  'Phylip/random.phy', 10),
    ("phylip",    True,  'Phylip/interlaced.phy', 3),
    ("phylip",    True,  'Phylip/interlaced2.phy', 4),
    ]

def records_match(record_one, record_two) :
    """This is meant to be a strict comparison for exact agreement"""
    if record_one.id <> record_two.id :
        return False
    if record_one.name <> record_two.name :
        return False
    if record_one.description <> record_two.description :
        return False
    #Close enough... should I check for features, annotation etc?
    return True

def record_summary(record, indent=" ") :
    """Returns a concise summary of a SeqRecord object as a string"""
    if record.id == record.name :
        answer = "%sID and Name='%s',\n%sSeq='" % (indent, record.id, indent)
    else :
        answer = "%sID = '%s', Name='%s',\n%sSeq='" % (indent, record.id, record.name, indent)
    if len(record.seq) > 50 :
        answer += record.seq[:40].tostring() + "..." + record.seq[-7:].tostring()
    else :
        answer += record.seq.tostring()
    answer += "', length=%i" % (len(record.seq))
    return answer

def col_summary(col_text) :
    if len(col_text) < 65 :
        return col_text
    else :
        return col_text[:60] + "..." + col_text[-5:]

def alignment_summary(alignment, index=" ") :
    """Returns a concise summary of an Alignment object as a string"""
    answer = []
    alignment_len = alignment.get_alignment_length()
    rec_count = len(alignment.get_all_seqs())
    for i in range(min(5,alignment_len)) :
        answer.append(index + col_summary(alignment.get_column(i)) \
                            + " alignment column %i" % i)
    if alignment_len > 5 :
        i = alignment_len - 1
        answer.append(index + col_summary("|" * rec_count) \
                            + " ...")
        answer.append(index + col_summary(alignment.get_column(i)) \
                            + " alignment column %i" % i)
    return "\n".join(answer)


def check_simple_write_read(records, formats, indent=" ") :
    #print indent+"Checking we can write and then read back these records"
    for format in formats :
        print indent+"Checking can write/read as '%s' format" % format
        
        #Going to write to a handle...
        handle = StringIO()
        WriteSequences(sequences=records, handle=handle, format=format)
        handle.flush()
        handle.seek(0)
        #Now ready to read back from the handle...
        records2 = list(SequenceIterator(handle=handle, format=format))

        assert len(records2) == t_count
        for i in range(t_count) :
            #Check the bare minimum (ID and sequence) as
            #many formats can't store more than that.
            #
            #Note some formats allow spaces, others don't.
            #Assume spaces are turned into underscores.
            records[i].id.replace(" ","_") == records2[i].id.replace(" ","_")
            records[i].seq.data == records2[i].seq.data

for (t_format, t_alignment, t_filename, t_count) in test_files :
    print "Testing reading %s format file %s" % (t_format, t_filename)
    assert os.path.isfile(t_filename)

    #Try as an iterator using handle
    records  = list(SequenceIterator(handle=open(t_filename,"rU"), format=t_format))
    assert len(records)  == t_count

    #Try using the iterator with a for loop
    records2 = []
    for record in SequenceIterator(handle=open(t_filename,"rU"), format=t_format) :
        records2.append(record)
    assert len(records2) == t_count

    #Try using the iterator with the next() method
    records3 = []
    seq_iterator = SequenceIterator(handle=open(t_filename,"rU"), format=t_format)
    while True :
        try :
            record = seq_iterator.next()
        except StopIteration :
            record = None
        if record :
            records3.append(record)
        else :
            break

    #Try a mixture of next() and list (a torture test!)
    seq_iterator = SequenceIterator(handle=open(t_filename,"rU"), format=t_format)
    try :
        record = seq_iterator.next()
    except StopIteration :
        record = None
    if record is not None :
        records4 = [record]
        records4.extend(list(seq_iterator))
    else :
        records4 = []
    assert len(records4) == t_count

    #Try a mixture of next() and for loop (a torture test!)
    seq_iterator = SequenceIterator(handle=open(t_filename,"rU"), format=t_format)
    try :
        record = seq_iterator.next()
    except StopIteration :
        record = None
    if record is not None :
        records5 = [record]
        for record in seq_iterator :
            records5.append(record)
    else :
        records5 = []
    assert len(records5) == t_count

    for i in range(t_count) :
        record = records[i]

        #Check returned expected object type
        assert isinstance(record, SeqRecord)
        assert isinstance(record.seq, Seq)
        assert isinstance(record.id, basestring)
        assert isinstance(record.name, basestring)
        assert isinstance(record.description, basestring)
        assert record.id <> ""

        #Check the lists obtained by the different methods agree
        assert records_match(record, records2[i])
        assert records_match(record, records3[i])
        assert records_match(record, records4[i])
        assert records_match(record, records5[i])

        if i < 10 :
            print record_summary(record)
    if t_count > 10 :
        print " ..."
        print record_summary(records[-1])


    if t_alignment :
        print "Testing reading %s format file %s as an alignment" \
              % (t_format, t_filename)

        #Using SequencesToAlignment(SequenceIterator(...))
        alignment = SequencesToAlignment(SequenceIterator( \
                    handle=open(t_filename,"rU"), format=t_format))
        assert len(alignment.get_all_seqs()) == t_count

        alignment_len = alignment.get_alignment_length()

        #Check the record order agrees, and double check the
        #sequence lengths all agree too.
        for i in range(t_count) :
            assert records_match(records[i], alignment.get_all_seqs()[i])
            assert len(records[i].seq) == alignment_len

        print alignment_summary(alignment)

    #Some alignment file formats have magic characters which mean
    #use the letter in this position in the first sequence.
    #They should all have been converted by the parser, but if
    #not reversing the record order might expose an error.  Maybe.
    records.reverse()
    if t_alignment :
        check_simple_write_read(records, test_write_read_alignment_formats)
    else :
        check_simple_write_read(records, test_write_read_non_alignment_formats)

print "Finished tested reading files with Bio.SeqIO"
