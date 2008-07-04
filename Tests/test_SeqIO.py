# Copyright 2007-2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from sets import Set
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from StringIO import StringIO
from Bio.Alphabet import generic_protein, generic_rna, generic_dna

#List of formats including alignment only file formats we can read AND write.
#The list is initially hard coded to preserve the original order of the unit
#test output, with any new formats added since appended to the end.
test_write_read_alignment_formats = ["fasta","clustal","phylip","stockholm"]
for format in SeqIO._FormatToWriter :
    if format not in test_write_read_alignment_formats :
        test_write_read_alignment_formats.append(format)
for format in AlignIO._FormatToWriter :
    if format not in test_write_read_alignment_formats :
        test_write_read_alignment_formats.append(format)

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
#Following example is also used in test_registry.py
    ("fasta",  False, 'Registry/seqs.fasta', 2), #contains blank line
#Following example is also used in test_Nexus.py
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
    ("swiss",  False, 'SwissProt/sp016', 1),
#Following example is also used in test_registry.py
    ("swiss",  False, 'Registry/EDD_RAT.dat', 1),
#Following examples are also used in test_GenBank.py
    ("genbank",False, 'GenBank/noref.gb', 1),
    ("genbank",False, 'GenBank/cor6_6.gb', 6),
    ("genbank",False, 'GenBank/iro.gb', 1),
    ("genbank",False, 'GenBank/pri1.gb', 1),
    ("genbank",False, 'GenBank/arab1.gb', 1),
    ("genbank",False, 'GenBank/protein_refseq.gb', 1), #Old version
    ("genbank",False, 'GenBank/protein_refseq2.gb', 1), #Revised version
    ("genbank",False, 'GenBank/extra_keywords.gb', 1),
    ("genbank",False, 'GenBank/one_of.gb', 1),
    ("genbank",False, 'GenBank/NT_019265.gb', 1),
    ("genbank",False, 'GenBank/origin_line.gb', 1),
    ("genbank",False, 'GenBank/blank_seq.gb', 1),
    ("genbank",False, 'GenBank/dbsource_wrap.gb', 1),
# The next example is a truncated copy of gbvrl1.seq from
# ftp://ftp.ncbi.nih.gov/genbank/gbvrl1.seq.gz
# This includes an NCBI header, and the first three records:
    ("genbank",False, 'GenBank/gbvrl1_start.seq', 3),
#Following files are also used in test_GFF.py
    ("genbank",False, 'GFF/NC_001422.gbk', 1),
#Following files are currently only used here:
    ("embl",      False, 'EMBL/TRBG361.embl', 1),
    ("embl",      False, 'EMBL/DD231055_edited.embl', 1),
    ("embl",      False, 'EMBL/SC10H5.embl', 1), # Pre 2006 style ID line
    ("embl",      False, 'EMBL/U87107.embl', 1), # Old ID line with SV line
    ("embl",      False, 'EMBL/AAA03323.embl', 1), # 2008, PA line but no AC
    ("stockholm", True,  'Stockholm/simple.sth', 2),
    ("stockholm", True,  'Stockholm/funny.sth', 5),
#Following PHYLIP files are currently only used here and in test_AlignIO.py,
#and are mostly from Joseph Felsenstein's PHYLIP v3.6 documentation:
    ("phylip",    True,  'Phylip/reference_dna.phy', 6),
    ("phylip",    True,  'Phylip/reference_dna2.phy', 6),
    ("phylip",    True,  'Phylip/hennigian.phy', 10),
    ("phylip",    True,  'Phylip/horses.phy', 10),
    ("phylip",    True,  'Phylip/random.phy', 10),
    ("phylip",    True,  'Phylip/interlaced.phy', 3),
    ("phylip",    True,  'Phylip/interlaced2.phy', 4),
#Following are EMBOSS simple or pairs format alignments
    ("emboss",    True,  'Emboss/alignret.txt', 4),
    ("emboss",    False, 'Emboss/needle.txt', 10),
    ("emboss",    True,  'Emboss/water.txt', 2),
#Following PHD (PHRAP) sequencing files are also used in test_Phd.py
    ("phd",       False, 'Phd/phd1', 3),
    ("phd",       False, 'Phd/phd2', 1),
#Following ACE assembly files are also used in test_Phd.py
    ("ace",       False, 'Ace/contig1.ace', 2),
    ("ace",       False, 'Ace/consed_sample.ace', 1),
    ("ace",       False, 'Ace/seq.cap.ace', 1),
#Following IntelliGenetics / MASE files are also used in test_intelligenetics.py
    ("ig",        False, 'IntelliGenetics/TAT_mase_nuc.txt', 17),
    ("ig",        True,  'IntelliGenetics/VIF_mase-pro.txt', 16),
    #This next file is a MASE alignment but sequence O_ANT70 is shorter than
    #the others (so the to_alignment() call will fail).  Perhaps MASE doesn't
    #write trailing gaps?
    ("ig",        False,  'IntelliGenetics/vpu_nucaligned.txt', 9),
    ]

# This is a list of two-tuples.  Each tuple contains a
# list of SeqRecord objects and a description (string)
test_records = [
    ([], "zero records"),
    ([SeqRecord(Seq("CHSMAIKLSSEHNIPSGIANAL",generic_protein), id="Alpha"),
      SeqRecord(Seq("HNGFTALEGEIHHLTHGEKVAF",generic_protein), id="Gamma"),
      SeqRecord(Seq("DITHGVG",generic_protein), id="delta")],
     "three peptides of different lengths"),
    ([SeqRecord(Seq("CHSMAIKLSSEHNIPSGIANAL",generic_protein), id="Alpha"),
      SeqRecord(Seq("VHGMAHPLGAFYNTPHGVANAI",generic_protein), id="Beta"),
      SeqRecord(Seq("HNGFTALEGEIHHLTHGEKVAF",generic_protein), id="Gamma")],
     "three proteins alignment"),
    ([SeqRecord(Seq("AATAAACCTTGCTGGCCATTGTGATCCATCCA",generic_dna), id="X"),
      SeqRecord(Seq("ACTCAACCTTGCTGGTCATTGTGACCCCAGCA",generic_dna), id="Y"),
      SeqRecord(Seq("TTTCCTCGGAGGCCAATCTGGATCAAGACCAT",generic_dna), id="Z")],
     "three DNA sequence alignment"),
    ([SeqRecord(Seq("AATAAACCTTGCTGGCCATTGTGATCCATCCA",generic_dna), id="X",
                name="The\nMystery\rSequece:%sX" % os.linesep),
      SeqRecord(Seq("ACTCAACCTTGCTGGTCATTGTGACCCCAGCA",generic_dna), id="Y",
                description="an%sevil\rdescription right\nhere" % os.linesep),
      SeqRecord(Seq("TTTCCTCGGAGGCCAATCTGGATCAAGACCAT",generic_dna), id="Z")],
     "3 DNA seq alignment with CR/LF in name/descr"),
    ([SeqRecord(Seq("CHSMAIKLSSEHNIPSGIANAL",generic_protein), id="Alpha"),
      SeqRecord(Seq("VHGMAHPLGAFYNTPHGVANAI",generic_protein), id="Beta"),
      SeqRecord(Seq("VHGMAHPLGAFYNTPHGVANAI",generic_protein), id="Beta"),
      SeqRecord(Seq("HNGFTALEGEIHHLTHGEKVAF",generic_protein), id="Gamma")],
     "alignment with repeated record"),
    ]
# Meddle with the annotation too:
assert test_records[4][1] == "3 DNA seq alignment with CR/LF in name/descr"
# Add a list of strings,
test_records[4][0][2].annotations["note"] = ["Note%salso" % os.linesep \
                                    + "\r\nhas\n evil line\rbreaks!", "Wow"]
# Add a simple string
test_records[4][0][2].annotations["comment"] = "More%sof" % os.linesep \
                                          + "\r\nthese\n evil line\rbreaks!"
# Add a float too:
test_records[4][0][2].annotations["weight"] = 2.5

def records_match(record_one, record_two) :
    """This is meant to be a strict comparison for exact agreement"""
    if record_one.id <> record_two.id :
        return False
    if record_one.name <> record_two.name :
        return False
    if record_one.description <> record_two.description :
        return False
    if record_one.seq.tostring() <> record_two.seq.tostring() :
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


def check_simple_write_read(records, indent=" ") :
    #print indent+"Checking we can write and then read back these records"
    for format in test_write_read_alignment_formats :
        print indent+"Checking can write/read as '%s' format" % format
        
        #Going to write to a handle...
        handle = StringIO()
        
        try :
            SeqIO.write(sequences=records, handle=handle, format=format)
        except ValueError, e :
            #This is often expected to happen, for example when we try and
            #write sequences of different lengths to an alignment file.
            print indent+"Failed: %s" % str(e)
            #Carry on to the next format:
            continue

        handle.flush()
        handle.seek(0)
        #Now ready to read back from the handle...
        try :
            records2 = list(SeqIO.parse(handle=handle, format=format))
        except ValueError, e :
            #This is BAD.  We can't read our own output.
            #I want to see the output when called from the test harness,
            #run_tests.py (which can be funny about new lines on Windows)
            handle.seek(0)
            raise ValueError("%s\n\n%s\n\n%s" \
                              % (str(e), repr(handle.read()), repr(records)))

        assert len(records2) == t_count
        for r1, r2 in zip(records, records2) :
            #Check the bare minimum (ID and sequence) as
            #many formats can't store more than that.

            #Check the sequence
            assert r1.seq.tostring() == r2.seq.tostring()
            #Beware of different quirks and limitations in the
            #valid character sets and the identifier lengths!
            if format=="phylip" :
                assert r1.id.replace("[","").replace("]","")[:10] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format=="clustal" :
                assert r1.id.replace(" ","_")[:30] == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format=="stockholm" :
                assert r1.id.replace(" ","_") == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)
            elif format=="fasta" :
                assert r1.id.split()[0] == r2.id
            else :
                assert r1.id == r2.id, \
                       "'%s' vs '%s'" % (r1.id, r2.id)


#Check parsers can cope with an empty file
for t_format in SeqIO._FormatToIterator :
    handle = StringIO()
    records = list(SeqIO.parse(handle, t_format))
    assert len(records) == 0

for (t_format, t_alignment, t_filename, t_count) in test_files :
    print "Testing reading %s format file %s" % (t_format, t_filename)
    assert os.path.isfile(t_filename), t_filename

    #Try as an iterator using handle
    records  = list(SeqIO.parse(handle=open(t_filename,"r"), format=t_format))
    assert len(records)  == t_count, \
         "Found %i records but expected %i" % (len(records), t_count)

    #Try using the iterator with a for loop
    records2 = []
    for record in SeqIO.parse(handle=open(t_filename,"r"), format=t_format) :
        records2.append(record)
    assert len(records2) == t_count

    #Try using the iterator with the next() method
    records3 = []
    seq_iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
    while True :
        try :
            record = seq_iterator.next()
        except StopIteration :
            record = None
        #Note that if the SeqRecord class has a __len__ method,
        #and it has a zero-length sequence, this would fail an
        #"if record" test.
        if record is not None :
            records3.append(record)
        else :
            break

    #Try a mixture of next() and list (a torture test!)
    seq_iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
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
    seq_iterator = SeqIO.parse(handle=open(t_filename,"r"), format=t_format)
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

        if "accessions" in record.annotations :
            accs = record.annotations["accessions"]
            #Check for blanks, or entries with leading/trailing spaces
            for acc in accs :
                assert acc and acc == acc.strip(), \
                    "Bad accession in annotations: %s" % repr(acc)
            assert len(Set(accs)) == len(accs), \
                   "Repeated accession in annotations: %s" % repr(accs)
        for ref in record.dbxrefs :
            assert ref and ref == ref.strip(), \
                "Bad cross reference in dbxrefs: %s" % repr(ref)
        assert len(record.dbxrefs) == len(record.dbxrefs), \
               "Repeated cross reference in dbxrefs: %s" % repr(record.dbxrefs)
                
            
        #Check the lists obtained by the different methods agree
        assert records_match(record, records2[i])
        assert records_match(record, records3[i])
        assert records_match(record, records4[i])
        assert records_match(record, records5[i])

        if i < 3 :
            print record_summary(record)
    # Only printed the only first three records: 0,1,2 
    if t_count > 4 :
        print " ..."
    if t_count > 3 :
        print record_summary(records[-1])

    # Check Bio.SeqIO.read(...)
    if t_count == 1 :
        record = SeqIO.read(handle=open(t_filename), format=t_format)
        assert isinstance(record, SeqRecord)
    else :
        try :
            record = SeqIO.read(open(t_filename), t_format)
            assert False, "Bio.SeqIO.read(...) should have failed"
        except ValueError :
            #Expected to fail
            pass


    if t_alignment :
        print "Testing reading %s format file %s as an alignment" \
              % (t_format, t_filename)

        #Using SeqIO.to_alignment(SeqIO.parse(...))
        alignment = SeqIO.to_alignment(SeqIO.parse( \
                    handle=open(t_filename,"r"), format=t_format))
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
    check_simple_write_read(records)

print "Finished tested reading files"
print
print "Starting testing writing records"
print "(Note that some of these are expected to 'fail' and say why)"
print
for (records, descr) in test_records :
    print "Testing can write/read %s" % descr
    for format in test_write_read_alignment_formats :
        print " Checking can write/read as '%s' format" % format

        #################
        # Write records #
        #################
        handle = StringIO()
        try :
            SeqIO.write(records, handle, format)
        except ValueError, e:
            #This is expected to happen on many of the examples.
            print " Failed: %s" % str(e)
            continue #goto next test

        #################
        # Read records  #
        #################
        handle.seek(0)
        try :
            new_records = list(SeqIO.parse(handle, format))
        except ValueError, e :
            #THIS INDICATES A SIGNIFICANT PROBLEM,
            #as we can't read the file we just wrote!
            print " FAILED: %s" % str(e)
            continue #goto next test

        #################
        # Check records #
        #################
        assert len(new_records) == len(records)
        for record, new_record in zip(records, new_records) :
            assert record.id == new_record.id
            assert record.seq.tostring() == new_record.seq.tostring()
            #Using records_match(record, new_record) is too strict

        #Close now, after checking, so that it can be used at the console for debugging
        handle.close()
        
print "Finished tested writing files"
