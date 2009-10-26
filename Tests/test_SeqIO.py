# Copyright 2007-2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os

from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq
from StringIO import StringIO
from Bio import Alphabet

import warnings
def send_warnings_to_stdout(message, category, filename, lineno,
                                file=None, line=None):
    #TODO - Have Biopython DataLossWarning?
    if category in [UserWarning] :
        print "%s - %s" % (category.__name__, message)
warnings.resetwarnings()
warnings.showwarning = send_warnings_to_stdout


protein_alphas = [Alphabet.generic_protein]
dna_alphas = [Alphabet.generic_dna]
rna_alphas = [Alphabet.generic_rna]
nucleotide_alphas = [Alphabet.generic_nucleotide,
                     Alphabet.Gapped(Alphabet.generic_nucleotide)]
no_alpha_formats = ["fasta","clustal","phylip","tab","ig","stockholm","emboss",
                    "fastq","fastq-solexa","fastq-illumina","qual"]
possible_unknown_seq_formats = ["qual", "genbank", "gb", "embl"]

#List of formats including alignment only file formats we can read AND write.
#The list is initially hard coded to preserve the original order of the unit
#test output, with any new formats added since appended to the end.
test_write_read_alignment_formats = ["fasta","clustal","phylip","stockholm"]
for format in sorted(SeqIO._FormatToWriter) :
    if format not in test_write_read_alignment_formats :
        test_write_read_alignment_formats.append(format)
for format in sorted(AlignIO._FormatToWriter) :
    if format not in test_write_read_alignment_formats :
        test_write_read_alignment_formats.append(format)
test_write_read_alignment_formats.remove("gb") #an alias for genbank
test_write_read_alignment_formats.remove("fastq-sanger") #an alias for fastq

# test_files is a list of tuples containing:
# - string:  file format
# - boolean: alignment (requires all seqs be same length)
# - string:  relative filename
# - integer: number of sequences

test_files = [ \
    ("sff",    False, 'Roche/E3MFGYR02_random_10_reads.sff', 10),
#Following examples are also used in test_Clustalw.py
    ("clustal",True,  'Clustalw/cw02.aln', 2),
    ("clustal",True,  'Clustalw/opuntia.aln', 7),
    ("clustal",True,  'Clustalw/hedgehog.aln', 5),
    ("clustal",True,  'Clustalw/odd_consensus.aln', 2),
#Following nucleic examples are also used in test_SeqIO_FastaIO.py
    ("fasta",  False, 'Nucleic/lupine.nu', 1),
    ("fasta",  False, 'Nucleic/elderberry.nu', 1),
    ("fasta",  False, 'Nucleic/phlox.nu', 1),
    ("fasta",  False, 'Nucleic/centaurea.nu', 1),
    ("fasta",  False, 'Nucleic/wisteria.nu', 1),
    ("fasta",  False, 'Nucleic/sweetpea.nu', 1),
    ("fasta",  False, 'Nucleic/lavender.nu', 1),
#Following protein examples are also used in test_SeqIO_FastaIO.py
    ("fasta",  False, 'Amino/aster.pro', 1),
    ("fasta",  False, 'Amino/loveliesbleeding.pro', 1),
    ("fasta",  False, 'Amino/rose.pro', 1),
    ("fasta",  False, 'Amino/rosemary.pro', 1),
#Following examples are also used in test_Fasta.py
    ("fasta",  False, 'Fasta/f001', 1), #Protein
    ("fasta",  False, 'Fasta/f002', 3), #DNA
    #("fasta", False, 'Fasta/f003', 2), #Protein with comments
    ("fasta",  False, 'Fasta/fa01', 2), #Protein with gaps
#Following are also used in test_SeqIO_features.py, see also NC_005816.gb
    ("fasta",  False, 'GenBank/NC_005816.fna', 1),
    ("fasta",  False, 'GenBank/NC_005816.ffn', 10),
    ("fasta",  False, 'GenBank/NC_005816.faa', 10),
    ("fasta",  False, 'GenBank/NC_000932.faa', 85),
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
    ("genbank",False, 'GenBank/NC_005816.gb', 1), #See also AE017046.embl
    ("genbank",False, 'GenBank/NC_000932.gb', 1),
    ("genbank",False, 'GenBank/pBAD30.gb', 1), #Odd LOCUS line from Vector NTI
# The next example is a truncated copy of gbvrl1.seq from
# ftp://ftp.ncbi.nih.gov/genbank/gbvrl1.seq.gz
# This includes an NCBI header, and the first three records:
    ("genbank",False, 'GenBank/gbvrl1_start.seq', 3),
#Following files are also used in test_GFF.py
    ("genbank",False, 'GFF/NC_001422.gbk', 1),
#Following files are currently only used here:
    ("embl",   False, 'EMBL/TRBG361.embl', 1),
    ("embl",   False, 'EMBL/DD231055_edited.embl', 1),
    ("embl",   False, 'EMBL/SC10H5.embl', 1), # Pre 2006 style ID line
    ("embl",   False, 'EMBL/U87107.embl', 1), # Old ID line with SV line
    ("embl",   False, 'EMBL/AAA03323.embl', 1), # 2008, PA line but no AC
    ("embl",   False, 'EMBL/AE017046.embl', 1), #See also NC_005816.gb
    ("stockholm", True,  'Stockholm/simple.sth', 2),
    ("stockholm", True,  'Stockholm/funny.sth', 5),
#Following PHYLIP files are currently only used here and in test_AlignIO.py,
#and are mostly from Joseph Felsenstein's PHYLIP v3.6 documentation:
    ("phylip", True,  'Phylip/reference_dna.phy', 6),
    ("phylip", True,  'Phylip/reference_dna2.phy', 6),
    ("phylip", True,  'Phylip/hennigian.phy', 10),
    ("phylip", True,  'Phylip/horses.phy', 10),
    ("phylip", True,  'Phylip/random.phy', 10),
    ("phylip", True,  'Phylip/interlaced.phy', 3),
    ("phylip", True,  'Phylip/interlaced2.phy', 4),
#Following are EMBOSS simple or pairs format alignments
    ("emboss", True,  'Emboss/alignret.txt', 4),
    ("emboss", False, 'Emboss/needle.txt', 10),
    ("emboss", True,  'Emboss/water.txt', 2),
#Following PHD (PHRAP) sequencing files are also used in test_Phd.py
    ("phd", False, 'Phd/phd1', 3),
    ("phd", False, 'Phd/phd2', 1),
    ("phd", False, 'Phd/phd_solexa', 2),
    ("phd", False, 'Phd/phd_454', 1),
#Following ACE assembly files are also used in test_Ace.py
    ("ace", False, 'Ace/contig1.ace', 2),
    ("ace", False, 'Ace/consed_sample.ace', 1),
    ("ace", False, 'Ace/seq.cap.ace', 1),
#Following IntelliGenetics / MASE files are also used in test_intelligenetics.py
    ("ig",  False, 'IntelliGenetics/TAT_mase_nuc.txt', 17),
    ("ig",  True,  'IntelliGenetics/VIF_mase-pro.txt', 16),
    #This next file is a MASE alignment but sequence O_ANT70 is shorter than
    #the others (so the to_alignment() call will fail).  Perhaps MASE doesn't
    #write trailing gaps?
    ("ig",  False,  'IntelliGenetics/vpu_nucaligned.txt', 9),
#Following NBRD-PIR files are used in test_nbrf.py
    ("pir", False, 'NBRF/B_nuc.pir', 444),
    ("pir", False, 'NBRF/Cw_prot.pir', 111),
    ("pir", False, 'NBRF/DMA_nuc.pir', 4),
    ("pir", False, 'NBRF/DMB_prot.pir', 6),
    ("pir", True,  'NBRF/clustalw.pir', 2),
#Following quality files are also used in the Bio.SeqIO.QualityIO doctests:
    ("fasta", True, 'Quality/example.fasta', 3),
    ("qual",  False,'Quality/example.qual',  3),
    ("fastq", True, 'Quality/example.fastq', 3),
    ("fastq", True, 'Quality/tricky.fastq', 4),
    ("fastq", False,'Quality/sanger_faked.fastq', 1),
    ("fastq", False,'Quality/sanger_93.fastq', 1),
    ("fastq-illumina", False,'Quality/illumina_faked.fastq', 1),
    ("fastq-solexa", False, 'Quality/solexa_faked.fastq', 1),
    ("fastq-solexa", True, 'Quality/solexa_example.fastq', 5),
    ]

# This is a list of two-tuples.  Each tuple contains a
# list of SeqRecord objects and a description (string)
test_records = [
    ([], "zero records"),
    ([SeqRecord(Seq("CHSMAIKLSSEHNIPSGIANAL",Alphabet.generic_protein), id="Alpha"),
      SeqRecord(Seq("HNGFTALEGEIHHLTHGEKVAF",Alphabet.generic_protein), id="Gamma"),
      SeqRecord(Seq("DITHGVG",Alphabet.generic_protein), id="delta")],
     "three peptides of different lengths"),
    ([SeqRecord(Seq("CHSMAIKLSSEHNIPSGIANAL",Alphabet.generic_protein), id="Alpha"),
      SeqRecord(Seq("VHGMAHPLGAFYNTPHGVANAI",Alphabet.generic_protein), id="Beta"),
      SeqRecord(Seq("HNGFTALEGEIHHLTHGEKVAF",Alphabet.generic_protein), id="Gamma")],
     "three proteins alignment"),
    ([SeqRecord(Seq("AATAAACCTTGCTGGCCATTGTGATCCATCCA",Alphabet.generic_dna), id="X"),
      SeqRecord(Seq("ACTCAACCTTGCTGGTCATTGTGACCCCAGCA",Alphabet.generic_dna), id="Y"),
      SeqRecord(Seq("TTTCCTCGGAGGCCAATCTGGATCAAGACCAT",Alphabet.generic_dna), id="Z")],
     "three DNA sequence alignment"),
    ([SeqRecord(Seq("AATAAACCTTGCTGGCCATTGTGATCCATCCA",Alphabet.generic_dna), id="X",
                name="The\nMystery\rSequece:\r\nX"),
      SeqRecord(Seq("ACTCAACCTTGCTGGTCATTGTGACCCCAGCA",Alphabet.generic_dna), id="Y",
                description="an%sevil\rdescription right\nhere" % os.linesep),
      SeqRecord(Seq("TTTCCTCGGAGGCCAATCTGGATCAAGACCAT",Alphabet.generic_dna), id="Z")],
     "3 DNA seq alignment with CR/LF in name/descr"),
    ([SeqRecord(Seq("CHSMAIKLSSEHNIPSGIANAL",Alphabet.generic_protein), id="Alpha"),
      SeqRecord(Seq("VHGMAHPLGAFYNTPHGVANAI",Alphabet.generic_protein), id="Beta"),
      SeqRecord(Seq("VHGMAHPLGAFYNTPHGVANAI",Alphabet.generic_protein), id="Beta"),
      SeqRecord(Seq("HNGFTALEGEIHHLTHGEKVAF",Alphabet.generic_protein), id="Gamma")],
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

def compare_record(record_one, record_two) :
    """This is meant to be a strict comparison for exact agreement..."""
    assert isinstance(record_one, SeqRecord)
    assert isinstance(record_two, SeqRecord)
    if record_one.id != record_two.id :
        return False
    if record_one.name != record_two.name :
        return False
    if record_one.description != record_two.description :
        return False
    if record_one.seq is not None and record_two.seq is not None \
    and record_one.seq.tostring() != record_two.seq.tostring() :
        return False
    #TODO - check features and annotation (see code for BioSQL tests)
    for key in set(record_one.letter_annotations).intersection( \
                   record_two.letter_annotations) :
        if record_one.letter_annotations[key] != \
           record_two.letter_annotations[key] :
            return False
    return True

def record_summary(record, indent=" ") :
    """Returns a concise summary of a SeqRecord object as a string"""
    if record.id == record.name :
        answer = "%sID and Name='%s',\n%sSeq='" % (indent, record.id, indent)
    else :
        answer = "%sID = '%s', Name='%s',\n%sSeq='" % (indent, record.id, record.name, indent)
    if record.seq is None :
        answer += "None"
    else :
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
        if format not in possible_unknown_seq_formats \
        and isinstance(records[0].seq, UnknownSeq) \
        and len(records[0].seq) > 100 :
           #Skipping for speed.  Some of the unknown sequences are
           #rather long, and it seems a bit pointless to record them.
           continue
        print indent+"Checking can write/read as '%s' format" % format
        
        #Going to write to a handle...
        handle = StringIO()
        
        try :
            c = SeqIO.write(sequences=records, handle=handle, format=format)
            assert c == len(records)
        except (TypeError, ValueError), e :
            #This is often expected to happen, for example when we try and
            #write sequences of different lengths to an alignment file.
            if "len()" in str(e) :
                #Python 2.4.3,
                #>>> len(None)
                #...
                #TypeError: len() of unsized object
                #
                #Python 2.5.2,
                #>>> len(None)
                #...
                #TypeError: object of type 'NoneType' has no len()
                print "Failed: Probably len() of None"
            else :
                print indent+"Failed: %s" % str(e)
            assert format != t_format, \
                   "Should be able to re-write in the original format!"
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
            if format in ["gb", "genbank"] :
                #The GenBank parser will convert everything to upper case.
                assert r1.seq.tostring().upper() == r2.seq.tostring()
            elif format == "qual" :
                assert isinstance(r2.seq, UnknownSeq)
                assert len(r2) == len(r1)
            else :
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
    if t_format in ["sff", "sff-trim"] :
        #Not allowed empty SFF files.
        continue
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
        if t_format in possible_unknown_seq_formats :
            assert isinstance(record.seq, Seq) or \
                   isinstance(record.seq, UnknownSeq)
        else :
            assert isinstance(record.seq, Seq)
        assert isinstance(record.id, basestring)
        assert isinstance(record.name, basestring)
        assert isinstance(record.description, basestring)
        assert record.id != ""

        if "accessions" in record.annotations :
            accs = record.annotations["accessions"]
            #Check for blanks, or entries with leading/trailing spaces
            for acc in accs :
                assert acc and acc == acc.strip(), \
                    "Bad accession in annotations: %s" % repr(acc)
            assert len(set(accs)) == len(accs), \
                   "Repeated accession in annotations: %s" % repr(accs)
        for ref in record.dbxrefs :
            assert ref and ref == ref.strip(), \
                "Bad cross reference in dbxrefs: %s" % repr(ref)
        assert len(record.dbxrefs) == len(record.dbxrefs), \
               "Repeated cross reference in dbxrefs: %s" % repr(record.dbxrefs)
                
            
        #Check the lists obtained by the different methods agree
        assert compare_record(record, records2[i])
        assert compare_record(record, records3[i])
        assert compare_record(record, records4[i])
        assert compare_record(record, records5[i])

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

    # Check alphabets
    for record in records :
        base_alpha = Alphabet._get_base_alphabet(record.seq.alphabet)
        assert isinstance(base_alpha, Alphabet.SingleLetterAlphabet)
        if t_format in no_alpha_formats :
            assert base_alpha == Alphabet.single_letter_alphabet # Too harsh?
    if base_alpha is None :
        good = []
        bad =[]
        given_alpha=None
    elif isinstance(base_alpha, Alphabet.ProteinAlphabet) :
        good = protein_alphas
        bad = dna_alphas + rna_alphas + nucleotide_alphas
    elif isinstance(base_alpha, Alphabet.RNAAlphabet) :
        good = nucleotide_alphas + rna_alphas
        bad = protein_alphas + dna_alphas
    elif isinstance(base_alpha, Alphabet.DNAAlphabet) :
        good = nucleotide_alphas + dna_alphas
        bad = protein_alphas + rna_alphas
    elif isinstance(base_alpha, Alphabet.NucleotideAlphabet) :
        good = nucleotide_alphas
        bad = protein_alphas
    else :
        assert t_format in no_alpha_formats, "Got %s from %s file" \
               % (repr(base_alpha), t_format)
        good = protein_alphas + dna_alphas + rna_alphas + nucleotide_alphas
        bad = []
    for given_alpha in good :
        #These should all work...
        given_base = Alphabet._get_base_alphabet(given_alpha)
        for record in SeqIO.parse(open(t_filename),t_format,given_alpha):
            base_alpha = Alphabet._get_base_alphabet(record.seq.alphabet)
            assert isinstance(base_alpha, given_base.__class__)
            assert base_alpha == given_base
        if t_count == 1 :
            record = SeqIO.read(open(t_filename),t_format,given_alpha)
            assert isinstance(base_alpha, given_base.__class__)
            assert base_alpha == given_base
    for given_alpha in bad :
        #These should all fail...
        try :
            print SeqIO.parse(open(t_filename),t_format,given_alpha).next()
            assert False, "Forcing wrong alphabet, %s, should fail (%s)" \
                   % (repr(given_alpha), t_filename)
        except ValueError :
            pass
    del good, bad, given_alpha, base_alpha

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
            assert compare_record(records[i], alignment.get_all_seqs()[i])
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
            c = SeqIO.write(records, handle, format)
            assert c == len(records)
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
            if format == "nexus" :
                #The nexus parser will dis-ambiguate repeated record ids.
                assert record.id == new_record.id or \
                       new_record.id.startswith(record.id+".copy")
            else :
                assert record.id == new_record.id
            assert record.seq.tostring() == new_record.seq.tostring()
            #Using compare_record(record, new_record) is too strict

        #Close now, after checking, so that it can be used at the console for debugging
        handle.close()

#Check writers can cope with no alignments
for format in SeqIO._FormatToWriter :
    handle = StringIO()
    try :
        assert 0 == SeqIO.write([], handle, format), \
               "Writing no records to %s format should work!" \
               % t_format
    except ValueError, err:
        print "Writing no records to %s format failed: %s" % (format, err)

print "Finished tested writing files"
