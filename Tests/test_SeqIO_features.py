# Copyright 2009-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SeqFeature related tests for SeqRecord objects from Bio.SeqIO.

Initially this takes matched tests of GenBank and FASTA files from the NCBI
and confirms they are consistent using our different parsers.
"""
import os
import unittest
from Bio.Alphabet import generic_dna, generic_rna, generic_protein
from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, MutableSeq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, \
                           BeforePosition, AfterPosition, OneOfPosition, \
                           WithinPosition
from StringIO import StringIO
from Bio.SeqIO.InsdcIO import _insdc_feature_location_string

#Top level function as this makes it easier to use for debugging:
def write_read(filename, in_format="gb", out_formats=["gb", "embl", "imgt"]):
    for out_format in out_formats:
        gb_records = list(SeqIO.parse(open(filename),in_format))
        #Write it out...
        handle = StringIO()
        SeqIO.write(gb_records, handle, out_format)
        handle.seek(0)
        #Now load it back and check it agrees,
        gb_records2 = list(SeqIO.parse(handle,out_format))
        compare_records(gb_records, gb_records2)

def compare_record(old, new, expect_minor_diffs=False):
    #Note the name matching is a bit fuzzy
    if not expect_minor_diffs \
    and old.id != new.id and old.name != new.name \
    and (old.id not in new.id) and (new.id not in old.id) \
    and (old.id.replace(" ","_") != new.id.replace(" ","_")):
        raise ValueError("'%s' or '%s' vs '%s' or '%s' records" \
                         % (old.id, old.name, new.id, new.name))
    if len(old.seq) != len(new.seq):
        raise ValueError("%i vs %i" % (len(old.seq), len(new.seq)))
    if isinstance(old.seq, UnknownSeq) \
    and isinstance(new.seq, UnknownSeq):
        #Jython didn't like us comparing the string of very long
        #UnknownSeq object (out of heap memory error)
        if old.seq._character.upper() != new.seq._character:
            raise ValueError("%s vs %s" % (repr(old.seq), repr(new.seq)))
    elif str(old.seq).upper() != str(new.seq).upper():
        if len(old.seq) < 200:
            raise ValueError("'%s' vs '%s'" % (old.seq, new.seq))
        else:
            raise ValueError("'%s...' vs '%s...'" % (old.seq[:100], new.seq[:100]))
    if old.features and new.features:
        if not compare_features(old.features, new.features):
            return False
    #Just insist on at least one word in common:
    if (old.description or new.description) \
    and not set(old.description.split()).intersection(new.description.split()):
        raise ValueError("%s versus %s" \
                         % (repr(old.description), repr(new.description)))
    #This only checks common annotation
    #Would a white list be easier?
    for key in set(old.annotations).intersection(new.annotations):
        if key in ["data_file_division", "accessions"]:
            #TODO - These are not yet supported on output, or
            #have other complications (e.g. different number of accessions
            #allowed in various file formats)
            continue
        if key == "comment":
            #Ignore whitespace
            if old.annotations[key].split() != new.annotations[key].split():
                raise ValueError("Annotation mis-match for comment:\n%s\n%s" \
                                % (old.annotations[key], new.annotations[key]))
            continue
        if key == "references":
            if expect_minor_diffs:
                #TODO - Implement EMBL output of references
                continue
            assert len(old.annotations[key]) == len(new.annotations[key])
            for r1, r2 in zip(old.annotations[key], new.annotations[key]):
                assert r1.title == r2.title
                assert r1.authors == r2.authors, \
                       "Old: '%s'\nNew: '%s'" % (r1.authors, r2.authors)
                assert r1.journal == r2.journal
                if r1.consrtm and r2.consrtm:
                    #Not held in EMBL files
                    assert r1.consrtm == r2.consrtm
                if r1.medline_id and r2.medline_id:
                    #Not held in EMBL files
                    assert r1.medline_id == r2.medline_id
                assert r1.pubmed_id == r2.pubmed_id
            continue
        if repr(old.annotations[key]) != repr(new.annotations[key]):
            raise ValueError("Annotation mis-match for %s:\n%s\n%s" \
                             % (key, old.annotations[key], new.annotations[key]))
    return True

def compare_records(old_list, new_list, expect_minor_diffs=False):
    """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list):
        raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list):
        if not compare_record(old,new,expect_minor_diffs):
            return False
    return True

def compare_feature(old, new, ignore_sub_features=False):
    """Check two SeqFeatures agree."""
    if old.type != new.type:
        raise ValueError("Type %s versus %s" % (repr(old.type), repr(new.type)))
    if old.location.nofuzzy_start != new.location.nofuzzy_start \
    or old.location.nofuzzy_end != new.location.nofuzzy_end:
        raise ValueError("%s versus %s:\n%s\nvs:\n%s" \
                         % (old.location, new.location, repr(old), repr(new)))
    if old.strand != new.strand:
        raise ValueError("Different strand:\n%s\nvs:\n%s" % (repr(old), repr(new)))
    if old.ref != new.ref:
        raise ValueError("Different ref:\n%s\nvs:\n%s" % (repr(old), repr(new)))
    if old.ref_db != new.ref_db:
        raise ValueError("Different ref_db:\n%s\nvs:\n%s" % (repr(old), repr(new)))
    if old.location_operator != new.location_operator:
        raise ValueError("Different location_operator:\n%s\nvs:\n%s" % (repr(old), repr(new)))
    if old.location.start != new.location.start \
    or str(old.location.start) != str(new.location.start):
        raise ValueError("Start %s versus %s:\n%s\nvs:\n%s" \
                         % (old.location.start, new.location.start, repr(old), repr(new)))
    if old.location.end != new.location.end \
    or str(old.location.end) != str(new.location.end):
        raise ValueError("End %s versus %s:\n%s\nvs:\n%s" \
                         % (old.location.end, new.location.end, repr(old), repr(new)))
    if not ignore_sub_features:
        if len(old.sub_features) != len(new.sub_features):
            raise ValueError("Different sub features")
        for a,b in zip(old.sub_features, new.sub_features):
            if not compare_feature(a,b):
                return False
    #This only checks key shared qualifiers
    #Would a white list be easier?
    #for key in ["name","gene","translation","codon_table","codon_start","locus_tag"]:
    for key in set(old.qualifiers).intersection(new.qualifiers):
        if key in ["db_xref","protein_id","product","note"]:
            #EMBL and GenBank files are use different references/notes/etc
            continue
        if old.qualifiers[key] != new.qualifiers[key]:
            raise ValueError("Qualifier mis-match for %s:\n%s\n%s" \
                             % (key, old.qualifiers[key], new.qualifiers[key]))
    return True

def compare_features(old_list, new_list, ignore_sub_features=False):
    """Check two lists of SeqFeatures agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list):
        raise ValueError("%i vs %i features" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list):
        #This assumes they are in the same order
        if not compare_feature(old,new,ignore_sub_features):
            return False
    return True

def make_join_feature(f_list, ftype="misc_feature"):
    #NOTE - Does NOT reorder the sub-features (which you may
    #want to do for reverse strand features...)
    if len(set(f.strand for f in f_list))==1:
        strand = f_list[0].strand
    else:
        strand = None
    for f in f_list:
        f.type=ftype
        f.location_operator="join"
    jf = SeqFeature(FeatureLocation(f_list[0].location.start,
                                    f_list[-1].location.end,
                                    strand),
                    type=ftype, location_operator="join")
    assert jf.location.strand == strand
    assert jf.strand == strand
    jf.sub_features = f_list
    return jf

#Prepare a single GenBank record with one feature with a %s place holder for
#the feature location
gbk_template = open("GenBank/iro.gb", "rU").read()
gbk_template = gbk_template.replace('     gene            341..756\n'
                                    '                     /gene="FTCD"\n',
                                    '     misc_feature    %s\n'
                                    '                     /note="Test case"\n')
gbk_template = gbk_template.replace('     exon            341..384\n'
                                    '                     /gene="FTCD"\n'
                                    '                     /number=1\n', '')
gbk_template = gbk_template.replace('     intron          385..617\n'
                                    '                     /gene="FTCD"\n'
                                    '                     /number=1\n', '')
gbk_template = gbk_template.replace('     exon            618..756\n'
                                    '                     /gene="FTCD"\n'
                                    '                     /number=2\n', '')
assert len(gbk_template)==4445
assert gbk_template.count("%") == 1, gbk_template

class SeqFeatureExtractionWritingReading(unittest.TestCase):
    """Tests for SeqFeature sequence extract method, writing, and reading."""

    def check(self, parent_seq, feature, answer_str, location_str):
        self.assertEqual(location_str,
            _insdc_feature_location_string(feature,len(parent_seq)))

        new = feature.extract(parent_seq)
        self.assertTrue(isinstance(new, Seq))
        self.assertEqual(str(new), answer_str)

        if not feature.sub_features:
            new = parent_seq[feature.location.start:feature.location.end]
            if feature.strand == -1:
                new = reverse_complement(new)
            self.assertEqual(str(new), answer_str)

        new = feature.extract(str(parent_seq))
        self.assertTrue(isinstance(new, str))
        self.assertEqual(new, answer_str)

        new = feature.extract(parent_seq.tomutable())
        self.assertTrue(isinstance(new, Seq)) #Not MutableSeq!
        self.assertEqual(str(new), answer_str)

        new = feature.extract(UnknownSeq(len(parent_seq), parent_seq.alphabet))
        self.assertTrue(isinstance(new, UnknownSeq))
        self.assertEqual(len(new), len(answer_str))
        
        if _insdc_feature_location_string(feature, 1326) != location_str:
            #This is to avoid issues with the N^1 between feature which only
            #makes sense at the end of the sequence
            return
        #This template is DNA, but that will still be OK for protein features
        #as they have no strand information... but see below for strand fun
        rec = SeqIO.read(StringIO(gbk_template % location_str), "gb")
        self.assertEqual(1326, len(rec))
        self.assertEqual(2, len(rec.features))
        self.assertEqual(rec.features[0].type, "source")
        self.assertEqual(rec.features[1].type, "misc_feature")
        new_f = rec.features[1]
        self.assertEqual(location_str,
                         _insdc_feature_location_string(new_f,1326))

        #Checking the strand is tricky - on parsing a GenBank file
        #strand +1 is assumed, but our constructed features for the
        #unit test have mostly defaulted to strand None.
        self.assertEqual(len(feature.sub_features), len(new_f.sub_features))
        for f1, f2 in zip(feature.sub_features, new_f.sub_features):
            f1.type = "misc_feature" #hack as may not be misc_feature
            if f1.strand is None:
                f1.strand = f2.strand #hack as described above
            self.assertEqual(f1.strand, f2.strand)
            self.assertTrue(compare_feature(f1,f2))
        feature.type = "misc_feature" #hack as may not be misc_feature
        if not feature.strand:
            feature.strand = new_f.strand #hack as above
        self.assertEqual(feature.strand, new_f.strand)
        self.assertTrue(compare_feature(feature, new_f))
        
        #Some feature method tests
        parent = "ACGT"*250
        s = feature.extract(parent)
        self.assertEqual(len(feature), len(s))
        for i in feature:
            self.assertTrue(i in feature)
        self.assertEqual(set(feature),
                         set(i for i in range(1000) if i in feature))
        if feature.strand == +1:
            self.assertEqual(s, "".join(parent[i] for i in feature))
            

    def test_simple_rna(self):
        """Feature on RNA (simple, default strand)"""
        s = Seq("GAUCRYWSMKHBVDN", generic_rna)
        f = SeqFeature(FeatureLocation(5,10))
        self.assertEqual(f.strand, None)
        self.assertEqual(f.location.strand, None)
        self.check(s, f, "YWSMK", "6..10")

    def test_simple_dna(self):
        """Feature on DNA (simple, default strand)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(5,10))
        self.check(s, f, "YWSMK", "6..10")

    def test_single_letter_dna(self):
        """Feature on DNA (single letter, default strand)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(5,6))
        self.check(s, f, "Y", "6")

    def test_zero_len_dna(self):
        """Feature on DNA (between location, zero length, default strand)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(5,5))
        self.check(s, f, "", "5^6")

    def test_zero_len_dna_end(self):
        """Feature on DNA (between location at end, zero length, default strand)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(15,15))
        self.check(s, f, "", "15^1")

    def test_simple_dna_strand0(self):
        """Feature on DNA (simple, strand 0)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(5,10), strand=0)
        self.check(s, f, "YWSMK", "6..10")

    def test_simple_dna_strand_none(self):
        """Feature on DNA (simple, strand None)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(5,10), strand=None)
        self.check(s, f, "YWSMK", "6..10")

    def test_simple_dna_strand1(self):
        """Feature on DNA (simple, strand +1)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(5,10), strand=1)
        self.assertEqual(f.strand, +1)
        self.assertEqual(f.location.strand, +1)
        self.check(s, f, "YWSMK", "6..10")
        
    def test_simple_dna_strand_minus(self):
        """Feature on DNA (simple, strand -1)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f = SeqFeature(FeatureLocation(5,10), strand=-1)
        self.assertEqual(f.strand, -1)
        self.assertEqual(f.location.strand, -1)
        self.check(s, f, "MKSWR", "complement(6..10)")

    def test_simple_dna_join(self):
        """Feature on DNA (join, strand +1)"""
        s = Seq("GATCRYWSMKHBVDN", generic_dna)
        f1 = SeqFeature(FeatureLocation(5,10), strand=1)
        f2 = SeqFeature(FeatureLocation(12,15), strand=1)
        f = make_join_feature([f1,f2])
        self.check(s, f, "YWSMKVDN", "join(6..10,13..15)")

    def test_simple_dna_join(self):
        """Feature on DNA (join, strand -1)"""
        s = Seq("AAAAACCCCCTTTTTGGGGG", generic_dna)
        f1 = SeqFeature(FeatureLocation(5,10), strand=-1)
        f2 = SeqFeature(FeatureLocation(12,15), strand=-1)
        f = make_join_feature([f1,f2])
        self.check(s, f, reverse_complement("CCCCC"+"TTT"),
                   "complement(join(6..10,13..15))")

    def test_simple_dna_join(self):
        """Feature on DNA (join, strand -1, before position)"""
        s = Seq("AAAAACCCCCTTTTTGGGGG", generic_dna)
        f1 = SeqFeature(FeatureLocation(BeforePosition(5),10), strand=-1)
        f2 = SeqFeature(FeatureLocation(12,15), strand=-1)
        f = make_join_feature([f1,f2])
        self.check(s, f, reverse_complement("CCCCC"+"TTT"),
                   "complement(join(<6..10,13..15))")

    def test_simple_dna_join_after(self):
        """Feature on DNA (join, strand -1, after position)"""
        s = Seq("AAAAACCCCCTTTTTGGGGG", generic_dna)
        f1 = SeqFeature(FeatureLocation(5,10), strand=-1)
        f2 = SeqFeature(FeatureLocation(12,AfterPosition(15)), strand=-1)
        f = make_join_feature([f1,f2])
        self.check(s, f, reverse_complement("CCCCC"+"TTT"),
                   "complement(join(6..10,13..>15))")

    def test_mixed_strand_dna_join(self):
        """Feature on DNA (join, mixed strand)"""
        s = Seq("AAAAACCCCCTTTTTGGGGG", generic_dna)
        f1 = SeqFeature(FeatureLocation(5,10), strand=+1)
        f2 = SeqFeature(FeatureLocation(12,15), strand=-1)
        f = make_join_feature([f1,f2])
        self.check(s, f, "CCCCC"+reverse_complement("TTT"),
                   "join(6..10,complement(13..15))")

    def test_mixed_strand_dna_multi_join(self):
        """Feature on DNA (multi-join, mixed strand)"""
        s = Seq("AAAAACCCCCTTTTTGGGGG", generic_dna)
        f1 = SeqFeature(FeatureLocation(5,10), strand=+1)
        f2 = SeqFeature(FeatureLocation(12,15), strand=-1)
        f3 = SeqFeature(FeatureLocation(BeforePosition(0),5), strand=+1)
        f = make_join_feature([f1,f2,f3])
        self.check(s, f, "CCCCC"+reverse_complement("TTT")+"AAAAA",
                   "join(6..10,complement(13..15),<1..5)")

    def test_protein_simple(self):
        """Feature on protein (simple)"""
        s = Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ", generic_protein)
        f = SeqFeature(FeatureLocation(5,10))
        self.check(s, f, "FGHIJ", "6..10")

    def test_protein_join(self):
        """Feature on protein (join)"""
        s = Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ", generic_protein)
        f1 = SeqFeature(FeatureLocation(5,10))
        f2 = SeqFeature(FeatureLocation(15,20))
        f = make_join_feature([f1,f2])
        self.check(s, f, "FGHIJ"+"PQRST", "join(6..10,16..20)")

    def test_protein_join_fuzzy(self):
        """Feature on protein (fuzzy join)"""
        s = Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ", generic_protein)
        f1 = SeqFeature(FeatureLocation(BeforePosition(5),10))
        f2 = SeqFeature(FeatureLocation(OneOfPosition(15, (ExactPosition(15),
                                                           ExactPosition(16))),
                                        AfterPosition(20)))
        f = make_join_feature([f1,f2])
        self.check(s, f, "FGHIJ"+"PQRST", "join(<6..10,one-of(16,17)..>20)")

    def test_protein_multi_join(self):
        """Feature on protein (multi-join)"""
        s = Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ", generic_protein)
        f1 = SeqFeature(FeatureLocation(1,2))
        f2 = SeqFeature(FeatureLocation(8,9))
        f3 = SeqFeature(FeatureLocation(14,16))
        f4 = SeqFeature(FeatureLocation(24,25))
        f5 = SeqFeature(FeatureLocation(19,20))
        f6 = SeqFeature(FeatureLocation(7,8))
        f7 = SeqFeature(FeatureLocation(14,15))
        f8 = SeqFeature(FeatureLocation(13,14))
        f = make_join_feature([f1,f2,f3,f4,f5,f6,f7,f8])
        self.check(s, f, "BIOPYTHON", "join(2,9,15..16,25,20,8,15,14)")

    def test_protein_between(self):
        """Feature on protein (between location, zero length)"""
        s = Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ", generic_protein)
        f = SeqFeature(FeatureLocation(5,5))
        self.check(s, f, "", "5^6")

    def test_protein_oneof(self):
        """Feature on protein (one-of positions)"""
        s = Seq("ABCDEFGHIJKLMNOPQRSTUVWXYZ", generic_protein)
        start = OneOfPosition(5, (ExactPosition(5),ExactPosition(7)))
        end = OneOfPosition(11, (ExactPosition(10),ExactPosition(11)))
        f = SeqFeature(FeatureLocation(start,10))
        self.check(s, f, "FGHIJ", "one-of(6,8)..10")
        f = SeqFeature(FeatureLocation(start,end))
        self.check(s, f, "FGHIJK", "one-of(6,8)..one-of(10,11)")
        f = SeqFeature(FeatureLocation(5,end))
        self.check(s, f, "FGHIJK", "6..one-of(10,11)")


class SeqFeatureCreation(unittest.TestCase):
    """Test basic creation of SeqFeatures.
    """

    def test_qualifiers(self):
        """Pass in qualifiers to SeqFeatures.
        """
        f = SeqFeature(FeatureLocation(10,20), strand=+1, type="CDS")
        self.assertEqual(f.qualifiers, {})
        f = SeqFeature(FeatureLocation(10,20), strand=+1, type="CDS",
                qualifiers={"test": ["a test"]})
        self.assertEqual(f.qualifiers["test"], ["a test"])

class FeatureWriting(unittest.TestCase):
    def setUp(self):
        self.record = SeqRecord(Seq("ACGT"*100, generic_dna),
                                id="Test", name="Test", description="Test")
    def write_read_check(self, format):
        handle = StringIO()
        SeqIO.write([self.record], handle, format)
        handle.seek(0)
        record2 = SeqIO.read(handle, format)
        compare_record(self.record, record2)
    
    def write_read_checks(self, formats=["gb", "embl", "imgt"]):
        for f in formats:
            self.write_read_check(f)

    def test_exact(self):
        """GenBank/EMBL write/read simple exact locations."""
        #Note we don't have to explicitly give an ExactPosition object,
        #an integer will also work:
        f = SeqFeature(FeatureLocation(10,20), strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "11..20")
        self.assertEqual(_insdc_feature_location_string(f._flip(20),20),
                         "complement(1..10)")
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "complement(81..90)")
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(30,40), strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(31..40)")
        self.assertEqual(_insdc_feature_location_string(f._flip(40),40),
                         "1..10")
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "61..70")
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)

        f = SeqFeature(FeatureLocation(ExactPosition(50),ExactPosition(60)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "51..60")
        self.assertEqual(_insdc_feature_location_string(f._flip(60),60),
                         "complement(1..10)")
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "complement(41..50)")
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)

        self.write_read_checks()
        #The write/check won't work on strandless features due to the
        #limitations of the GenBank (and EMBL) feature location scheme
        for s in [0, None] :
            #Check flipping of a simple strand 0 feature:
            f = SeqFeature(FeatureLocation(0,100), strand=s, type="source")
            self.assertEqual(_insdc_feature_location_string(f,100),
                             "1..100")
            self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                             "1..100")
            self.assertEqual(_insdc_feature_location_string(f._flip(200),200),
                             "101..200")
            self.assertEqual(f._flip(100).strand, f.strand)

    def test_between(self):
        """GenBank/EMBL write/read simple between locations."""
        #Note we don't use the BetweenPosition any more!
        f = SeqFeature(FeatureLocation(10,10), strand=+1, type="variation")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "10^11")
        self.assertEqual(_insdc_feature_location_string(f._flip(20),20),
                         "complement(10^11)")
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "complement(90^91)")
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(20,20), strand=-1, type="variation")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(20^21)")
        self.assertEqual(_insdc_feature_location_string(f._flip(40),40),
                         "20^21")
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "80^81")
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        self.write_read_checks()

    def test_join(self):
        """GenBank/EMBL write/read simple join locations."""
        f1 = SeqFeature(FeatureLocation(10,20), strand=+1)
        f2 = SeqFeature(FeatureLocation(25,40), strand=+1)
        f = make_join_feature([f1,f2])
        self.record.features.append(f)
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "join(11..20,26..40)")
        self.assertEqual(_insdc_feature_location_string(f._flip(60),60),
                         "complement(join(21..35,41..50))")
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "complement(join(61..75,81..90))")
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand, -1)
        self.assertEqual(f._flip(100).strand, -1)
        f1 = SeqFeature(FeatureLocation(110,120), strand=+1)
        f2 = SeqFeature(FeatureLocation(125,140), strand=+1)
        f3 = SeqFeature(FeatureLocation(145,150), strand=+1)
        f = make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "join(111..120,126..140,146..150)")
        self.assertEqual(_insdc_feature_location_string(f._flip(150),150),
                         "complement(join(1..5,11..25,31..40))")
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand,-1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        f1 = SeqFeature(FeatureLocation(210,220), strand=-1)
        f2 = SeqFeature(FeatureLocation(225,240), strand=-1)
        f = make_join_feature([f1,f2], ftype="gene")
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "complement(join(211..220,226..240))")
        self.assertEqual(_insdc_feature_location_string(f._flip(300),300),
                         "join(61..75,81..90)")
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand, +1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        f1 = SeqFeature(FeatureLocation(310,320), strand=-1)
        f2 = SeqFeature(FeatureLocation(325,340), strand=-1)
        f3 = SeqFeature(FeatureLocation(345,350), strand=-1)
        f = make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "complement(join(311..320,326..340,346..350))")
        self.assertEqual(_insdc_feature_location_string(f._flip(350),350),
                         "join(1..5,11..25,31..40)")
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand, +1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        self.write_read_checks()

    def test_fuzzy_join(self):
        """Features: write/read fuzzy join locations."""
        s = "N"*500
        f1 = SeqFeature(FeatureLocation(BeforePosition(10),20), strand=+1)
        f2 = SeqFeature(FeatureLocation(25,AfterPosition(40)), strand=+1)
        f = make_join_feature([f1,f2])
        self.record.features.append(f)
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "join(<11..20,26..>40)")
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "complement(join(<61..75,81..>90))")
        self.assertEqual(f.strand, +1)
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand, -1)
        self.assertEqual(f._flip(100).strand, -1)
        
        f1 = SeqFeature(FeatureLocation(OneOfPosition(107, [ExactPosition(107),
                                                            ExactPosition(110)]),
                                        120),
                        strand=+1)
        f2 = SeqFeature(FeatureLocation(125,140), strand=+1)
        f3 = SeqFeature(FeatureLocation(145,WithinPosition(160, left=150, right=160)), strand=+1)
        f = make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "join(one-of(108,111)..120,126..140,146..(150.160))")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(200),200),
                         "complement(join((41.51)..55,61..75,81..one-of(90,93)))")
        self.assertEqual(f.strand, +1)
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand,-1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f1 = SeqFeature(FeatureLocation(BeforePosition(210),220), strand=-1)
        f2 = SeqFeature(FeatureLocation(225,WithinPosition(244, left=240, right=244)), strand=-1)
        f = make_join_feature([f1,f2], "gene")
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "complement(join(<211..220,226..(240.244)))")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(300),300),
                         "join((57.61)..75,81..>90)")
        self.assertEqual(f.strand, -1)
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand, +1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        f1 = SeqFeature(FeatureLocation(AfterPosition(310),320), strand=-1)
        #Note - is one-of(340,337) allowed or should it be one-of(337,340)?
        f2 = SeqFeature(FeatureLocation(325,OneOfPosition(340, [ExactPosition(340),
                                                                ExactPosition(337)])),
                        strand=-1)
        f3 = SeqFeature(FeatureLocation(345,WithinPosition(355, left=350, right=355)), strand=-1)
        f = make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f,500),
                         "complement(join(>311..320,326..one-of(340,337),346..(350.355)))")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(400),400),
                         "join((46.51)..55,one-of(64,61)..75,81..<90)")
        self.assertEqual(f.strand, -1)
        for sub_f in f._flip(100).sub_features :
            self.assertEqual(sub_f.strand, +1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        self.write_read_checks()

    def test_before(self):
        """Features: write/read simple before locations."""
        s = "N"*200
        f = SeqFeature(FeatureLocation(BeforePosition(5),10), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "<6..10")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(20),20),
                         "complement(11..>15)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(BeforePosition(15),BeforePosition(20)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "<16..<20")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(20),20),
                         "complement(>1..>5)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(25,BeforePosition(30)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "26..<30")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(40),40),
                         "complement(>11..15)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(BeforePosition(35),40), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(<36..40)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(40),40),
                         "1..>5")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(BeforePosition(45),BeforePosition(50)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(<46..<50)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         ">51..>55")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(55,BeforePosition(60)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(56..<60)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         ">41..45")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        self.write_read_checks()
        
    def test_after(self):
        """Features: write/read simple after locations."""
        s = "N" * 200
        f = SeqFeature(FeatureLocation(AfterPosition(5),10), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         ">6..10")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "complement(91..<95)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)

        f = SeqFeature(FeatureLocation(AfterPosition(15),AfterPosition(20)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         ">16..>20")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(20),20),
                         "complement(<1..<5)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)

        f = SeqFeature(FeatureLocation(25,AfterPosition(30)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "26..>30")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(30),30),
                         "complement(<1..5)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)

        f = SeqFeature(FeatureLocation(AfterPosition(35),40), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(>36..40)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "61..<65")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)

        f = SeqFeature(FeatureLocation(AfterPosition(45),AfterPosition(50)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(>46..>50)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "<51..<55")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)

        f = SeqFeature(FeatureLocation(55,AfterPosition(60)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(56..>60)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "<41..45")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)

        self.write_read_checks()

    def test_oneof(self):
        """Features: write/read simple one-of locations."""
        s = "N" * 100
        start = OneOfPosition(0, [ExactPosition(0),ExactPosition(3),ExactPosition(6)])
        f = SeqFeature(FeatureLocation(start,21), strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "one-of(1,4,7)..21")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "complement(80..one-of(94,97,100))")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        start = OneOfPosition(10, [ExactPosition(x) for x in [10,13,16]])
        end = OneOfPosition(50, [ExactPosition(x) for x in [41,44,50]])
        f = SeqFeature(FeatureLocation(start,end), strand=+1, type="gene")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "one-of(11,14,17)..one-of(41,44,50)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(50),50),
                         "complement(one-of(1,7,10)..one-of(34,37,40))")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        end = OneOfPosition(33, [ExactPosition(x) for x in [30,33]])
        f = SeqFeature(FeatureLocation(27,end), strand=+1, type="gene")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "28..one-of(30,33)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(40),40),
                         "complement(one-of(8,11)..13)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        start = OneOfPosition(36, [ExactPosition(x) for x in [36,40]])
        f = SeqFeature(FeatureLocation(start,46), strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(one-of(37,41)..46)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(50),50),
                         "5..one-of(10,14)")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        start = OneOfPosition(45, [ExactPosition(x) for x in [45,60]])
        end = OneOfPosition(90, [ExactPosition(x) for x in [70,90]])
        f = SeqFeature(FeatureLocation(start,end), strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(one-of(46,61)..one-of(70,90))")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "one-of(11,31)..one-of(40,55)")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        end = OneOfPosition(63, [ExactPosition(x) for x in [60,63]])
        f = SeqFeature(FeatureLocation(55,end), strand=-1, type="tRNA")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(56..one-of(60,63))")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "one-of(38,41)..45")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        self.write_read_checks()

    def test_within(self):
        """Features: write/read simple within locations."""
        s = "N" * 100
        f = SeqFeature(FeatureLocation(WithinPosition(2, left=2, right=8),10), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "(3.9)..10")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(20),20),
                         "complement(11..(12.18))")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(WithinPosition(12, left=12, right=18),
                                       WithinPosition(28, left=20, right=28)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "(13.19)..(20.28)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(30),30),
                         "complement((3.11)..(12.18))")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(25,WithinPosition(33, left=30, right=33)), \
                       strand=+1, type="misc_feature")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "26..(30.33)")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(40),40),
                         "complement((8.11)..15)")
        self.assertEqual(f.strand, +1)
        self.assertEqual(f._flip(100).strand, -1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(WithinPosition(35, left=35, right=39),40), \
                       strand=-1, type="rRNA")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement((36.40)..40)")
        self.assertEqual(_insdc_feature_location_string(f._flip(40),40),
                         "1..(1.5)")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(WithinPosition(45, left=45, right=47),
                                       WithinPosition(53, left=50, right=53)), \
                       strand=-1, type="repeat_region")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement((46.48)..(50.53))")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(60),60),
                         "(8.11)..(13.15)")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        f = SeqFeature(FeatureLocation(55,WithinPosition(65, left=60, right=65)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f,100),
                         "complement(56..(60.65))")
        self.assertEqual(len(f), len(f.extract(s)))
        self.assertEqual(_insdc_feature_location_string(f._flip(100),100),
                         "(36.41)..45")
        self.assertEqual(f.strand, -1)
        self.assertEqual(f._flip(100).strand, +1)
        self.record.features.append(f)
        
        self.write_read_checks()
        
class NC_000932(unittest.TestCase):
    #This includes an evil dual strand gene
    basename = "NC_000932"
    emblname = None # "AP000423" has different annotation (e.g. more CDS)
    table = 11
    skip_trans_test = ["gi|7525080|ref|NP_051037.1|", #dual-strand
                       "gi|7525057|ref|NP_051038.1|", #dual-strand
                       "gi|90110725|ref|NP_051109.2|", #Invalid annotation? No start codon
                       ]
    __doc__ = "Tests using %s GenBank and FASTA files from the NCBI" % basename
    #TODO - neat way to change the docstrings...

    def setUp(self):
        self.gb_filename = os.path.join("GenBank",self.basename+".gb")
        self.ffn_filename = os.path.join("GenBank",self.basename+".ffn")
        self.faa_filename = os.path.join("GenBank",self.basename+".faa")
        self.fna_filename = os.path.join("GenBank",self.basename+".fna")
        if self.emblname:
            self.embl_filename = os.path.join("EMBL",self.emblname+".embl")

    #These tests only need the GenBank file and the FAA file:
    def test_CDS(self):
        #"""Checking GenBank CDS translations vs FASTA faa file."""
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        gb_cds = list(SeqIO.parse(open(self.gb_filename),"genbank-cds"))
        fasta = list(SeqIO.parse(open(self.faa_filename),"fasta"))
        compare_records(gb_cds, fasta)
        cds_features = [f for f in gb_record.features if f.type=="CDS"]
        self.assertEqual(len(cds_features), len(fasta))
        for f, r in zip(cds_features, fasta):
            if r.id in self.skip_trans_test:
                continue
            #Get the nucleotides and translate them
            nuc = f.extract(gb_record.seq)
            self.assertEqual(len(nuc), len(f))
            pro = nuc.translate(table=self.table, cds=True)
            #print r.id, nuc, pro, r.seq
            #print f
            if pro[-1] == "*":
                self.assertEqual(str(pro)[:-1], str(r.seq))
            else:
                self.assertEqual(str(pro), str(r.seq))

class NC_005816(NC_000932):
    basename = "NC_005816"
    emblname = "AE017046"
    table = 11
    skip_trans_test = []
    __doc__ = "Tests using %s GenBank and FASTA files from the NCBI" % basename

    def test_GenBank_vs_EMBL(self):
        if not self.emblname:
            return
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        embl_record = SeqIO.read(open(self.embl_filename),"embl")
        if len(embl_record.features) < len(gb_record.features):
            #Used to match, but I've update the GenBank files
            #which now has lots of misc_feature entries not in EMBL
            embl_record.features = [f for f in embl_record.features \
                                    if f.type != "misc_feature"]
            gb_record.features = [f for f in gb_record.features \
                                  if f.type != "misc_feature"]
        return compare_record(gb_record, embl_record, expect_minor_diffs=True)

    def test_Translations(self):
        #"""Checking translation of FASTA features (faa vs ffn)."""
        faa_records = list(SeqIO.parse(open(self.faa_filename),"fasta"))
        ffn_records = list(SeqIO.parse(open(self.ffn_filename),"fasta"))
        self.assertEqual(len(faa_records),len(ffn_records))
        for faa, fna in zip(faa_records, ffn_records):
            translation = fna.seq.translate(self.table, cds=True)
            if faa.id in self.skip_trans_test:
                continue
            if (str(translation) != str(faa.seq)) \
            and (str(translation) != str(faa.seq)+"*"):
                t = SeqRecord(translation, id="Translation",
                              description="Table %s" % self.table)
                raise ValueError("FAA vs FNA translation problem:\n%s\n%s\n%s\n" \
                                 % (fna.format("fasta"),
                                    t.format("fasta"),
                                    faa.format("fasta")))
    
    def test_Genome(self):
        #"""Checking GenBank sequence vs FASTA fna file."""
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        fa_record = SeqIO.read(open(self.fna_filename),"fasta")
        compare_record(gb_record, fa_record)
        if self.emblname is None:
            return
        embl_record = SeqIO.read(open(self.embl_filename),"embl")
        if len(embl_record.features) < len(gb_record.features):
            #Hack since now out of sync for NC_005816
            embl_record.features = [f for f in embl_record.features \
                                    if f.type != "misc_feature"]
            gb_record.features = [f for f in gb_record.features \
                                  if f.type != "misc_feature"]
        compare_record(gb_record, embl_record, expect_minor_diffs=True)

    def test_Features(self):
        #"""Checking GenBank features sequences vs FASTA ffn file."""
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        features = [f for f in gb_record.features if f.type=="CDS"]
        fa_records = list(SeqIO.parse(open(self.ffn_filename),"fasta"))
        self.assertEqual(len(fa_records), len(features))
        #This assumes they are in the same order...
        for fa_record, f in zip(fa_records, features):
            #TODO - check the FASTA ID line against the co-ordinates?
            f_seq = f.extract(gb_record.seq)
            self.assertEqual(len(fa_record.seq),
                             len(f_seq))
            self.assertEqual(str(fa_record.seq),
                             str(f_seq))
            self.assertEqual(len(f_seq), len(f))


class TestWriteRead(unittest.TestCase):
    """Test can write and read back files."""
    def test_NC_000932(self):
        """Write and read back NC_000932.gb"""
        write_read(os.path.join("GenBank", "NC_000932.gb"), "gb")

    def test_NC_005816(self):
        """Write and read back NC_005816.gb"""
        write_read(os.path.join("GenBank", "NC_005816.gb"), "gb")

    def test_gbvrl1_start(self):
        """Write and read back gbvrl1_start.seq"""
        write_read(os.path.join("GenBank", "gbvrl1_start.seq"), "gb")

    def test_NT_019265(self):
        """Write and read back NT_019265.gb"""
        write_read(os.path.join("GenBank", "NT_019265.gb"), "gb")

    def test_cor6(self):
        """Write and read back cor6_6.gb"""
        write_read(os.path.join("GenBank", "cor6_6.gb"), "gb")

    def test_arab1(self):
        """Write and read back arab1.gb"""
        write_read(os.path.join("GenBank", "arab1.gb"), "gb")

    def test_one_of(self):
        """Write and read back of_one.gb"""
        write_read(os.path.join("GenBank", "one_of.gb"), "gb")

    def test_pri1(self):
        """Write and read back pri1.gb"""
        write_read(os.path.join("GenBank", "pri1.gb"), "gb")

    def test_noref(self):
        """Write and read back noref.gb"""
        write_read(os.path.join("GenBank", "noref.gb"), "gb")

    def test_origin_line(self):
        """Write and read back origin_line.gb"""
        write_read(os.path.join("GenBank", "origin_line.gb"), "gb")

    def test_dbsource_wrap(self):
        """Write and read back dbsource_wrap.gb"""
        write_read(os.path.join("GenBank", "dbsource_wrap.gb"), "gb", ["gb"])
        #Protein so can't convert this to EMBL format

    def test_blank_seq(self):
        """Write and read back blank_seq.gb"""
        write_read(os.path.join("GenBank", "blank_seq.gb"), "gb", ["gb"])
        #Protein so can't convert this to EMBL format

    def test_extra_keywords(self):
        """Write and read back extra_keywords.gb"""
        write_read(os.path.join("GenBank", "extra_keywords.gb"), "gb")

    def test_protein_refseq(self):
        """Write and read back protein_refseq.gb"""
        write_read(os.path.join("GenBank", "protein_refseq.gb"), "gb", ["gb"])
        #Protein so can't convert this to EMBL format

    def test_protein_refseq2(self):
        """Write and read back protein_refseq2.gb"""
        write_read(os.path.join("GenBank", "protein_refseq2.gb"), "gb", ["gb"])
        #Protein so can't convert this to EMBL format

    def test_AAA03323(self):
        """Write and read back AAA03323.embl"""
        write_read(os.path.join("EMBL", "AAA03323.embl"), "embl")

    def test_AE017046(self):
        """Write and read back AE017046.embl"""
        write_read(os.path.join("EMBL", "AE017046.embl"), "embl")

    def test_DD231055_edited(self):
        """Write and read back DD231055_edited.embl"""
        write_read(os.path.join("EMBL", "DD231055_edited.embl"), "embl")

    def test_Human_contigs(self):
        """Write and read back Human_contigs.embl"""
        write_read(os.path.join("EMBL", "Human_contigs.embl"), "embl")

    def test_SC10H5(self):
        """Write and read back SC10H5.embl"""
        write_read(os.path.join("EMBL", "SC10H5.embl"), "embl")

    def test_TRBG361(self):
        """Write and read back TRBG361.embl"""
        write_read(os.path.join("EMBL", "TRBG361.embl"), "embl")

    def test_U87107(self):
        """Write and read back U87107.embl"""
        write_read(os.path.join("EMBL", "U87107.embl"), "embl")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
