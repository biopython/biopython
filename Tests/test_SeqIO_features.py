# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""SeqFeature related tests for SeqRecord objects from Bio.SeqIO.

Initially this takes matched tests of GenBank and FASTA files from the NCBI
and confirms they are consistent using our different parsers.
"""
import os
import unittest
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, \
                           BeforePosition, AfterPosition, OneOfPosition, \
                           WithinPosition
from StringIO import StringIO
from Bio.SeqIO.InsdcIO import _insdc_feature_location_string

#Top level function as this makes it easier to use for debugging:
def write_read(filename, in_format="gb", out_format="gb") :
    gb_records = list(SeqIO.parse(open(filename),in_format))
    #Write it out...
    handle = StringIO()
    SeqIO.write(gb_records, handle, out_format)
    handle.seek(0)
    #Now load it back and check it agrees,
    gb_records2 = list(SeqIO.parse(handle,out_format))
    compare_records(gb_records, gb_records2)

def compare_record(old, new, ignore_name=False) :
    #Note the name matching is a bit fuzzy
    if not ignore_name \
    and old.id != new.id and old.name != new.name \
    and (old.id not in new.id) and (new.id not in old.id) \
    and (old.id.replace(" ","_") != new.id.replace(" ","_")) :
        raise ValueError("'%s' or '%s' vs '%s' or '%s' records" \
                         % (old.id, old.name, new.id, new.name))
    if len(old.seq) != len(new.seq) :
        raise ValueError("%i vs %i" % (len(old.seq), len(new.seq)))
    if str(old.seq).upper() != str(new.seq).upper() :
        if len(old.seq) < 200 :
            raise ValueError("'%s' vs '%s'" % (old.seq, new.seq))
        else :
            raise ValueError("'%s...' vs '%s...'" % (old.seq[:100], new.seq[:100]))
    if old.features and new.features :
        return compare_features(old.features, new.features)
    #Just insist on at least one word in common:
    if (old.description or new.description) \
    and not set(old.description.split()).intersection(new.description.split()):
        raise ValueError("%s versus %s" \
                         % (repr(old.description), repr(new.description)))
    #TODO - check annotation
    return True

def compare_records(old_list, new_list, ignore_name=False) :
    """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list) :
        raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list) :
        if not compare_record(old,new,ignore_name) :
            return False
    return True

def compare_feature(old, new, ignore_sub_features=False) :
    """Check two SeqFeatures agree."""
    if old.type != new.type :
        raise ValueError("Type %s versus %s" % (old.type, new.type))
    if old.location.nofuzzy_start != new.location.nofuzzy_start \
    or old.location.nofuzzy_end != new.location.nofuzzy_end :
        raise ValueError("%s versus %s:\n%s\nvs:\n%s" \
                         % (old.location, new.location, str(old), str(new)))
    if old.strand != new.strand :
        raise ValueError("Different strand:\n%s\nvs:\n%s" % (str(old), str(new)))
    if old.ref != new.ref :
        raise ValueError("Different ref:\n%s\nvs:\n%s" % (str(old), str(new)))
    if old.ref_db != new.ref_db :
        raise ValueError("Different ref:\n%s\nvs:\n%s" % (str(old), str(new)))
    if old.location.start != new.location.start \
    or str(old.location.start) != str(new.location.start) :
        raise ValueError("Start %s versus %s:\n%s\nvs:\n%s" \
                         % (old.location.start, new.location.start, str(old), str(new)))
    if old.location.end != new.location.end \
    or str(old.location.end) != str(new.location.end) :
        raise ValueError("End %s versus %s:\n%s\nvs:\n%s" \
                         % (old.location.end, new.location.end, str(old), str(new)))
    if not ignore_sub_features :
        if len(old.sub_features) != len(new.sub_features) :
            raise ValueError("Different sub features")
        for a,b in zip(old.sub_features, new.sub_features) :
            if not compare_feature(a,b) :
                return False
    #This only checks key shared qualifiers
    #Would a white list be easier?
    #for key in ["name","gene","translation","codon_table","codon_start","locus_tag"] :
    for key in set(old.qualifiers.keys()).intersection(new.qualifiers.keys()):
        if key in ["db_xref","protein_id","product","note"] :
            #EMBL and GenBank files are use different references/notes/etc
            continue
        if old.qualifiers[key] != new.qualifiers[key] :
            raise ValueError("Qualifier mis-match for %s:\n%s\n%s" \
                             % (key, old.qualifiers[key], new.qualifiers[key]))
    return True

def compare_features(old_list, new_list, ignore_sub_features=False) :
    """Check two lists of SeqFeatures agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list) :
        raise ValueError("%i vs %i features" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list) :
        #This assumes they are in the same order
        if not compare_feature(old,new,ignore_sub_features) :
            return False
    return True

#TODO - Add this functionality to Biopython itself...
def get_feature_nuc(f, parent_seq) :
    if f.sub_features :
        if f.location_operator!="join":
            raise ValueError(f.location_operator)
        #TODO - This should recurse to cope with join(complement(...),...) properly
        #for mixed-strand features, BUT that is impossible with the current GenBank
        #parser due to how the strand is recorded on both the parent and subfeatures.
        f_subs = [parent_seq[f_sub.location.nofuzzy_start:f_sub.location.nofuzzy_end] \
                  for f_sub in f.sub_features]
        #f_subs = [get_feature_nuc(f_sub, parent_seq) for f_sub in f.sub_features]
        #TODO - Join support in Seq object?  But how to deal with alphabets...
        f_seq = Seq("".join(map(str,f_subs)),f_subs[0].alphabet)
    else :
        f_seq = parent_seq[f.location.nofuzzy_start:f.location.nofuzzy_end]
    if f.strand == -1 : f_seq = f_seq.reverse_complement()
    return f_seq

class SeqRecordCreation(unittest.TestCase):
    """Test basic creation of SeqRecords.
    """
    def test_annotations(self):
        """Pass in annotations to SeqRecords.
        """
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test")
        assert rec.annotations == {}
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test",
                        annotations={"test" : ["a test"]})
        assert rec.annotations.get("test", "") == ["a test"]

    def test_letter_annotations(self):
        """Pass in letter annotations to SeqRecords.
        """
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test")
        assert rec.letter_annotations == {}
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test",
                        letter_annotations={"test" : [1, 2, 3, 4]})
        assert rec.letter_annotations.get("test", []) == [1, 2, 3, 4]
        # XXX should raise an error with incorrect number of letter annotations
        rec = SeqRecord(Seq("ACGT", generic_dna),
                        id="Test", name="Test", description="Test",
                        letter_annotations={"test" : [1, 2, 3]})

    def test_qualifiers(self):
        """Pass in qualifiers to SeqFeatures.
        """
        f = SeqFeature(FeatureLocation(10,20), strand=+1, type="CDS")
        assert f.qualifiers == {}
        f = SeqFeature(FeatureLocation(10,20), strand=+1, type="CDS",
                qualifiers={"test": ["a test"]})
        assert f.qualifiers.get("test", "") == ["a test"]

class FeatureWriting(unittest.TestCase) :
    def setUp(self) :
        self.record = SeqRecord(Seq("ACGT"*100, generic_dna),
                                id="Test", name="Test", description="Test")
    def write_read_check(self) :
        handle = StringIO()
        SeqIO.write([self.record], handle, "gb")
        handle.seek(0)
        record2 = SeqIO.read(handle, "gb")
        return compare_record(self.record, record2)

    def test_exact(self) :
        """Features: write/read simple exact locations."""
        #Note we don't have to explicitly give an ExactPosition object,
        #and integer will also work:
        f = SeqFeature(FeatureLocation(10,20), strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "11..20")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(30,40), strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(31..40)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(ExactPosition(50),ExactPosition(60)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "51..60")
        self.record.features.append(f)
        self.write_read_check()

    def test_between(self) :
        """Features: write/read simple between locations."""
        #Note we don't use the BetweenPosition any more!
        f = SeqFeature(FeatureLocation(10,10), strand=+1, type="variation")
        self.assertEqual(_insdc_feature_location_string(f),
                         "10^11")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(20,20), strand=-1, type="variation")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(20^21)")
        self.record.features.append(f)
        self.write_read_check()

    def make_join_feature(self, f_list, ftype="misc_feature"):
        #NOTE - Does NOT reorder the sub-features (which you may
        #want to do for reverse strand features...)
        strands = set(f.strand for f in f_list)
        if len(strands)==1 :
            strand = f_list[0].strand
        else :
            strand = None
        for f in f_list :
            f.type=ftype
        jf = SeqFeature(FeatureLocation(f_list[0].location.start,
                                        f_list[-1].location.end),
                        type=ftype, strand=strand, location_operator="join")
        jf.sub_features = f_list
        return jf
        
    def test_join(self):
        """Features: write/read simple join locations."""
        f1 = SeqFeature(FeatureLocation(10,20), strand=+1)
        f2 = SeqFeature(FeatureLocation(25,40), strand=+1)
        f = self.make_join_feature([f1,f2])
        self.record.features.append(f)
        self.assertEqual(_insdc_feature_location_string(f),
                         "join(11..20,26..40)")
        f1 = SeqFeature(FeatureLocation(110,120), strand=+1)
        f2 = SeqFeature(FeatureLocation(125,140), strand=+1)
        f3 = SeqFeature(FeatureLocation(145,150), strand=+1)
        f = self.make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "join(111..120,126..140,146..150)")
        self.record.features.append(f)
        f1 = SeqFeature(FeatureLocation(210,220), strand=-1)
        f2 = SeqFeature(FeatureLocation(225,240), strand=-1)
        f = self.make_join_feature([f1,f2], ftype="gene")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(join(211..220,226..240))")
        self.record.features.append(f)
        f1 = SeqFeature(FeatureLocation(310,320), strand=-1)
        f2 = SeqFeature(FeatureLocation(325,340), strand=-1)
        f3 = SeqFeature(FeatureLocation(345,350), strand=-1)
        f = self.make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(join(311..320,326..340,346..350))")
        self.record.features.append(f)
        self.write_read_check()

    def test_fuzzy_join(self):
        """Features: write/read fuzzy join locations."""
        f1 = SeqFeature(FeatureLocation(BeforePosition(10),20), strand=+1)
        f2 = SeqFeature(FeatureLocation(25,AfterPosition(40)), strand=+1)
        f = self.make_join_feature([f1,f2])
        self.record.features.append(f)
        self.assertEqual(_insdc_feature_location_string(f),
                         "join(<11..20,26..>40)")
        f1 = SeqFeature(FeatureLocation(OneOfPosition([ExactPosition(107),
                                                       ExactPosition(110)]),120),
                        strand=+1)
        f2 = SeqFeature(FeatureLocation(125,140), strand=+1)
        f3 = SeqFeature(FeatureLocation(145,WithinPosition(150,10)), strand=+1)
        f = self.make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "join(one-of(108,111)..120,126..140,146..(150.160))")
        self.record.features.append(f)
        f1 = SeqFeature(FeatureLocation(BeforePosition(210),220), strand=-1)
        f2 = SeqFeature(FeatureLocation(225,WithinPosition(240,4)), strand=-1)
        f = self.make_join_feature([f1,f2], "gene")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(join(<211..220,226..(240.244)))")
        self.record.features.append(f)
        f1 = SeqFeature(FeatureLocation(AfterPosition(310),320), strand=-1)
        f2 = SeqFeature(FeatureLocation(325,OneOfPosition([ExactPosition(340),
                                                           ExactPosition(337)])),
                        strand=-1)
        f3 = SeqFeature(FeatureLocation(345,WithinPosition(350,5)), strand=-1)
        f = self.make_join_feature([f1,f2,f3], "CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(join(>311..320,326..one-of(340,337),346..(350.355)))")
        self.record.features.append(f)
        self.write_read_check()


    def test_before(self) :
        """Features: write/read simple before locations."""
        f = SeqFeature(FeatureLocation(BeforePosition(5),10), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "<6..10")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(BeforePosition(15),BeforePosition(20)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "<16..<20")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(25,BeforePosition(30)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "26..<30")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(BeforePosition(35),40), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(<36..40)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(BeforePosition(45),BeforePosition(50)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(<46..<50)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(55,BeforePosition(60)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(56..<60)")
        self.record.features.append(f)
        self.write_read_check()
        
    def test_after(self) :
        """Features: write/read simple after locations."""
        f = SeqFeature(FeatureLocation(AfterPosition(5),10), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         ">6..10")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(AfterPosition(15),AfterPosition(20)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         ">16..>20")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(25,AfterPosition(30)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "26..>30")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(AfterPosition(35),40), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(>36..40)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(AfterPosition(45),AfterPosition(50)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(>46..>50)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(55,AfterPosition(60)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(56..>60)")
        self.record.features.append(f)
        self.write_read_check()

    def test_oneof(self) :
        """Features: write/read simple one-of locations."""
        start = OneOfPosition([ExactPosition(0),ExactPosition(3),ExactPosition(6)])
        f = SeqFeature(FeatureLocation(start,21), strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "one-of(1,4,7)..21")
        self.record.features.append(f)
        start = OneOfPosition([ExactPosition(x) for x in [10,13,16]])
        end = OneOfPosition([ExactPosition(x) for x in [41,44,50]])
        f = SeqFeature(FeatureLocation(start,end), strand=+1, type="gene")
        self.assertEqual(_insdc_feature_location_string(f),
                         "one-of(11,14,17)..one-of(41,44,50)")
        self.record.features.append(f)
        end = OneOfPosition([ExactPosition(x) for x in [30,33]])
        f = SeqFeature(FeatureLocation(27,end), strand=+1, type="gene")
        self.assertEqual(_insdc_feature_location_string(f),
                         "28..one-of(30,33)")
        self.record.features.append(f)
        start = OneOfPosition([ExactPosition(x) for x in [36,40]])
        f = SeqFeature(FeatureLocation(start,46), strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(one-of(37,41)..46)")
        self.record.features.append(f)
        start = OneOfPosition([ExactPosition(x) for x in [45,60]])
        end = OneOfPosition([ExactPosition(x) for x in [70,90]])
        f = SeqFeature(FeatureLocation(start,end), strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(one-of(46,61)..one-of(70,90))")
        self.record.features.append(f)
        end = OneOfPosition([ExactPosition(x) for x in [60,63]])
        f = SeqFeature(FeatureLocation(55,end), strand=-1, type="tRNA")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(56..one-of(60,63))")
        self.record.features.append(f)
        self.write_read_check()

    def test_within(self):
        """Features: write/read simple within locations."""
        f = SeqFeature(FeatureLocation(WithinPosition(2,6),10), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "(3.9)..10")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(WithinPosition(12,6),
                                       WithinPosition(20,8)), \
                       strand=+1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "(13.19)..(20.28)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(25,WithinPosition(30,3)), \
                       strand=+1, type="misc_feature")
        self.assertEqual(_insdc_feature_location_string(f),
                         "26..(30.33)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(WithinPosition(35,4),40), \
                       strand=-1, type="rRNA")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement((36.40)..40)")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(WithinPosition(45,2),
                                       WithinPosition(50,3)), \
                       strand=-1, type="repeat_region")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement((46.48)..(50.53))")
        self.record.features.append(f)
        f = SeqFeature(FeatureLocation(55,WithinPosition(60,5)), \
                       strand=-1, type="CDS")
        self.assertEqual(_insdc_feature_location_string(f),
                         "complement(56..(60.65))")
        self.record.features.append(f)
        self.write_read_check()
        
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

    def setUp(self) :
        self.gb_filename = os.path.join("GenBank",self.basename+".gb")
        self.ffn_filename = os.path.join("GenBank",self.basename+".ffn")
        self.faa_filename = os.path.join("GenBank",self.basename+".faa")
        self.fna_filename = os.path.join("GenBank",self.basename+".fna")
        if self.emblname :
            self.embl_filename = os.path.join("EMBL",self.emblname+".embl")

    #These tests only need the GenBank file and the FAA file:
    def test_CDS(self) :
        #"""Checking GenBank CDS translations vs FASTA faa file."""
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        gb_cds = list(SeqIO.parse(open(self.gb_filename),"genbank-cds"))
        fasta = list(SeqIO.parse(open(self.faa_filename),"fasta"))
        compare_records(gb_cds, fasta)
        cds_features = [f for f in gb_record.features if f.type=="CDS"]
        self.assertEqual(len(cds_features), len(fasta))
        for f, r in zip(cds_features, fasta) :
            if r.id in self.skip_trans_test :
                continue
            #Get the nucleotides and translate them
            nuc = get_feature_nuc(f, gb_record.seq)
            pro = nuc.translate(table=self.table, cds=True)
            if pro[-1] == "*" :
                self.assertEqual(str(pro)[:-1], str(r.seq))
            else :
                self.assertEqual(str(pro), str(r.seq))

class NC_005816(NC_000932):
    basename = "NC_005816"
    emblname = "AE017046"
    table = 11
    skip_trans_test = []
    __doc__ = "Tests using %s GenBank and FASTA files from the NCBI" % basename

    def test_GenBank_vs_EMBL(self) :
        if not self.emblname :
            return
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        embl_record = SeqIO.read(open(self.embl_filename),"embl")
        return compare_record(gb_record, embl_record, ignore_name=True)

    def test_Translations(self):
        #"""Checking translation of FASTA features (faa vs ffn)."""
        faa_records = list(SeqIO.parse(open(self.faa_filename),"fasta"))
        ffn_records = list(SeqIO.parse(open(self.ffn_filename),"fasta"))
        self.assertEqual(len(faa_records),len(ffn_records))
        for faa, fna in zip(faa_records, ffn_records) :
            translation = fna.seq.translate(self.table, cds=True)
            if faa.id in self.skip_trans_test :
                continue
            if (str(translation) != str(faa.seq)) \
            and (str(translation) != str(faa.seq)+"*") :
                t = SeqRecord(translation, id="Translation",
                              description="Table %s" % self.table)
                raise ValueError("FAA vs FNA translation problem:\n%s\n%s\n%s\n" \
                                 % (fna.format("fasta"),
                                    t.format("fasta"),
                                    faa.format("fasta")))
    
    def test_Genome(self) :
        #"""Checking GenBank sequence vs FASTA fna file."""
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        fa_record = SeqIO.read(open(self.fna_filename),"fasta")
        compare_record(gb_record, fa_record)
        if self.emblname is None :
            return
        embl_record = SeqIO.read(open(self.embl_filename),"embl")
        compare_record(gb_record, embl_record, ignore_name=True)

    def test_Features(self) :
        #"""Checking GenBank features sequences vs FASTA ffn file."""
        gb_record = SeqIO.read(open(self.gb_filename),"genbank")
        features = [f for f in gb_record.features if f.type=="CDS"]
        fa_records = list(SeqIO.parse(open(self.ffn_filename),"fasta"))
        self.assertEqual(len(fa_records), len(features))
        #This assumes they are in the same order...
        for fa_record, f in zip(fa_records, features) :
            #TODO - check the FASTA ID line against the co-ordinates?
            f_seq = get_feature_nuc(f, gb_record.seq)
            self.assertEqual(len(fa_record.seq),
                             len(f_seq))
            self.assertEqual(str(fa_record.seq),
                             str(f_seq))


class TestWriteRead(unittest.TestCase) :
    """Test can write and read back files."""

    def test_NC_000932(self) :
        """Write and read back NC_000932.gb"""
        write_read(os.path.join("GenBank", "NC_000932.gb"), "gb", "gb")

    def test_NC_005816(self) :
        """Write and read back NC_005816.gb"""
        write_read(os.path.join("GenBank", "NC_005816.gb"), "gb", "gb")

    def test_gbvrl1_start(self) :
        """Write and read back gbvrl1_start.seq"""
        write_read(os.path.join("GenBank", "gbvrl1_start.seq"), "gb", "gb")

    def test_NT_019265(self) :
        """Write and read back NT_019265.gb"""
        write_read(os.path.join("GenBank", "NT_019265.gb"), "gb", "gb")

    def test_cor6(self) :
        """Write and read back cor6_6.gb"""
        write_read(os.path.join("GenBank", "cor6_6.gb"), "gb", "gb")

    def test_arab1(self) :
        """Write and read back arab1.gb"""
        write_read(os.path.join("GenBank", "arab1.gb"), "gb", "gb")

    def test_one_of(self) :
        """Write and read back of_one.gb"""
        write_read(os.path.join("GenBank", "one_of.gb"), "gb", "gb")

    def test_pri1(self) :
        """Write and read back pri1.gb"""
        write_read(os.path.join("GenBank", "pri1.gb"), "gb", "gb")

    def test_noref(self) :
        """Write and read back noref.gb"""
        write_read(os.path.join("GenBank", "noref.gb"), "gb", "gb")

    def test_origin_line(self) :
        """Write and read back origin_line.gb"""
        write_read(os.path.join("GenBank", "origin_line.gb"), "gb", "gb")

    def test_dbsource_wrap(self) :
        """Write and read back dbsource_wrap.gb"""
        write_read(os.path.join("GenBank", "dbsource_wrap.gb"), "gb", "gb")

    def test_blank_seq(self) :
        """Write and read back blank_seq.gb"""
        write_read(os.path.join("GenBank", "blank_seq.gb"), "gb", "gb")

    def test_extra_keywords(self) :
        """Write and read back extra_keywords.gb"""
        write_read(os.path.join("GenBank", "extra_keywords.gb"), "gb", "gb")

    def test_protein_refseq(self) :
        """Write and read back protein_refseq.gb"""
        write_read(os.path.join("GenBank", "protein_refseq.gb"), "gb", "gb")

    def test_protein_refseq2(self) :
        """Write and read back protein_refseq2.gb"""
        write_read(os.path.join("GenBank", "protein_refseq2.gb"), "gb", "gb")

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
