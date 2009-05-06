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
from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, MutableSeq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO

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
    if old.location.start != new.location.start :
        raise ValueError("Start %s versus %s:\n%s\nvs:\n%s" \
                         % (old.location.start, new.location.start, str(old), str(new)))
    if old.location.end != new.location.end :
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

from Bio.Data import CodonTable
def methionine_translate(nuc, table) :
    """Hack until we fix Bug 2783."""
    translation = nuc.translate(table)
    start_codons = CodonTable.ambiguous_dna_by_id[table].start_codons
    #There may not be a start codon, for example:
    #Consider NC_006980, protein ID XP_627884.1, <58180..59604, RVSSSSLLF...
    if str(nuc)[:3] in start_codons and translation[0] != "M" :
        translation = Seq("M", translation.alphabet) + translation[1:]
    return translation

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

        
class NC_000932(unittest.TestCase):
    #This includes an evil dual strand gene (trans-splicing!)
    basename = "NC_000932"
    emblname = None # "AP000423" has different annotation (e.g. more CDS)
    table = 11
    skip_trans_test = ["gi|7525080|ref|NP_051037.1|", #Trans-splicing
                       "gi|7525057|ref|NP_051038.1|", #Trans-splicing
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
            pro = methionine_translate(nuc, self.table)
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
            translation = methionine_translate(fna.seq, self.table)
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
