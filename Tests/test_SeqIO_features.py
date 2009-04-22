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

#Top level function as this makes it easier to use for debugging:
def compare_record(old, new) :
    #Note the name matching is a bit fuzzy, e.g. truncation and
    #no spaces in PHYLIP files.
    if old.id != new.id and old.name != new.name \
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
    if old.features and new.features \
    and len(old.features) != len(new.features) :
        raise ValueError("%i vs %i features" \
                         % (len(old.features, len(new.features))))
    #TODO - check annotation

def compare_records(old_list, new_list) :
    """Check two lists of SeqRecords agree, raises a ValueError if mismatch."""
    if len(old_list) != len(new_list) :
        raise ValueError("%i vs %i records" % (len(old_list), len(new_list)))
    for old, new in zip(old_list, new_list) :
        if not compare_record(old,new) :
            return False
    return True

class GenBankTests(unittest.TestCase):
    basename = "NC_005816"
    table = 11
    skip_trans_test = []
    
    def test_CDS(self) :
        """Checking GenBank CDS translations vs FASTA faa file."""
        filename = os.path.join("GenBank",self.basename+".gb")
        gb_cds = list(SeqIO.parse(open(filename),"genbank-cds"))
        filename = os.path.join("GenBank",self.basename+".faa")
        fasta = list(SeqIO.parse(open(filename),"fasta"))
        compare_records(gb_cds, fasta)

    def test_Translations(self):
        """Checking translation of FASTA features (faa vs ffn)."""
        filename = os.path.join("GenBank",self.basename+".faa")
        faa_records = list(SeqIO.parse(open(filename),"fasta"))
        filename = os.path.join("GenBank",self.basename+".ffn")
        ffn_records = list(SeqIO.parse(open(filename),"fasta"))
        self.assertEqual(len(faa_records),len(ffn_records))
        for faa, fna in zip(faa_records, ffn_records) :
            translation = fna.seq.translate(self.table)
            if translation[0] != "M" and faa.seq[0] == "M" :
                #Hack until we fix Bug 2783
                translation = Seq("M", translation.alphabet) + translation[1:]
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
        """Checking GenBank sequence vs FASTA fna file."""
        filename = os.path.join("GenBank",self.basename+".gb")
        gb_record = SeqIO.read(open(filename),"genbank")
        filename = os.path.join("GenBank",self.basename+".fna")
        fa_record = SeqIO.read(open(filename),"fasta")
        compare_record(gb_record, fa_record)

    def test_Features(self) :
        """Checking GenBank features sequences vs FASTA ffn file."""
        filename = os.path.join("GenBank",self.basename+".gb")
        gb_record = SeqIO.read(open(filename),"genbank")
        features = [f for f in gb_record.features if f.type=="CDS"]
        filename = os.path.join("GenBank",self.basename+".ffn")
        fa_records = list(SeqIO.parse(open(filename),"fasta"))
        self.assertEqual(len(fa_records), len(features))
        #This assumes they are in the same order...
        for fa_record, f in zip(fa_records, features) :
            #TODO - check the FASTA ID line against the co-ordinates?
            if f.sub_features :
                #TODO - Add this functionality to Biopython itself...
                self.assertEqual(f.location_operator,"join")
                f_subs = [gb_record.seq[f_sub.location.nofuzzy_start:f_sub.location.nofuzzy_end] \
                          for f_sub in f.sub_features]
                f_seq = Seq("".join(map(str,f_subs)),f_subs[0].alphabet)
            else :
                f_seq = gb_record.seq[f.location.nofuzzy_start:f.location.nofuzzy_end]
            if f.strand == -1 : f_seq = f_seq.reverse_complement()
            self.assertEqual(len(fa_record.seq),
                             len(f_seq))
            self.assertEqual(str(fa_record.seq),
                             str(f_seq))

class GenBankTests2(GenBankTests):
    #This includes several joins and fuzzy joins :)
    basename = "NC_006980"
    table = 1
    skip_trans_test = ["gi|126643907|ref|XP_001388140.1|"]
    #TODO - neat way to change the docstrings...

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
