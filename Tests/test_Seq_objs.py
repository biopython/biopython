# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unittests for the Seq objects."""
import unittest
import sys
if sys.version_info[0] == 3:
   maketrans = str.maketrans
else:
   from string import maketrans

from Bio.Alphabet import generic_protein, generic_nucleotide, \
                         generic_dna, generic_rna
from Bio.Alphabet.IUPAC import protein, extended_protein
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna, ambiguous_rna
from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_rna_values
from Bio.Seq import Seq, UnknownSeq, MutableSeq, translate
from Bio.Data.CodonTable import TranslationError, CodonTable

#This is just the standard table with less stop codons
#(replaced with coding for O as an artifical example)
special_table = CodonTable(forward_table={
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': 'O',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'O', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
    start_codons=['TAA', 'TAG', 'TGA'],
    stop_codons=['TAG'])

Chilodonella_uncinata_table = CodonTable(forward_table={
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y',             'TAG': 'Q', 
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
    start_codons = [ 'ATG'],
    stop_codons = ['TAA' ])

class StringMethodTests(unittest.TestCase):
    _examples = [ \
        Seq("ACGTGGGGT", generic_protein),
        Seq("ACGTGGGGT", generic_nucleotide),
        Seq("ACGTGGGGT", generic_dna),
        Seq("ACGUGGGGU", generic_rna),
        Seq("GG", generic_protein),
        Seq("GG", generic_nucleotide),
        Seq("GG", generic_dna),
        Seq("GG", generic_rna),
        Seq("A", generic_protein),
        Seq("A", generic_nucleotide),
        Seq("A", generic_dna),
        Seq("A", generic_rna),
        UnknownSeq(1),
        UnknownSeq(1, character="n"),
        UnknownSeq(1, generic_rna),
        UnknownSeq(1, generic_rna, "n"),
        UnknownSeq(1, generic_rna, "N"),
        UnknownSeq(10, generic_rna, "N"),
        UnknownSeq(10, generic_dna, "N"),
        UnknownSeq(10, generic_nucleotide, "N"),
        UnknownSeq(10, generic_protein, "X"),
        UnknownSeq(10, character="X"),
        UnknownSeq(10),
        ]
    for seq in _examples[:]:
        if isinstance(seq, Seq):
            _examples.append(seq.tomutable())
    _start_end_values = [0, 1, 2, 1000, -1, -2, -999]


    def _test_method(self, method_name, pre_comp_function=None, start_end=False):
        """Check this method matches the plain string's method."""
        self.assertTrue(isinstance(method_name, str))
        for example1 in self._examples:
            if not hasattr(example1, method_name):
                #e.g. MutableSeq does not support find
                continue
            str1 = str(example1)

            for example2 in self._examples:
                if not hasattr(example2, method_name):
                    #e.g. MutableSeq does not support find
                    continue
                str2 = str(example2)

                i = getattr(example1,method_name)(str2)
                j = getattr(str1,method_name)(str2)
                if pre_comp_function:
                    i = pre_comp_function(i)
                    j = pre_comp_function(j)
                if i != j:
                    raise ValueError("%s.%s(%s) = %i, not %i" \
                                     % (repr(example1),
                                        method_name,
                                        repr(str2),
                                        i,
                                        j))

                try:
                    i = getattr(example1,method_name)(example2)
                    j = getattr(str1,method_name)(str2)
                    if pre_comp_function:
                        i = pre_comp_function(i)
                        j = pre_comp_function(j)
                    if i != j:
                        raise ValueError("%s.%s(%s) = %i, not %i" \
                                         % (repr(example1),
                                            method_name,
                                            repr(example2),
                                            i,
                                            j))
                except TypeError:
                    #TODO - Check the alphabets do clash!
                    pass

                if start_end:
                    for start in self._start_end_values:
                        i = getattr(example1,method_name)(str2, start)
                        j = getattr(str1,method_name)(str2, start)
                        if pre_comp_function:
                            i = pre_comp_function(i)
                            j = pre_comp_function(j)
                        if i != j:
                            raise ValueError("%s.%s(%s, %i) = %i, not %i" \
                                             % (repr(example1),
                                                method_name,
                                                repr(str2),
                                                start,
                                                i,
                                                j))
                        
                        for end in self._start_end_values:
                            i = getattr(example1,method_name)(str2, start, end)
                            j = getattr(str1,method_name)(str2, start, end)
                            if pre_comp_function:
                                i = pre_comp_function(i)
                                j = pre_comp_function(j)
                            if i != j:
                                raise ValueError("%s.%s(%s, %i, %i) = %i, not %i" \
                                                 % (repr(example1),
                                                    method_name,
                                                    repr(str2),
                                                    start,
                                                    end,
                                                    i,
                                                    j))

    def test_str_count(self):
        """Check matches the python string count method."""
        self._test_method("count", start_end=True)

    def test_str_find(self):
        """Check matches the python string find method."""
        self._test_method("find", start_end=True)

    def test_str_rfind(self):
        """Check matches the python string rfind method."""
        self._test_method("rfind", start_end=True)

    def test_str_startswith(self):
        """Check matches the python string startswith method."""
        self._test_method("startswith", start_end=True)

        try:
            self.assertTrue("ABCDE".startswith(("ABE","OBE","ABC")))
        except TypeError:
            #Base string only supports this on Python 2.5+, skip this
            return
        
        #Now check with a tuple of sub sequences
        for example1 in self._examples:
            if not hasattr(example1, "startswith"):
                #e.g. MutableSeq does not support this
                continue
            subs = tuple([example1[start:start+2] for start \
                          in range(0, len(example1)-2,3)])
            subs_str = tuple([str(s) for s in subs])

            self.assertEqual(str(example1).startswith(subs_str),
                             example1.startswith(subs))
            self.assertEqual(str(example1).startswith(subs_str),
                             example1.startswith(subs_str)) #strings!
            self.assertEqual(str(example1).startswith(subs_str,3),
                             example1.startswith(subs,3))
            self.assertEqual(str(example1).startswith(subs_str,2,6),
                             example1.startswith(subs,2,6))        

    def test_str_endswith(self):
        """Check matches the python string endswith method."""
        self._test_method("endswith", start_end=True)

        try:
            self.assertTrue("ABCDE".endswith(("ABE","OBE","CDE")))
        except TypeError:
            #Base string only supports this on Python 2.5+, skip this
            return

        #Now check with a tuple of sub sequences
        for example1 in self._examples:
            if not hasattr(example1, "endswith"):
                #e.g. MutableSeq does not support this
                continue
            subs = tuple([example1[start:start+2] for start \
                          in range(0, len(example1)-2,3)])
            subs_str = tuple([str(s) for s in subs])

            self.assertEqual(str(example1).endswith(subs_str),
                             example1.endswith(subs))
            self.assertEqual(str(example1).startswith(subs_str),
                             example1.startswith(subs_str)) #strings!
            self.assertEqual(str(example1).endswith(subs_str,3),
                             example1.endswith(subs,3))
            self.assertEqual(str(example1).endswith(subs_str,2,6),
                             example1.endswith(subs,2,6))

    def test_str_strip(self):
        """Check matches the python string strip method."""
        self._test_method("strip", pre_comp_function=str)

    def test_str_rstrip(self):
        """Check matches the python string rstrip method."""
        self._test_method("rstrip", pre_comp_function=str)

    def test_str_split(self):
        """Check matches the python string rstrip method."""
        #Calling (r)split should return a list of Seq-like objects, we'll
        #just apply str() to each of them so it matches the string method
        self._test_method("rstrip", pre_comp_function=lambda x : map(str,x))

    def test_str_rsplit(self):
        """Check matches the python string rstrip method."""
        #Calling (r)split should return a list of Seq-like objects, we'll
        #just apply str() to each of them so it matches the string method
        self._test_method("rstrip", pre_comp_function=lambda x : map(str,x))

    def test_str_lsplit(self):
        """Check matches the python string rstrip method."""
        #Calling (r)split should return a list of Seq-like objects, we'll
        #just apply str() to each of them so it matches the string method
        self._test_method("rstrip", pre_comp_function=lambda x : map(str,x))

    def test_str_length(self):
        """Check matches the python string __len__ method."""
        for example1 in self._examples:
            str1 = str(example1)
            self.assertEqual(len(example1), len(str1))

    def test_str_upper(self):
        """Check matches the python string upper method."""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq) : continue
            str1 = str(example1)
            self.assertEqual(str(example1.upper()), str1.upper())

    def test_str_upper(self):
        """Check matches the python string lower method."""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq) : continue
            str1 = str(example1)
            self.assertEqual(str(example1.lower()), str1.lower())

    def test_str_getitem(self):
        """Check slicing and indexing works like a string."""
        for example1 in self._examples:
            str1 = str(example1)
            for i in self._start_end_values:
                if abs(i) < len(example1):
                    self.assertEqual(str(example1[i]), str1[i])
                self.assertEqual(str(example1[:i]), str1[:i])
                self.assertEqual(str(example1[i:]), str1[i:])
                for j in self._start_end_values:
                    self.assertEqual(str(example1[i:j]), str1[i:j])
                    for step in range(-3,4):
                        if step == 0:
                            try:
                                print example1[i:j:step]
                                self._assert(False) #Should fail!
                            except ValueError:
                                pass
                        else:
                            self.assertEqual(str(example1[i:j:step]), \
                                             str1[i:j:step])

    def test_tostring(self):
        """Check str(obj) and obj.tostring() match."""
        for example1 in self._examples:
            str1 = str(example1)
            self.assertEqual(example1.tostring(), str1)

    def test_tomutable(self):
        """Check obj.tomutable() method."""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq) : continue
            mut = example1.tomutable()
            self.assertTrue(isinstance(mut, MutableSeq))
            self.assertEqual(str(mut), str(example1))
            self.assertEqual(mut.alphabet, example1.alphabet)

    def test_toseq(self):
        """Check obj.toseq() method."""
        for example1 in self._examples:
            try :
                seq = example1.toseq()
            except AttributeError :
                self.assertTrue(isinstance(example1, Seq))
                continue
            self.assertTrue(isinstance(seq, Seq))
            self.assertEqual(str(seq), str(example1))
            self.assertEqual(seq.alphabet, example1.alphabet)

    def test_the_complement(self):
        """Check obj.complement() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq) : continue
            try :
                comp = example1.complement()
            except ValueError, e:
                self.assertEqual(str(e), "Proteins do not have complements!")
                continue
            str1 = str(example1)
            #This only does the unambiguous cases
            if "U" in str1 or "u" in str1 \
            or example1.alphabet==generic_rna:
                mapping = maketrans("ACGUacgu","UGCAugca")
            elif "T" in str1 or "t" in str1 \
            or example1.alphabet==generic_dna \
            or example1.alphabet==generic_nucleotide:
                mapping = maketrans("ACGTacgt","TGCAtgca")
            elif "A" not in str1 and "a" not in str1:
                mapping = maketrans("CGcg","GCgc")
            else :
                #TODO - look at alphabet?
                raise ValueError(example1)
            self.assertEqual(str1.translate(mapping), str(comp))
            self.assertEqual(comp.alphabet, example1.alphabet)
                
    def test_the_reverse_complement(self):
        """Check obj.reverse_complement() method."""
        mapping = ""
        for example1 in self._examples:
            if isinstance(example1, MutableSeq) : continue
            try :
                comp = example1.reverse_complement()
            except ValueError, e:
                self.assertEqual(str(e), "Proteins do not have complements!")
                continue
            str1 = str(example1)
            #This only does the unambiguous cases
            if "U" in str1 or "u" in str1 \
            or example1.alphabet==generic_rna:
                mapping = maketrans("ACGUacgu","UGCAugca")
            elif "T" in str1 or "t" in str1 \
            or example1.alphabet==generic_dna \
            or example1.alphabet==generic_nucleotide:
                mapping = maketrans("ACGTacgt","TGCAtgca")
            elif "A" not in str1 and "a" not in str1:
                mapping = maketrans("CGcg","GCgc")
            else :
                #TODO - look at alphabet?
                continue
            self.assertEqual(str1.translate(mapping)[::-1], str(comp))
            self.assertEqual(comp.alphabet, example1.alphabet)

    def test_the_transcription(self):
            """Check obj.transcribe() method."""
            mapping = ""
            for example1 in self._examples:
                if isinstance(example1, MutableSeq) : continue
                try :
                    tran = example1.transcribe()
                except ValueError, e:
                    if str(e) == "Proteins cannot be transcribed!" : continue
                    if str(e) == "RNA cannot be transcribed!" : continue
                    raise e
                str1 = str(example1)
                self.assertEqual(str1.replace("T","U").replace("t","u"), str(tran))
                self.assertEqual(tran.alphabet, generic_rna) #based on limited examples             

    def test_the_back_transcription(self):
            """Check obj.back_transcribe() method."""
            mapping = ""
            for example1 in self._examples:
                if isinstance(example1, MutableSeq) : continue
                try :
                    tran = example1.back_transcribe()
                except ValueError, e:
                    if str(e) == "Proteins cannot be back transcribed!" : continue
                    if str(e) == "DNA cannot be back transcribed!" : continue
                    raise e
                str1 = str(example1)
                self.assertEqual(str1.replace("U","T").replace("u","t"), str(tran))
                self.assertEqual(tran.alphabet, generic_dna) #based on limited examples             

    def test_the_translate(self):
            """Check obj.translate() method."""
            mapping = ""
            for example1 in self._examples:
                if isinstance(example1, MutableSeq) : continue
                try :
                    tran = example1.translate()
                except ValueError, e:
                    if str(e) == "Proteins cannot be translated!" : continue
                    raise e
                #This is based on the limited example not having stop codons:
                if tran.alphabet not in [extended_protein, protein, generic_protein]:
                    print tran.alphabet
                    self.assertTrue(False)
                #TODO - check the actual translation, and all the optional args

    def test_the_translation_of_stops(self):
        """Check obj.translate() method with stop codons."""
        misc_stops = "TAATAGTGAAGAAGG"
        for nuc in [Seq(misc_stops),
                    Seq(misc_stops, generic_nucleotide),
                    Seq(misc_stops, generic_dna),
                    Seq(misc_stops, unambiguous_dna)]:
            self.assertEqual("***RR", str(nuc.translate()))
            self.assertEqual("***RR", str(nuc.translate(1)))
            self.assertEqual("***RR", str(nuc.translate("SGC0")))
            self.assertEqual("**W**", str(nuc.translate(table=2)))
            self.assertEqual("**WRR", str(nuc.translate(table='Yeast Mitochondrial')))
            self.assertEqual("**WSS", str(nuc.translate(table=5)))
            self.assertEqual("**WSS", str(nuc.translate(table=9)))
            self.assertEqual("**CRR", str(nuc.translate(table='Euplotid Nuclear')))
            self.assertEqual("***RR", str(nuc.translate(table=11)))
            self.assertEqual("***RR", str(nuc.translate(table='11')))
            self.assertEqual("***RR", str(nuc.translate(table='Bacterial')))
            self.assertEqual("", str(nuc.translate(to_stop=True)))
            self.assertEqual("O*ORR", str(nuc.translate(table=special_table)))
            self.assertEqual("*QWRR", str(nuc.translate(table=Chilodonella_uncinata_table)))
            #These test the Bio.Seq.translate() function - move these?:
            self.assertEqual("*QWRR", translate(str(nuc), table=Chilodonella_uncinata_table))
            self.assertEqual("O*ORR", translate(str(nuc), table=special_table))
            self.assertEqual("", translate(str(nuc), to_stop=True))
            self.assertEqual("***RR", translate(str(nuc), table='Bacterial'))
            self.assertEqual("***RR", translate(str(nuc), table='11'))
            self.assertEqual("***RR", translate(str(nuc), table=11))
            self.assertEqual("**W**", translate(str(nuc), table=2))
        self.assertEqual(str(Seq("TAT").translate()), "Y")
        self.assertEqual(str(Seq("TAR").translate()), "*")
        self.assertEqual(str(Seq("TAN").translate()), "X")
        self.assertEqual(str(Seq("NNN").translate()), "X")
        self.assertEqual(str(Seq("TAt").translate()), "Y")
        self.assertEqual(str(Seq("TaR").translate()), "*")
        self.assertEqual(str(Seq("TaN").translate()), "X")
        self.assertEqual(str(Seq("nnN").translate()), "X")
        self.assertEqual(str(Seq("tat").translate()), "Y")
        self.assertEqual(str(Seq("tar").translate()), "*")
        self.assertEqual(str(Seq("tan").translate()), "X")
        self.assertEqual(str(Seq("nnn").translate()), "X")


    def test_the_translation_of_invalid_codons(self):
        """Check obj.translate() method with invalid codons."""
        for codon in ["TA?", "N-N", "AC_", "Ac_"]:
            for nuc in [Seq(codon),
                        Seq(codon, generic_nucleotide),
                        Seq(codon, generic_dna),
                        Seq(codon, unambiguous_dna)]:
                try :
                    print nuc.translate()
                    self.assertTrue(False, "Transating %s should fail" % codon)
                except TranslationError :
                    pass

    def test_the_translation_of_ambig_codons(self):
        """Check obj.translate() method with ambiguous codons."""
        for letters, ambig_values in [(ambiguous_dna.letters, ambiguous_dna_values),
                                      (ambiguous_rna.letters, ambiguous_rna_values)] :
            ambig = set(letters)
            for c1 in ambig:
                for c2 in ambig:
                    for c3 in ambig:
                        values = set([str(Seq(a+b+c).translate()) \
                                      for a in ambig_values[c1] \
                                      for b in ambig_values[c2] \
                                      for c in ambig_values[c3]])
                        t = str(Seq(c1+c2+c3).translate())
                        if t=="*":
                            self.assertEqual(values, set("*"))
                        elif t=="X":
                            self.assertTrue(len(values) > 1, \
                                "translate('%s') = '%s' not '%s'" \
                                % (c1+c2+c3, t, ",".join(values)))
                        elif t=="Z":
                            self.assertEqual(values, set("EQ"))
                        elif t=="B":
                            self.assertEqual(values, set("DN"))
                        elif t=="J":
                            self.assertEqual(values, set("LI"))
                        else:
                            self.assertEqual(values, set(t))
                        #TODO - Use the Bio.Data.IUPACData module for the
                        #ambiguous protein mappings?

    def test_init_typeerror(self):
        """Check Seq __init__ gives TypeError exceptions."""
        #Only expect it to take strings and unicode - not Seq objects!
        self.assertRaises(TypeError, Seq, (1066))
        self.assertRaises(TypeError, Seq, (Seq("ACGT", generic_dna)))

    #TODO - Addition...

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
