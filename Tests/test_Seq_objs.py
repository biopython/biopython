# Copyright 2009 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Unittests for the Seq objects."""
import unittest

from Bio.Alphabet import generic_protein, generic_nucleotide, \
                         generic_dna, generic_rna
from Bio.Seq import Seq, UnknownSeq, MutableSeq

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
        assert isinstance(method_name, str)
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

    def test_count(self):
        """Check matches the python string count method."""
        self._test_method("count", start_end=True)

    def test_find(self):
        """Check matches the python string find method."""
        self._test_method("find", start_end=True)

    def test_rfind(self):
        """Check matches the python string rfind method."""
        self._test_method("rfind", start_end=True)

    def test_startswith(self):
        """Check matches the python string startswith method."""
        self._test_method("startswith", start_end=True)

        try:
            self.assert_("ABCDE".startswith(("ABE","OBE","ABC")))
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

    def test_endswith(self):
        """Check matches the python string endswith method."""
        self._test_method("endswith", start_end=True)

        try:
            self.assert_("ABCDE".endswith(("ABE","OBE","CDE")))
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

    def test_strip(self):
        """Check matches the python string strip method."""
        self._test_method("strip", pre_comp_function=str)

    def test_rstrip(self):
        """Check matches the python string rstrip method."""
        self._test_method("rstrip", pre_comp_function=str)

    def test_split(self):
        """Check matches the python string rstrip method."""
        #Calling (r)split should return a list of Seq-like objects, we'll
        #just apply str() to each of them so it matches the string method
        self._test_method("rstrip", pre_comp_function=lambda x : map(str,x))

    def test_rsplit(self):
        """Check matches the python string rstrip method."""
        #Calling (r)split should return a list of Seq-like objects, we'll
        #just apply str() to each of them so it matches the string method
        self._test_method("rstrip", pre_comp_function=lambda x : map(str,x))

    def test_lsplit(self):
        """Check matches the python string rstrip method."""
        #Calling (r)split should return a list of Seq-like objects, we'll
        #just apply str() to each of them so it matches the string method
        self._test_method("rstrip", pre_comp_function=lambda x : map(str,x))

    def test_length(self):
        """Check matches the python string __len__ method."""
        for example1 in self._examples:
            str1 = str(example1)
            self.assertEqual(len(example1), len(str1))

    def test_getitem(self):
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
        """Check str(obj) and obj.tostring() match.

        Also check the obj.data attribute for non-MutableSeq objects."""
        for example1 in self._examples:
            str1 = str(example1)
            self.assertEqual(example1.tostring(), str1)
            if not isinstance(example1, MutableSeq):
                self.assertEqual(example1.data, str1)

    #TODO - Addition...

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
