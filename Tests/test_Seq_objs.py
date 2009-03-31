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
    for seq in _examples[:] :
        if isinstance(seq, Seq) :
            _examples.append(seq.tomutable())
    _start_end_values = [0, 1, 2, 1000, -1, -2, -999]


    def _test_method(self, method_name, apply_str=False, start_end=False) :
        """Check this method matches the plain string's method."""
        assert isinstance(method_name, str)
        for example1 in self._examples :
            if not hasattr(example1, method_name) :
                #e.g. MutableSeq does not support find
                continue
            str1 = str(example1)

            for example2 in self._examples :
                if not hasattr(example2, method_name) :
                    #e.g. MutableSeq does not support find
                    continue
                str2 = str(example2)

                i = getattr(example1,method_name)(str2)
                j = getattr(str1,method_name)(str2)
                if apply_str :
                    i = str(i)
                    j = str(j)
                if i != j :
                    raise ValueError("%s.%s(%s) = %i, not %i" \
                                     % (repr(example1),
                                        method_name,
                                        repr(str2),
                                        i,
                                        j))

                try :
                    i = getattr(example1,method_name)(example2)
                    j = getattr(str1,method_name)(str2)
                    if apply_str :
                        i = str(i)
                        j = str(j)
                    if i != j :
                        raise ValueError("%s.%s(%s) = %i, not %i" \
                                         % (repr(example1),
                                            method_name,
                                            repr(example2),
                                            i,
                                            j))
                except TypeError :
                    #TODO - Check the alphabets do clash!
                    pass

                if start_end :
                    for start in self._start_end_values :
                        i = getattr(example1,method_name)(str2, start)
                        j = getattr(str1,method_name)(str2, start)
                        if apply_str :
                            i = str(i)
                            j = str(j)
                        if i != j :
                            raise ValueError("%s.%s(%s, %i) = %i, not %i" \
                                             % (repr(example1),
                                                method_name,
                                                repr(str2),
                                                start,
                                                i,
                                                j))
                        
                        for end in self._start_end_values :
                            i = getattr(example1,method_name)(str2, start, end)
                            j = getattr(str1,method_name)(str2, start, end)
                            if apply_str :
                                i = str(i)
                                j = str(j)
                            if i != j :
                                raise ValueError("%s.%s(%s, %i, %i) = %i, not %i" \
                                                 % (repr(example1),
                                                    method_name,
                                                    repr(str2),
                                                    start,
                                                    end,
                                                    i,
                                                    j))

    def test_count(self) :
        """Check matches the python string count method."""
        self._test_method("count", start_end=True)

    def test_find(self) :
        """Check matches the python string find method."""
        self._test_method("find", start_end=True)

    def test_rfind(self) :
        """Check matches the python string rfind method."""
        self._test_method("rfind", start_end=True)

    def test_strip(self) :
        """Check matches the python string strip method."""
        self._test_method("strip", apply_str=True)

    def test_rstrip(self) :
        """Check matches the python string rstrip method."""
        self._test_method("rstrip", apply_str=True)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
