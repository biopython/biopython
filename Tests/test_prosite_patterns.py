# Copyright 1999 by Cayte Lindner.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
import unittest

if os.name == 'java':
    try :
        buffer
    except NameError :
        from Bio import MissingExternalDependencyError
        #This is a slight miss-use of MissingExternalDependencyError,
        #but it will do in the short term to skip this unit test on Jython
        raise MissingExternalDependencyError(\
            "The (deprecated) Bio.Prosite.Pattern module uses Python "
            "function buffer which is not supported on Jython, see "
            "http://bugs.jython.org/issue1521")

if sys.version_info[0] >= 3:
    #This is a slight miss-use of MissingExternalDependencyError,
    #but it will do in the short term to skip this unit test on Python 3
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(\
        "This deprecated module doesn't work on Python 3.")
        
import warnings
from Bio import BiopythonDeprecationWarning
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
from Bio.Prosite import Pattern
warnings.filters.pop()

from Bio import Seq

class TestPrositePattern(unittest.TestCase):

    def test_pattern01(self):
        "Testing Prosite pattern 'A.'"
        p = Pattern.Prosite(pattern = "A.")
        self.assertEqual(repr(p.re.pattern), "'A'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)'")
        self.assertEqual(p.tostring(), "A.")
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.letters, 'A')
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, 'A')
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)
        # Bad pattern BBP
        self.assertEqual(p.re.search("BBP"), None)
        self.assertEqual(p.search(Seq.Seq("BBP")), None)

    def test_pattern02(self):
        "Testing Prosite pattern 'x.'"
        p = Pattern.Prosite(pattern = "x.")
        self.assertEqual(repr(p.re.pattern), "'.'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(.)'")
        self.assertEqual(p.tostring(), "x.")
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, 'x')
        # Good pattern BA
        self.assertNotEqual(p.re.search('BA'), None)
        m = p.search(Seq.Seq('BA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, 'x')
        # Good pattern B
        self.assertNotEqual(p.re.search('B'), None)
        m = p.search(Seq.Seq('B'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, 'x')
        # Bad pattern (empty string)
        self.assertEqual(p.re.search(""), None)
        self.assertEqual(p.search(Seq.Seq("")), None)

    def test_pattern03(self):
        "Testing Prosite pattern '{A}.'"
        p = Pattern.Prosite(pattern = "{A}.")
        self.assertEqual(repr(p.re.pattern), "'[^A]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^A])'")
        self.assertEqual(p.tostring(), "{A}.")
        # Good pattern B
        self.assertNotEqual(p.re.search("B"), None)
        m = p.search(Seq.Seq("B"))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1)
        self.assertEqual(mpi.letters, "A")
        # Good pattern BBP
        self.assertNotEqual(p.re.search("BBP"), None)
        m = p.search(Seq.Seq("BBP"))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1)
        self.assertEqual(mpi.letters, "A")
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern AAA
        self.assertEqual(p.re.search("AAA"), None)
        self.assertEqual(p.search(Seq.Seq("AAA")), None)
        # Bad pattern (empty string)
        self.assertEqual(p.re.search(""), None)
        self.assertEqual(p.search(Seq.Seq("")), None)

    def test_pattern04(self):
        "Testing Prosite pattern '[AB].'"
        p = Pattern.Prosite(pattern = "[AB].")
        self.assertEqual(repr(p.re.pattern), "'[AB]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])'")
        self.assertEqual(p.tostring(), "[AB].")
        # Good pattern AB
        self.assertNotEqual(p.re.search("AB"), None)
        m = p.search(Seq.Seq("AB"))
        self.assertNotEqual(m, None)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, "AB")
        # Good pattern PA
        self.assertNotEqual(p.re.search("PA"), None)
        m = p.search(Seq.Seq("PA"))
        self.assertNotEqual(m, None)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, "AB")
        # Good pattern PPAB
        self.assertNotEqual(p.re.search("PPAB"), None)
        m = p.search(Seq.Seq("PPAB"))
        self.assertNotEqual(m, None)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, "AB")
        # Good pattern PPABP
        self.assertNotEqual(p.re.search("PPABP"), None)
        m = p.search(Seq.Seq("PPABP"))
        self.assertNotEqual(m, None)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, "AB")
        # Good pattern B
        self.assertNotEqual(p.re.search("B"), None)
        m = p.search(Seq.Seq("B"))
        self.assertNotEqual(m, None)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0)
        self.assertEqual(mpi.letters, "AB")
        # Bad pattern C
        self.assertEqual(p.re.search("C"), None)
        self.assertEqual(p.search(Seq.Seq("C")), None)
        # Bad pattern CCCCCCFDS
        self.assertEqual(p.re.search("CCCCCCFDS"), None)
        self.assertEqual(p.search(Seq.Seq("CCCCCCFDS")), None)

    def test_pattern05(self):
        "Testing Prosite pattern '{AB}.'"
        p = Pattern.Prosite(pattern = "{AB}.")
        self.assertEqual(repr(p.re.pattern), "'[^AB]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^AB])'")
        self.assertEqual(p.tostring(), "{AB}.")
        # Good pattern C
        self.assertNotEqual(p.re.search("C"), None)
        m = p.search(Seq.Seq("C"))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1)
        self.assertEqual(mpi.letters, "AB")
        # Good pattern CDFDC
        self.assertNotEqual(p.re.search("CDFDC"), None)
        m = p.search(Seq.Seq("CDFDC"))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1)
        self.assertEqual(mpi.letters, "AB")
        # Bad pattern (empty string)
        self.assertEqual(p.re.search(""), None)
        self.assertEqual(p.search(Seq.Seq("")), None)
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)
        # Bad pattern BAB
        self.assertEqual(p.re.search("BAB"), None)
        self.assertEqual(p.search(Seq.Seq("BAB")), None)

    def test_pattern06(self):
        "Testing Prosite pattern 'A-B.'"
        p = Pattern.Prosite(pattern = "A-B.")
        self.assertEqual(repr(p.re.pattern), "'AB'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)(B)'")
        self.assertEqual(p.tostring(), "A-B.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern AACAB
        self.assertNotEqual(p.re.search('AACAB'), None)
        m = p.search(Seq.Seq('AACAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 3)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 3)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern ABAAA
        self.assertNotEqual(p.re.search('ABAAA'), None)
        m = p.search(Seq.Seq('ABAAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)
        # Bad pattern PPDAAPB
        self.assertEqual(p.re.search("PPDAAPB"), None)
        self.assertEqual(p.search(Seq.Seq("PPDAAPB")), None)
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)

    def test_pattern07(self):
        "Testing Prosite pattern 'A-x.'"
        p = Pattern.Prosite(pattern = "A-x.")
        self.assertEqual(repr(p.re.pattern), "'A.'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)(.)'")
        self.assertEqual(p.tostring(), "A-x.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        # Good pattern AACAB
        self.assertNotEqual(p.re.search('AACAB'), None)
        m = p.search(Seq.Seq('AACAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        # Good pattern ABAAA
        self.assertNotEqual(p.re.search('ABAAA'), None)
        m = p.search(Seq.Seq('ABAAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)

    def test_pattern08(self):
        "Testing Prosite pattern '[AB]-[BC].'"
        p = Pattern.Prosite(pattern = "[AB]-[BC].")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])'")
        self.assertEqual(p.tostring(), "[AB]-[BC].")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BB
        self.assertNotEqual(p.re.search('BB'), None)
        m = p.search(Seq.Seq('BB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BC
        self.assertNotEqual(p.re.search('BC'), None)
        m = p.search(Seq.Seq('BC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)
        # Bad pattern CA
        self.assertEqual(p.re.search("CA"), None)
        self.assertEqual(p.search(Seq.Seq("CA")), None)
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern C
        self.assertEqual(p.re.search("C"), None)
        self.assertEqual(p.search(Seq.Seq("C")), None)
        # Bad pattern PPQ
        self.assertEqual(p.re.search("PPQ"), None)
        self.assertEqual(p.search(Seq.Seq("PPQ")), None)

    def test_pattern09(self):
        "Testing Prosite pattern '[AB]-[BC]-[CD].'"
        p = Pattern.Prosite(pattern = "[AB]-[BC]-[CD].")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC][CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])([CD])'")
        self.assertEqual(p.tostring(), "[AB]-[BC]-[CD].")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QACD
        self.assertNotEqual(p.re.search('QACD'), None)
        m = p.search(Seq.Seq('QACD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QBBCQ
        self.assertNotEqual(p.re.search('QBBCQ'), None)
        m = p.search(Seq.Seq('QBBCQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABA
        self.assertEqual(p.re.search("ABA"), None)
        self.assertEqual(p.search(Seq.Seq("ABA")), None)
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern DCB
        self.assertEqual(p.re.search("DCB"), None)
        self.assertEqual(p.search(Seq.Seq("DCB")), None)

    def test_pattern10(self):
        "Testing Prosite pattern '[AB]-[BC]-[CDEFGHIKLMN].'"
        p = Pattern.Prosite(pattern = "[AB]-[BC]-[CDEFGHIKLMN].")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC][CDEFGHIKLMN]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])([CDEFGHIKLMN])'")
        self.assertEqual(p.tostring(), "[AB]-[BC]-[CDEFGHIKLMN].")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Good pattern ACN
        self.assertNotEqual(p.re.search('ACN'), None)
        m = p.search(Seq.Seq('ACN'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)
        # Bad pattern PERPPWPTW
        self.assertEqual(p.re.search("PERPPWPTW"), None)
        self.assertEqual(p.search(Seq.Seq("PERPPWPTW")), None)

    def test_pattern11(self):
        "Testing Prosite pattern '{AB}-[BC]-[CD].'"
        p = Pattern.Prosite(pattern = "{AB}-[BC]-[CD].")
        self.assertEqual(repr(p.re.pattern), "'[^AB][BC][CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^AB])([BC])([CD])'")
        self.assertEqual(p.tostring(), "{AB}-[BC]-[CD].")
        # Good pattern CCC
        self.assertNotEqual(p.re.search('CCC'), None)
        m = p.search(Seq.Seq('CCC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern CCD
        self.assertNotEqual(p.re.search('CCD'), None)
        m = p.search(Seq.Seq('CCD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QCCCP
        self.assertNotEqual(p.re.search('QCCCP'), None)
        m = p.search(Seq.Seq('QCCCP'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABC
        self.assertEqual(p.re.search("ABC"), None)
        self.assertEqual(p.search(Seq.Seq("ABC")), None)

    def test_pattern12(self):
        "Testing Prosite pattern '[AB]-{BC}-[CD].'"
        p = Pattern.Prosite(pattern = "[AB]-{BC}-[CD].")
        self.assertEqual(repr(p.re.pattern), "'[AB][^BC][CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([^BC])([CD])'")
        self.assertEqual(p.tostring(), "[AB]-{BC}-[CD].")
        # Good pattern AAD
        self.assertNotEqual(p.re.search('AAD'), None)
        m = p.search(Seq.Seq('AAD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern APD
        self.assertNotEqual(p.re.search('APD'), None)
        m = p.search(Seq.Seq('APD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QBACP
        self.assertNotEqual(p.re.search('QBACP'), None)
        m = p.search(Seq.Seq('QBACP'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABC
        self.assertEqual(p.re.search("ABC"), None)
        self.assertEqual(p.search(Seq.Seq("ABC")), None)
        # Bad pattern PAAEQ
        self.assertEqual(p.re.search("PAAEQ"), None)
        self.assertEqual(p.search(Seq.Seq("PAAEQ")), None)

    def test_pattern13(self):
        "Testing Prosite pattern '[AB]-[BC]-{CD}.'"
        p = Pattern.Prosite(pattern = "[AB]-[BC]-{CD}.")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC][^CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])([^CD])'")
        self.assertEqual(p.tostring(), "[AB]-[BC]-{CD}.")
        # Good pattern ABA
        self.assertNotEqual(p.re.search('ABA'), None)
        m = p.search(Seq.Seq('ABA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern ABP
        self.assertNotEqual(p.re.search('ABP'), None)
        m = p.search(Seq.Seq('ABP'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QBBAP
        self.assertNotEqual(p.re.search('QBBAP'), None)
        m = p.search(Seq.Seq('QBBAP'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABC
        self.assertEqual(p.re.search("ABC"), None)
        self.assertEqual(p.search(Seq.Seq("ABC")), None)
        # Bad pattern PAAE
        self.assertEqual(p.re.search("PAAE"), None)
        self.assertEqual(p.search(Seq.Seq("PAAE")), None)

    def test_pattern14(self):
        "Testing Prosite pattern '{AB}-[BC]-{CD}.'"
        p = Pattern.Prosite(pattern = "{AB}-[BC]-{CD}.")
        self.assertEqual(repr(p.re.pattern), "'[^AB][BC][^CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^AB])([BC])([^CD])'")
        self.assertEqual(p.tostring(), "{AB}-[BC]-{CD}.")
        # Good pattern CCB
        self.assertNotEqual(p.re.search('CCB'), None)
        m = p.search(Seq.Seq('CCB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern PBP
        self.assertNotEqual(p.re.search('PBP'), None)
        m = p.search(Seq.Seq('PBP'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABA
        self.assertEqual(p.re.search("ABA"), None)
        self.assertEqual(p.search(Seq.Seq("ABA")), None)
        # Bad pattern PACR
        self.assertEqual(p.re.search("PACR"), None)
        self.assertEqual(p.search(Seq.Seq("PACR")), None)

    def test_pattern15(self):
        "Testing Prosite pattern 'A-B-C-D.'"
        p = Pattern.Prosite(pattern = "A-B-C-D.")
        self.assertEqual(repr(p.re.pattern), "'ABCD'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)(B)(C)(D)'")
        self.assertEqual(p.tostring(), "A-B-C-D.")
        # Good pattern QABCDQ
        self.assertNotEqual(p.re.search('QABCDQ'), None)
        m = p.search(Seq.Seq('QABCDQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 4)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'C') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'D') 
        # Bad pattern ABCP
        self.assertEqual(p.re.search("ABCP"), None)
        self.assertEqual(p.search(Seq.Seq("ABCP")), None)

    def test_pattern16(self):
        "Testing Prosite pattern '<A.'"
        p = Pattern.Prosite(pattern = "<A.")
        self.assertEqual(repr(p.re.pattern), "'^A'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^(A)'")
        self.assertEqual(p.tostring(), "<A.")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern PAB
        self.assertEqual(p.re.search("PAB"), None)
        self.assertEqual(p.search(Seq.Seq("PAB")), None)
        # Bad pattern BAAAAAA
        self.assertEqual(p.re.search("BAAAAAA"), None)
        self.assertEqual(p.search(Seq.Seq("BAAAAAA")), None)

    def test_pattern17(self):
        "Testing Prosite pattern '<{A}.'"
        p = Pattern.Prosite(pattern = "<{A}.")
        self.assertEqual(repr(p.re.pattern), "'^[^A]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^A])'")
        self.assertEqual(p.tostring(), "<{A}.")
        # Good pattern PAB
        self.assertNotEqual(p.re.search('PAB'), None)
        m = p.search(Seq.Seq('PAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern BAAAAAA
        self.assertNotEqual(p.re.search('BAAAAAA'), None)
        m = p.search(Seq.Seq('BAAAAAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern ABC
        self.assertEqual(p.re.search("ABC"), None)
        self.assertEqual(p.search(Seq.Seq("ABC")), None)
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)

    def test_pattern18(self):
        "Testing Prosite pattern '<[AB].'"
        p = Pattern.Prosite(pattern = "<[AB].")
        self.assertEqual(repr(p.re.pattern), "'^[AB]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])'")
        self.assertEqual(p.tostring(), "<[AB].")
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern B
        self.assertNotEqual(p.re.search('B'), None)
        m = p.search(Seq.Seq('B'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern AP
        self.assertNotEqual(p.re.search('AP'), None)
        m = p.search(Seq.Seq('AP'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern BBQ
        self.assertNotEqual(p.re.search('BBQ'), None)
        m = p.search(Seq.Seq('BBQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Bad pattern QBB
        self.assertEqual(p.re.search("QBB"), None)
        self.assertEqual(p.search(Seq.Seq("QBB")), None)
        # Bad pattern PA
        self.assertEqual(p.re.search("PA"), None)
        self.assertEqual(p.search(Seq.Seq("PA")), None)

    def test_pattern19(self):
        "Testing Prosite pattern '<{AB}.'"
        p = Pattern.Prosite(pattern = "<{AB}.")
        self.assertEqual(repr(p.re.pattern), "'^[^AB]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^AB])'")
        self.assertEqual(p.tostring(), "<{AB}.")
        # Good pattern QBB
        self.assertNotEqual(p.re.search('QBB'), None)
        m = p.search(Seq.Seq('QBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern PA
        self.assertNotEqual(p.re.search('PA'), None)
        m = p.search(Seq.Seq('PA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)
        # Bad pattern AP
        self.assertEqual(p.re.search("AP"), None)
        self.assertEqual(p.search(Seq.Seq("AP")), None)
        # Bad pattern BBQ
        self.assertEqual(p.re.search("BBQ"), None)
        self.assertEqual(p.search(Seq.Seq("BBQ")), None)

    def test_pattern20(self):
        "Testing Prosite pattern '<A-B.'"
        p = Pattern.Prosite(pattern = "<A-B.")
        self.assertEqual(repr(p.re.pattern), "'^AB'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^(A)(B)'")
        self.assertEqual(p.tostring(), "<A-B.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern ABB
        self.assertNotEqual(p.re.search('ABB'), None)
        m = p.search(Seq.Seq('ABB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Bad pattern BAB
        self.assertEqual(p.re.search("BAB"), None)
        self.assertEqual(p.search(Seq.Seq("BAB")), None)
        # Bad pattern PABQ
        self.assertEqual(p.re.search("PABQ"), None)
        self.assertEqual(p.search(Seq.Seq("PABQ")), None)

    def test_pattern21(self):
        "Testing Prosite pattern '<[AB]-[BC].'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC].")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])'")
        self.assertEqual(p.tostring(), "<[AB]-[BC].")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BB
        self.assertNotEqual(p.re.search('BB'), None)
        m = p.search(Seq.Seq('BB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BC
        self.assertNotEqual(p.re.search('BC'), None)
        m = p.search(Seq.Seq('BC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Bad pattern PAB
        self.assertEqual(p.re.search("PAB"), None)
        self.assertEqual(p.search(Seq.Seq("PAB")), None)
        # Bad pattern QBB
        self.assertEqual(p.re.search("QBB"), None)
        self.assertEqual(p.search(Seq.Seq("QBB")), None)
        # Bad pattern PPPPPQABBCAC
        self.assertEqual(p.re.search("PPPPPQABBCAC"), None)
        self.assertEqual(p.search(Seq.Seq("PPPPPQABBCAC")), None)

    def test_pattern22(self):
        "Testing Prosite pattern '<[AB]-[BC]-[CD].'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC]-[CD].")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC][CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])([CD])'")
        self.assertEqual(p.tostring(), "<[AB]-[BC]-[CD].")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern ACC
        self.assertNotEqual(p.re.search('ACC'), None)
        m = p.search(Seq.Seq('ACC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern CCD
        self.assertEqual(p.re.search("CCD"), None)
        self.assertEqual(p.search(Seq.Seq("CCD")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)
        # Bad pattern QPPAABC
        self.assertEqual(p.re.search("QPPAABC"), None)
        self.assertEqual(p.search(Seq.Seq("QPPAABC")), None)

    def test_pattern23(self):
        "Testing Prosite pattern '<[AB]-[BC]-[CDEFGHIKLMN].'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC]-[CDEFGHIKLMN].")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC][CDEFGHIKLMN]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])([CDEFGHIKLMN])'")
        self.assertEqual(p.tostring(), "<[AB]-[BC]-[CDEFGHIKLMN].")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Good pattern BBN
        self.assertNotEqual(p.re.search('BBN'), None)
        m = p.search(Seq.Seq('BBN'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Bad pattern QABC
        self.assertEqual(p.re.search("QABC"), None)
        self.assertEqual(p.search(Seq.Seq("QABC")), None)
        # Bad pattern PBBN
        self.assertEqual(p.re.search("PBBN"), None)
        self.assertEqual(p.search(Seq.Seq("PBBN")), None)
        # Bad pattern PPQ
        self.assertEqual(p.re.search("PPQ"), None)
        self.assertEqual(p.search(Seq.Seq("PPQ")), None)

    def test_pattern24(self):
        "Testing Prosite pattern '<{AB}-[BC]-[CD].'"
        p = Pattern.Prosite(pattern = "<{AB}-[BC]-[CD].")
        self.assertEqual(repr(p.re.pattern), "'^[^AB][BC][CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^AB])([BC])([CD])'")
        self.assertEqual(p.tostring(), "<{AB}-[BC]-[CD].")
        # Good pattern CCC
        self.assertNotEqual(p.re.search('CCC'), None)
        m = p.search(Seq.Seq('CCC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QBDQ
        self.assertNotEqual(p.re.search('QBDQ'), None)
        m = p.search(Seq.Seq('QBDQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABD
        self.assertEqual(p.re.search("ABD"), None)
        self.assertEqual(p.search(Seq.Seq("ABD")), None)
        # Bad pattern PDD
        self.assertEqual(p.re.search("PDD"), None)
        self.assertEqual(p.search(Seq.Seq("PDD")), None)

    def test_pattern25(self):
        "Testing Prosite pattern '<[AB]-{BC}-[CD].'"
        p = Pattern.Prosite(pattern = "<[AB]-{BC}-[CD].")
        self.assertEqual(repr(p.re.pattern), "'^[AB][^BC][CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([^BC])([CD])'")
        self.assertEqual(p.tostring(), "<[AB]-{BC}-[CD].")
        # Good pattern AEC
        self.assertNotEqual(p.re.search('AEC'), None)
        m = p.search(Seq.Seq('AEC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern APD
        self.assertNotEqual(p.re.search('APD'), None)
        m = p.search(Seq.Seq('APD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern CDD
        self.assertEqual(p.re.search("CDD"), None)
        self.assertEqual(p.search(Seq.Seq("CDD")), None)
        # Bad pattern AEE
        self.assertEqual(p.re.search("AEE"), None)
        self.assertEqual(p.search(Seq.Seq("AEE")), None)
        # Bad pattern ACC
        self.assertEqual(p.re.search("ACC"), None)
        self.assertEqual(p.search(Seq.Seq("ACC")), None)
        # Bad pattern QAEC
        self.assertEqual(p.re.search("QAEC"), None)
        self.assertEqual(p.search(Seq.Seq("QAEC")), None)

    def test_pattern26(self):
        "Testing Prosite pattern '<[AB]-[BC]-{CD}.'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC]-{CD}.")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC][^CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])([^CD])'")
        self.assertEqual(p.tostring(), "<[AB]-[BC]-{CD}.")
        # Good pattern ABB
        self.assertNotEqual(p.re.search('ABB'), None)
        m = p.search(Seq.Seq('ABB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABC
        self.assertEqual(p.re.search("ABC"), None)
        self.assertEqual(p.search(Seq.Seq("ABC")), None)
        # Bad pattern QABB
        self.assertEqual(p.re.search("QABB"), None)
        self.assertEqual(p.search(Seq.Seq("QABB")), None)

    def test_pattern27(self):
        "Testing Prosite pattern '<{AB}-[BC]-{CD}.'"
        p = Pattern.Prosite(pattern = "<{AB}-[BC]-{CD}.")
        self.assertEqual(repr(p.re.pattern), "'^[^AB][BC][^CD]'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^AB])([BC])([^CD])'")
        self.assertEqual(p.tostring(), "<{AB}-[BC]-{CD}.")
        # Good pattern CBB
        self.assertNotEqual(p.re.search('CBB'), None)
        m = p.search(Seq.Seq('CBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABA
        self.assertEqual(p.re.search("ABA"), None)
        self.assertEqual(p.search(Seq.Seq("ABA")), None)
        # Bad pattern QCDB
        self.assertEqual(p.re.search("QCDB"), None)
        self.assertEqual(p.search(Seq.Seq("QCDB")), None)

    def test_pattern28(self):
        "Testing Prosite pattern 'A>.'"
        p = Pattern.Prosite(pattern = "A>.")
        self.assertEqual(repr(p.re.pattern), "'A$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)$'")
        self.assertEqual(p.tostring(), "A>.")
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern BA
        self.assertNotEqual(p.re.search('BA'), None)
        m = p.search(Seq.Seq('BA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern AA
        self.assertNotEqual(p.re.search('AA'), None)
        m = p.search(Seq.Seq('AA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern QQQA
        self.assertNotEqual(p.re.search('QQQA'), None)
        m = p.search(Seq.Seq('QQQA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 3)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 3)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)
        # Bad pattern AAAAB
        self.assertEqual(p.re.search("AAAAB"), None)
        self.assertEqual(p.search(Seq.Seq("AAAAB")), None)

    def test_pattern29(self):
        "Testing Prosite pattern '{A}>.'"
        p = Pattern.Prosite(pattern = "{A}>.")
        self.assertEqual(repr(p.re.pattern), "'[^A]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^A])$'")
        self.assertEqual(p.tostring(), "{A}>.")
        # Good pattern B
        self.assertNotEqual(p.re.search('B'), None)
        m = p.search(Seq.Seq('B'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern AAAQ
        self.assertNotEqual(p.re.search('AAAQ'), None)
        m = p.search(Seq.Seq('AAAQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 3)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 3)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)
        # Bad pattern QQQQQQA
        self.assertEqual(p.re.search("QQQQQQA"), None)
        self.assertEqual(p.search(Seq.Seq("QQQQQQA")), None)

    def test_pattern30(self):
        "Testing Prosite pattern '[AB]>.'"
        p = Pattern.Prosite(pattern = "[AB]>.")
        self.assertEqual(repr(p.re.pattern), "'[AB]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])$'")
        self.assertEqual(p.tostring(), "[AB]>.")
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern B
        self.assertNotEqual(p.re.search('B'), None)
        m = p.search(Seq.Seq('B'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern QQA
        self.assertNotEqual(p.re.search('QQA'), None)
        m = p.search(Seq.Seq('QQA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 2)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)
        # Bad pattern AQ
        self.assertEqual(p.re.search("AQ"), None)
        self.assertEqual(p.search(Seq.Seq("AQ")), None)
        # Bad pattern BQ
        self.assertEqual(p.re.search("BQ"), None)
        self.assertEqual(p.search(Seq.Seq("BQ")), None)
        # Bad pattern ABQ
        self.assertEqual(p.re.search("ABQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABQ")), None)

    def test_pattern31(self):
        "Testing Prosite pattern '{AB}>.'"
        p = Pattern.Prosite(pattern = "{AB}>.")
        self.assertEqual(repr(p.re.pattern), "'[^AB]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^AB])$'")
        self.assertEqual(p.tostring(), "{AB}>.")
        # Good pattern Q
        self.assertNotEqual(p.re.search('Q'), None)
        m = p.search(Seq.Seq('Q'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern AQ
        self.assertNotEqual(p.re.search('AQ'), None)
        m = p.search(Seq.Seq('AQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern BQ
        self.assertNotEqual(p.re.search('BQ'), None)
        m = p.search(Seq.Seq('BQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern ABQ
        self.assertNotEqual(p.re.search('ABQ'), None)
        m = p.search(Seq.Seq('ABQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 2)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)
        # Bad pattern QQA
        self.assertEqual(p.re.search("QQA"), None)
        self.assertEqual(p.search(Seq.Seq("QQA")), None)

    def test_pattern32(self):
        "Testing Prosite pattern 'A-B>.'"
        p = Pattern.Prosite(pattern = "A-B>.")
        self.assertEqual(repr(p.re.pattern), "'AB$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)(B)$'")
        self.assertEqual(p.tostring(), "A-B>.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern QAB
        self.assertNotEqual(p.re.search('QAB'), None)
        m = p.search(Seq.Seq('QAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Bad pattern ABA
        self.assertEqual(p.re.search("ABA"), None)
        self.assertEqual(p.search(Seq.Seq("ABA")), None)
        # Bad pattern ABQ
        self.assertEqual(p.re.search("ABQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABQ")), None)
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)

    def test_pattern33(self):
        "Testing Prosite pattern '[AB]-[BC]>.'"
        p = Pattern.Prosite(pattern = "[AB]-[BC]>.")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])$'")
        self.assertEqual(p.tostring(), "[AB]-[BC]>.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BB
        self.assertNotEqual(p.re.search('BB'), None)
        m = p.search(Seq.Seq('BB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BC
        self.assertNotEqual(p.re.search('BC'), None)
        m = p.search(Seq.Seq('BC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Bad pattern ABQ
        self.assertEqual(p.re.search("ABQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABQ")), None)
        # Bad pattern ACQ
        self.assertEqual(p.re.search("ACQ"), None)
        self.assertEqual(p.search(Seq.Seq("ACQ")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)
        # Bad pattern QQQ
        self.assertEqual(p.re.search("QQQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQQ")), None)

    def test_pattern34(self):
        "Testing Prosite pattern '[AB]-[BC]-[CD]>.'"
        p = Pattern.Prosite(pattern = "[AB]-[BC]-[CD]>.")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC][CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])([CD])$'")
        self.assertEqual(p.tostring(), "[AB]-[BC]-[CD]>.")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern BBD
        self.assertNotEqual(p.re.search('BBD'), None)
        m = p.search(Seq.Seq('BBD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QABC
        self.assertNotEqual(p.re.search('QABC'), None)
        m = p.search(Seq.Seq('QABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ABCQ
        self.assertEqual(p.re.search("ABCQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABCQ")), None)
        # Bad pattern BBDQ
        self.assertEqual(p.re.search("BBDQ"), None)
        self.assertEqual(p.search(Seq.Seq("BBDQ")), None)
        # Bad pattern QQQ
        self.assertEqual(p.re.search("QQQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQQ")), None)

    def test_pattern35(self):
        "Testing Prosite pattern '[AB]-[BC]-[CDEFGHIKLMN]>.'"
        p = Pattern.Prosite(pattern = "[AB]-[BC]-[CDEFGHIKLMN]>.")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC][CDEFGHIKLMN]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])([CDEFGHIKLMN])$'")
        self.assertEqual(p.tostring(), "[AB]-[BC]-[CDEFGHIKLMN]>.")
        # Good pattern ABN
        self.assertNotEqual(p.re.search('ABN'), None)
        m = p.search(Seq.Seq('ABN'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Good pattern QQABN
        self.assertNotEqual(p.re.search('QQABN'), None)
        m = p.search(Seq.Seq('QQABN'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 2)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Bad pattern ABB
        self.assertEqual(p.re.search("ABB"), None)
        self.assertEqual(p.search(Seq.Seq("ABB")), None)
        # Bad pattern BBB
        self.assertEqual(p.re.search("BBB"), None)
        self.assertEqual(p.search(Seq.Seq("BBB")), None)
        # Bad pattern ACCD
        self.assertEqual(p.re.search("ACCD"), None)
        self.assertEqual(p.search(Seq.Seq("ACCD")), None)

    def test_pattern36(self):
        "Testing Prosite pattern '{AB}-[BC]-[CD]>.'"
        p = Pattern.Prosite(pattern = "{AB}-[BC]-[CD]>.")
        self.assertEqual(repr(p.re.pattern), "'[^AB][BC][CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^AB])([BC])([CD])$'")
        self.assertEqual(p.tostring(), "{AB}-[BC]-[CD]>.")
        # Good pattern CCC
        self.assertNotEqual(p.re.search('CCC'), None)
        m = p.search(Seq.Seq('CCC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern CCD
        self.assertNotEqual(p.re.search('CCD'), None)
        m = p.search(Seq.Seq('CCD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern CCDQ
        self.assertEqual(p.re.search("CCDQ"), None)
        self.assertEqual(p.search(Seq.Seq("CCDQ")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)

    def test_pattern37(self):
        "Testing Prosite pattern '[AB]-{BC}-[CD]>.'"
        p = Pattern.Prosite(pattern = "[AB]-{BC}-[CD]>.")
        self.assertEqual(repr(p.re.pattern), "'[AB][^BC][CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([^BC])([CD])$'")
        self.assertEqual(p.tostring(), "[AB]-{BC}-[CD]>.")
        # Good pattern AAC
        self.assertNotEqual(p.re.search('AAC'), None)
        m = p.search(Seq.Seq('AAC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QBQD
        self.assertNotEqual(p.re.search('QBQD'), None)
        m = p.search(Seq.Seq('QBQD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern AACQ
        self.assertEqual(p.re.search("AACQ"), None)
        self.assertEqual(p.search(Seq.Seq("AACQ")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)

    def test_pattern38(self):
        "Testing Prosite pattern '[AB]-[BC]-{CD}>.'"
        p = Pattern.Prosite(pattern = "[AB]-[BC]-{CD}>.")
        self.assertEqual(repr(p.re.pattern), "'[AB][BC][^CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB])([BC])([^CD])$'")
        self.assertEqual(p.tostring(), "[AB]-[BC]-{CD}>.")
        # Good pattern ABB
        self.assertNotEqual(p.re.search('ABB'), None)
        m = p.search(Seq.Seq('ABB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QABE
        self.assertNotEqual(p.re.search('QABE'), None)
        m = p.search(Seq.Seq('QABE'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern ACBQ
        self.assertEqual(p.re.search("ACBQ"), None)
        self.assertEqual(p.search(Seq.Seq("ACBQ")), None)
        # Bad pattern QABEQ
        self.assertEqual(p.re.search("QABEQ"), None)
        self.assertEqual(p.search(Seq.Seq("QABEQ")), None)
        # Bad pattern P
        self.assertEqual(p.re.search("P"), None)
        self.assertEqual(p.search(Seq.Seq("P")), None)

    def test_pattern39(self):
        "Testing Prosite pattern '{AB}-[BC]-{CD}>.'"
        p = Pattern.Prosite(pattern = "{AB}-[BC]-{CD}>.")
        self.assertEqual(repr(p.re.pattern), "'[^AB][BC][^CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^AB])([BC])([^CD])$'")
        self.assertEqual(p.tostring(), "{AB}-[BC]-{CD}>.")
        # Good pattern CCA
        self.assertNotEqual(p.re.search('CCA'), None)
        m = p.search(Seq.Seq('CCA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern AAACCA
        self.assertNotEqual(p.re.search('AAACCA'), None)
        m = p.search(Seq.Seq('AAACCA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 3)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 3)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern CCAC
        self.assertEqual(p.re.search("CCAC"), None)
        self.assertEqual(p.search(Seq.Seq("CCAC")), None)
        # Bad pattern AAACCAD
        self.assertEqual(p.re.search("AAACCAD"), None)
        self.assertEqual(p.search(Seq.Seq("AAACCAD")), None)

    def test_pattern40(self):
        "Testing Prosite pattern '<A>.'"
        p = Pattern.Prosite(pattern = "<A>.")
        self.assertEqual(repr(p.re.pattern), "'^A$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^(A)$'")
        self.assertEqual(p.tostring(), "<A>.")
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)
        # Bad pattern AA
        self.assertEqual(p.re.search("AA"), None)
        self.assertEqual(p.search(Seq.Seq("AA")), None)
        # Bad pattern BAB
        self.assertEqual(p.re.search("BAB"), None)
        self.assertEqual(p.search(Seq.Seq("BAB")), None)
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)

    def test_pattern41(self):
        "Testing Prosite pattern '<{A}>.'"
        p = Pattern.Prosite(pattern = "<{A}>.")
        self.assertEqual(repr(p.re.pattern), "'^[^A]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^A])$'")
        self.assertEqual(p.tostring(), "<{A}>.")
        # Good pattern B
        self.assertNotEqual(p.re.search('B'), None)
        m = p.search(Seq.Seq('B'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern C
        self.assertNotEqual(p.re.search('C'), None)
        m = p.search(Seq.Seq('C'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern BC
        self.assertEqual(p.re.search("BC"), None)
        self.assertEqual(p.search(Seq.Seq("BC")), None)
        # Bad pattern BC
        self.assertEqual(p.re.search("BC"), None)
        self.assertEqual(p.search(Seq.Seq("BC")), None)
        # Bad pattern AA
        self.assertEqual(p.re.search("AA"), None)
        self.assertEqual(p.search(Seq.Seq("AA")), None)

    def test_pattern42(self):
        "Testing Prosite pattern '<[AB]>.'"
        p = Pattern.Prosite(pattern = "<[AB]>.")
        self.assertEqual(repr(p.re.pattern), "'^[AB]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])$'")
        self.assertEqual(p.tostring(), "<[AB]>.")
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern B
        self.assertNotEqual(p.re.search('B'), None)
        m = p.search(Seq.Seq('B'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Bad pattern C
        self.assertEqual(p.re.search("C"), None)
        self.assertEqual(p.search(Seq.Seq("C")), None)
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)
        # Bad pattern AA
        self.assertEqual(p.re.search("AA"), None)
        self.assertEqual(p.search(Seq.Seq("AA")), None)
        # Bad pattern BB
        self.assertEqual(p.re.search("BB"), None)
        self.assertEqual(p.search(Seq.Seq("BB")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)

    def test_pattern43(self):
        "Testing Prosite pattern '<{AB}>.'"
        p = Pattern.Prosite(pattern = "<{AB}>.")
        self.assertEqual(repr(p.re.pattern), "'^[^AB]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^AB])$'")
        self.assertEqual(p.tostring(), "<{AB}>.")
        # Good pattern C
        self.assertNotEqual(p.re.search('C'), None)
        m = p.search(Seq.Seq('C'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern Q
        self.assertNotEqual(p.re.search('Q'), None)
        m = p.search(Seq.Seq('Q'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        # Bad pattern A
        self.assertEqual(p.re.search("A"), None)
        self.assertEqual(p.search(Seq.Seq("A")), None)
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)
        # Bad pattern AC
        self.assertEqual(p.re.search("AC"), None)
        self.assertEqual(p.search(Seq.Seq("AC")), None)
        # Bad pattern BA
        self.assertEqual(p.re.search("BA"), None)
        self.assertEqual(p.search(Seq.Seq("BA")), None)

    def test_pattern44(self):
        "Testing Prosite pattern '<A-B>.'"
        p = Pattern.Prosite(pattern = "<A-B>.")
        self.assertEqual(repr(p.re.pattern), "'^AB$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^(A)(B)$'")
        self.assertEqual(p.tostring(), "<A-B>.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Bad pattern AAB
        self.assertEqual(p.re.search("AAB"), None)
        self.assertEqual(p.search(Seq.Seq("AAB")), None)
        # Bad pattern ABB
        self.assertEqual(p.re.search("ABB"), None)
        self.assertEqual(p.search(Seq.Seq("ABB")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)
        # Bad pattern QABQ
        self.assertEqual(p.re.search("QABQ"), None)
        self.assertEqual(p.search(Seq.Seq("QABQ")), None)

    def test_pattern45(self):
        "Testing Prosite pattern '<[AB]-[BC]>.'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC]>.")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])$'")
        self.assertEqual(p.tostring(), "<[AB]-[BC]>.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BB
        self.assertNotEqual(p.re.search('BB'), None)
        m = p.search(Seq.Seq('BB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern BC
        self.assertNotEqual(p.re.search('BC'), None)
        m = p.search(Seq.Seq('BC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Bad pattern PAB
        self.assertEqual(p.re.search("PAB"), None)
        self.assertEqual(p.search(Seq.Seq("PAB")), None)
        # Bad pattern ABQ
        self.assertEqual(p.re.search("ABQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABQ")), None)
        # Bad pattern QABQ
        self.assertEqual(p.re.search("QABQ"), None)
        self.assertEqual(p.search(Seq.Seq("QABQ")), None)

    def test_pattern46(self):
        "Testing Prosite pattern '<[AB]-[BC]-[CD]>.'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC]-[CD]>.")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC][CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])([CD])$'")
        self.assertEqual(p.tostring(), "<[AB]-[BC]-[CD]>.")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern BBD
        self.assertNotEqual(p.re.search('BBD'), None)
        m = p.search(Seq.Seq('BBD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern QABC
        self.assertEqual(p.re.search("QABC"), None)
        self.assertEqual(p.search(Seq.Seq("QABC")), None)
        # Bad pattern ABCQ
        self.assertEqual(p.re.search("ABCQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABCQ")), None)
        # Bad pattern QABCQ
        self.assertEqual(p.re.search("QABCQ"), None)
        self.assertEqual(p.search(Seq.Seq("QABCQ")), None)
        # Bad pattern QQ
        self.assertEqual(p.re.search("QQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQ")), None)

    def test_pattern47(self):
        "Testing Prosite pattern '<[AB]-[BC]-[CDEFGHIKLMN]>.'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC]-[CDEFGHIKLMN]>.")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC][CDEFGHIKLMN]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])([CDEFGHIKLMN])$'")
        self.assertEqual(p.tostring(), "<[AB]-[BC]-[CDEFGHIKLMN]>.")
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Good pattern ACI
        self.assertNotEqual(p.re.search('ACI'), None)
        m = p.search(Seq.Seq('ACI'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Good pattern BCN
        self.assertNotEqual(p.re.search('BCN'), None)
        m = p.search(Seq.Seq('BCN'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CDEFGHIKLMN') 
        # Bad pattern QABC
        self.assertEqual(p.re.search("QABC"), None)
        self.assertEqual(p.search(Seq.Seq("QABC")), None)
        # Bad pattern ABCQ
        self.assertEqual(p.re.search("ABCQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABCQ")), None)
        # Bad pattern QACI
        self.assertEqual(p.re.search("QACI"), None)
        self.assertEqual(p.search(Seq.Seq("QACI")), None)
        # Bad pattern ACIQ
        self.assertEqual(p.re.search("ACIQ"), None)
        self.assertEqual(p.search(Seq.Seq("ACIQ")), None)
        # Bad pattern QACIQ
        self.assertEqual(p.re.search("QACIQ"), None)
        self.assertEqual(p.search(Seq.Seq("QACIQ")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)
        # Bad pattern QQQ
        self.assertEqual(p.re.search("QQQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQQ")), None)

    def test_pattern48(self):
        "Testing Prosite pattern '<{AB}-[BC]-[CD]>.'"
        p = Pattern.Prosite(pattern = "<{AB}-[BC]-[CD]>.")
        self.assertEqual(repr(p.re.pattern), "'^[^AB][BC][CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^AB])([BC])([CD])$'")
        self.assertEqual(p.tostring(), "<{AB}-[BC]-[CD]>.")
        # Good pattern CCC
        self.assertNotEqual(p.re.search('CCC'), None)
        m = p.search(Seq.Seq('CCC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern CBC
        self.assertNotEqual(p.re.search('CBC'), None)
        m = p.search(Seq.Seq('CBC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QCD
        self.assertNotEqual(p.re.search('QCD'), None)
        m = p.search(Seq.Seq('QCD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern CCCC
        self.assertEqual(p.re.search("CCCC"), None)
        self.assertEqual(p.search(Seq.Seq("CCCC")), None)
        # Bad pattern QCDC
        self.assertEqual(p.re.search("QCDC"), None)
        self.assertEqual(p.search(Seq.Seq("QCDC")), None)
        # Bad pattern QCDQ
        self.assertEqual(p.re.search("QCDQ"), None)
        self.assertEqual(p.search(Seq.Seq("QCDQ")), None)

    def test_pattern49(self):
        "Testing Prosite pattern '<[AB]-{BC}-[CD]>.'"
        p = Pattern.Prosite(pattern = "<[AB]-{BC}-[CD]>.")
        self.assertEqual(repr(p.re.pattern), "'^[AB][^BC][CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([^BC])([CD])$'")
        self.assertEqual(p.tostring(), "<[AB]-{BC}-[CD]>.")
        # Good pattern AAD
        self.assertNotEqual(p.re.search('AAD'), None)
        m = p.search(Seq.Seq('AAD'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern AQC
        self.assertNotEqual(p.re.search('AQC'), None)
        m = p.search(Seq.Seq('AQC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern QAAD
        self.assertEqual(p.re.search("QAAD"), None)
        self.assertEqual(p.search(Seq.Seq("QAAD")), None)
        # Bad pattern AADQ
        self.assertEqual(p.re.search("AADQ"), None)
        self.assertEqual(p.search(Seq.Seq("AADQ")), None)
        # Bad pattern QAAD
        self.assertEqual(p.re.search("QAAD"), None)
        self.assertEqual(p.search(Seq.Seq("QAAD")), None)
        # Bad pattern QQQ
        self.assertEqual(p.re.search("QQQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQQ")), None)

    def test_pattern50(self):
        "Testing Prosite pattern '<[AB]-[BC]-{CD}>.'"
        p = Pattern.Prosite(pattern = "<[AB]-[BC]-{CD}>.")
        self.assertEqual(repr(p.re.pattern), "'^[AB][BC][^CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([AB])([BC])([^CD])$'")
        self.assertEqual(p.tostring(), "<[AB]-[BC]-{CD}>.")
        # Good pattern ABB
        self.assertNotEqual(p.re.search('ABB'), None)
        m = p.search(Seq.Seq('ABB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern BBB
        self.assertNotEqual(p.re.search('BBB'), None)
        m = p.search(Seq.Seq('BBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern QABB
        self.assertEqual(p.re.search("QABB"), None)
        self.assertEqual(p.search(Seq.Seq("QABB")), None)
        # Bad pattern BABQ
        self.assertEqual(p.re.search("BABQ"), None)
        self.assertEqual(p.search(Seq.Seq("BABQ")), None)

    def test_pattern51(self):
        "Testing Prosite pattern '<{AB}-[BC]-{CD}>.'"
        p = Pattern.Prosite(pattern = "<{AB}-[BC]-{CD}>.")
        self.assertEqual(repr(p.re.pattern), "'^[^AB][BC][^CD]$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^([^AB])([BC])([^CD])$'")
        self.assertEqual(p.tostring(), "<{AB}-[BC]-{CD}>.")
        # Good pattern CCA
        self.assertNotEqual(p.re.search('CCA'), None)
        m = p.search(Seq.Seq('CCA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern FCF
        self.assertNotEqual(p.re.search('FCF'), None)
        m = p.search(Seq.Seq('FCF'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern QQQQ
        self.assertEqual(p.re.search("QQQQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQQQ")), None)
        # Bad pattern QCCAQ
        self.assertEqual(p.re.search("QCCAQ"), None)
        self.assertEqual(p.search(Seq.Seq("QCCAQ")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)

    def test_pattern52(self):
        "Testing Prosite pattern 'A-[BC>].'"
        p = Pattern.Prosite(pattern = "A-[BC>].")
        self.assertEqual(repr(p.re.pattern), "'A(?:[BC]|$)'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)((?:[BC]|$))'")
        self.assertEqual(p.tostring(), "A-[BC>].")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern AQ
        self.assertEqual(p.re.search("AQ"), None)
        self.assertEqual(p.search(Seq.Seq("AQ")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)

    def test_pattern53(self):
        "Testing Prosite pattern '<A-[BC>].'"
        p = Pattern.Prosite(pattern = "<A-[BC>].")
        self.assertEqual(repr(p.re.pattern), "'^A(?:[BC]|$)'")
        self.assertEqual(repr(p.grouped_re.pattern), "'^(A)((?:[BC]|$))'")
        self.assertEqual(p.tostring(), "<A-[BC>].")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern AA
        self.assertEqual(p.re.search("AA"), None)
        self.assertEqual(p.search(Seq.Seq("AA")), None)
        # Bad pattern AQ
        self.assertEqual(p.re.search("AQ"), None)
        self.assertEqual(p.search(Seq.Seq("AQ")), None)
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)

    def test_pattern54(self):
        "Testing Prosite pattern 'A-[BC>]>.'"
        p = Pattern.Prosite(pattern = "A-[BC>]>.")
        self.assertEqual(repr(p.re.pattern), "'A(?:[BC]|$)$'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)((?:[BC]|$))$'")
        self.assertEqual(p.tostring(), "A-[BC>]>.")
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'BC') 
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern ABQ
        self.assertEqual(p.re.search("ABQ"), None)
        self.assertEqual(p.search(Seq.Seq("ABQ")), None)

    def test_pattern55(self):
        "Testing Prosite pattern '[<AB]-C.'"
        p = Pattern.Prosite(pattern = "[<AB]-C.")
        self.assertEqual(repr(p.re.pattern), "'(?:^|[AB])C'")
        self.assertEqual(repr(p.grouped_re.pattern), "'((?:^|[AB]))(C)'")
        self.assertEqual(p.tostring(), "[<AB]-C.")
        # Good pattern C
        self.assertNotEqual(p.re.search('C'), None)
        m = p.search(Seq.Seq('C'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'C') 
        # Good pattern AC
        self.assertNotEqual(p.re.search('AC'), None)
        m = p.search(Seq.Seq('AC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'C') 
        # Good pattern BC
        self.assertNotEqual(p.re.search('BC'), None)
        m = p.search(Seq.Seq('BC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'C') 
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'C') 
        # Bad pattern QCC
        self.assertEqual(p.re.search("QCC"), None)
        self.assertEqual(p.search(Seq.Seq("QCC")), None)
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)

    def test_pattern56(self):
        "Testing Prosite pattern '[<AB]-[CD>].'"
        p = Pattern.Prosite(pattern = "[<AB]-[CD>].")
        self.assertEqual(repr(p.re.pattern), "'(?:^|[AB])(?:[CD]|$)'")
        self.assertEqual(repr(p.grouped_re.pattern), "'((?:^|[AB]))((?:[CD]|$))'")
        self.assertEqual(p.tostring(), "[<AB]-[CD>].")
        # Good pattern 
        self.assertNotEqual(p.re.search(''), None)
        m = p.search(Seq.Seq(''))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 0)
        self.assertEqual(m.start(), 0)
        # Good pattern QA
        self.assertNotEqual(p.re.search('QA'), None)
        m = p.search(Seq.Seq('QA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern AB
        self.assertNotEqual(p.re.search('AB'), None)
        m = p.search(Seq.Seq('AB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern QADQ
        self.assertNotEqual(p.re.search('QADQ'), None)
        m = p.search(Seq.Seq('QADQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QAB
        self.assertNotEqual(p.re.search('QAB'), None)
        m = p.search(Seq.Seq('QAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 2)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern ADQ
        self.assertNotEqual(p.re.search('ADQ'), None)
        m = p.search(Seq.Seq('ADQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern QQB
        self.assertNotEqual(p.re.search('QQB'), None)
        m = p.search(Seq.Seq('QQB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 2)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern QQBDQ
        self.assertNotEqual(p.re.search('QQBDQ'), None)
        m = p.search(Seq.Seq('QQBDQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 2)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 2)
        self.assertEqual(m.start(), 2)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern C
        self.assertNotEqual(p.re.search('C'), None)
        m = p.search(Seq.Seq('C'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Good pattern B
        self.assertNotEqual(p.re.search('B'), None)
        m = p.search(Seq.Seq('B'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern CQ
        self.assertNotEqual(p.re.search('CQ'), None)
        m = p.search(Seq.Seq('CQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 1)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'CD') 
        # Bad pattern Q
        self.assertEqual(p.re.search("Q"), None)
        self.assertEqual(p.search(Seq.Seq("Q")), None)
        # Bad pattern QQ
        self.assertEqual(p.re.search("QQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQ")), None)
        # Bad pattern QAQ
        self.assertEqual(p.re.search("QAQ"), None)
        self.assertEqual(p.search(Seq.Seq("QAQ")), None)
        # Bad pattern QC
        self.assertEqual(p.re.search("QC"), None)
        self.assertEqual(p.search(Seq.Seq("QC")), None)
        # Bad pattern AQ
        self.assertEqual(p.re.search("AQ"), None)
        self.assertEqual(p.search(Seq.Seq("AQ")), None)
        # Bad pattern QCQ
        self.assertEqual(p.re.search("QCQ"), None)
        self.assertEqual(p.search(Seq.Seq("QCQ")), None)

    def test_pattern57(self):
        "Testing Prosite pattern '[<ABC>].'"
        p = Pattern.Prosite(pattern = "[<ABC>].")
        self.assertEqual(repr(p.re.pattern), "'(?:^|[ABC$])'")
        self.assertEqual(repr(p.grouped_re.pattern), "'((?:^|[ABC$]))'")
        self.assertEqual(p.tostring(), "[<ABC>].")
        # Good pattern 
        self.assertNotEqual(p.re.search(''), None)
        m = p.search(Seq.Seq(''))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 0)
        self.assertEqual(m.start(), 0)
        # Good pattern A
        self.assertNotEqual(p.re.search('A'), None)
        m = p.search(Seq.Seq('A'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 0)
        self.assertEqual(m.start(), 0)
        # Good pattern AQ
        self.assertNotEqual(p.re.search('AQ'), None)
        m = p.search(Seq.Seq('AQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 0)
        self.assertEqual(m.start(), 0)
        # Good pattern QA
        self.assertNotEqual(p.re.search('QA'), None)
        m = p.search(Seq.Seq('QA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 0)
        self.assertEqual(m.start(), 0)
        # Good pattern QAQ
        self.assertNotEqual(p.re.search('QAQ'), None)
        m = p.search(Seq.Seq('QAQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 0)
        self.assertEqual(m.start(), 0)

    def test_pattern58(self):
        "Testing Prosite pattern 'A(3).'"
        p = Pattern.Prosite(pattern = "A(3).")
        self.assertEqual(repr(p.re.pattern), "'A{3}'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A{3})'")
        self.assertEqual(p.tostring(), "A(3).")
        # Good pattern AAA
        self.assertNotEqual(p.re.search('AAA'), None)
        m = p.search(Seq.Seq('AAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern AAAB
        self.assertNotEqual(p.re.search('AAAB'), None)
        m = p.search(Seq.Seq('AAAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern AAAAA
        self.assertNotEqual(p.re.search('AAAAA'), None)
        m = p.search(Seq.Seq('AAAAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern BAAA
        self.assertNotEqual(p.re.search('BAAA'), None)
        m = p.search(Seq.Seq('BAAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 1)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 1)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern AA
        self.assertEqual(p.re.search("AA"), None)
        self.assertEqual(p.search(Seq.Seq("AA")), None)
        # Bad pattern AABA
        self.assertEqual(p.re.search("AABA"), None)
        self.assertEqual(p.search(Seq.Seq("AABA")), None)
        # Bad pattern BBB
        self.assertEqual(p.re.search("BBB"), None)
        self.assertEqual(p.search(Seq.Seq("BBB")), None)
        # Bad pattern (empty string)
        self.assertEqual(p.re.search(""), None)
        self.assertEqual(p.search(Seq.Seq("")), None)

    def test_pattern59(self):
        "Testing Prosite pattern 'x(3).'"
        p = Pattern.Prosite(pattern = "x(3).")
        self.assertEqual(repr(p.re.pattern), "'.{3}'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(.{3})'")
        self.assertEqual(p.tostring(), "x(3).")
        # Good pattern AAA
        self.assertNotEqual(p.re.search('AAA'), None)
        m = p.search(Seq.Seq('AAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        # Good pattern BBB
        self.assertNotEqual(p.re.search('BBB'), None)
        m = p.search(Seq.Seq('BBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        # Good pattern ABC
        self.assertNotEqual(p.re.search('ABC'), None)
        m = p.search(Seq.Seq('ABC'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        # Bad pattern AA
        self.assertEqual(p.re.search("AA"), None)
        self.assertEqual(p.search(Seq.Seq("AA")), None)
        # Bad pattern B
        self.assertEqual(p.re.search("B"), None)
        self.assertEqual(p.search(Seq.Seq("B")), None)
        # Bad pattern (empty string)
        self.assertEqual(p.re.search(""), None)
        self.assertEqual(p.search(Seq.Seq("")), None)

    def test_pattern60(self):
        "Testing Prosite pattern 'A(3)-B(3).'"
        p = Pattern.Prosite(pattern = "A(3)-B(3).")
        self.assertEqual(repr(p.re.pattern), "'A{3}B{3}'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A{3})(B{3})'")
        self.assertEqual(p.tostring(), "A(3)-B(3).")
        # Good pattern AAABBB
        self.assertNotEqual(p.re.search('AAABBB'), None)
        m = p.search(Seq.Seq('AAABBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 6)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[4]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[5]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern BBBAAABBB
        self.assertNotEqual(p.re.search('BBBAAABBB'), None)
        m = p.search(Seq.Seq('BBBAAABBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 3)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 6)
        self.assertEqual(m.start(), 3)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[4]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[5]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Bad pattern ABABABAB
        self.assertEqual(p.re.search("ABABABAB"), None)
        self.assertEqual(p.search(Seq.Seq("ABABABAB")), None)
        # Bad pattern AABBB
        self.assertEqual(p.re.search("AABBB"), None)
        self.assertEqual(p.search(Seq.Seq("AABBB")), None)

    def test_pattern61(self):
        "Testing Prosite pattern 'A(2,3)-B(1,3).'"
        p = Pattern.Prosite(pattern = "A(2,3)-B(1,3).")
        self.assertEqual(repr(p.re.pattern), "'A{2,3}B{1,3}'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A{2,3})(B{1,3})'")
        self.assertEqual(p.tostring(), "A(2,3)-B(1,3).")
        # Good pattern AABB
        self.assertNotEqual(p.re.search('AABB'), None)
        m = p.search(Seq.Seq('AABB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 4)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern AAB
        self.assertNotEqual(p.re.search('AAB'), None)
        m = p.search(Seq.Seq('AAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern AABBB
        self.assertNotEqual(p.re.search('AABBB'), None)
        m = p.search(Seq.Seq('AABBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 5)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        mpi = mapped_pattern[4]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern AAAB
        self.assertNotEqual(p.re.search('AAAB'), None)
        m = p.search(Seq.Seq('AAAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 4)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Bad pattern ABBB
        self.assertEqual(p.re.search("ABBB"), None)
        self.assertEqual(p.search(Seq.Seq("ABBB")), None)
        # Bad pattern QQQQQ
        self.assertEqual(p.re.search("QQQQQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQQQQ")), None)

    def test_pattern62(self):
        "Testing Prosite pattern 'A-x(1,5)-B.'"
        p = Pattern.Prosite(pattern = "A-x(1,5)-B.")
        self.assertEqual(repr(p.re.pattern), "'A.{1,5}B'")
        self.assertEqual(repr(p.grouped_re.pattern), "'(A)(.{1,5})(B)'")
        self.assertEqual(p.tostring(), "A-x(1,5)-B.")
        # Good pattern ABB
        self.assertNotEqual(p.re.search('ABB'), None)
        m = p.search(Seq.Seq('ABB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern ABBBBBB
        self.assertNotEqual(p.re.search('ABBBBBB'), None)
        m = p.search(Seq.Seq('ABBBBBB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 7)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[4]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[5]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[6]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern ACCCCCB
        self.assertNotEqual(p.re.search('ACCCCCB'), None)
        m = p.search(Seq.Seq('ACCCCCB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 7)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[4]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[5]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[6]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Good pattern ACB
        self.assertNotEqual(p.re.search('ACB'), None)
        m = p.search(Seq.Seq('ACB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'x') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'B') 
        # Bad pattern ACCCCCCB
        self.assertEqual(p.re.search("ACCCCCCB"), None)
        self.assertEqual(p.search(Seq.Seq("ACCCCCCB")), None)
        # Bad pattern AB
        self.assertEqual(p.re.search("AB"), None)
        self.assertEqual(p.search(Seq.Seq("AB")), None)

    def test_pattern63(self):
        "Testing Prosite pattern '[AB](3,4).'"
        p = Pattern.Prosite(pattern = "[AB](3,4).")
        self.assertEqual(repr(p.re.pattern), "'[AB]{3,4}'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([AB]{3,4})'")
        self.assertEqual(p.tostring(), "[AB](3,4).")
        # Good pattern AAA
        self.assertNotEqual(p.re.search('AAA'), None)
        m = p.search(Seq.Seq('AAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 3)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern AAAA
        self.assertNotEqual(p.re.search('AAAA'), None)
        m = p.search(Seq.Seq('AAAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 4)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern ABAB
        self.assertNotEqual(p.re.search('ABAB'), None)
        m = p.search(Seq.Seq('ABAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 4)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern BABA
        self.assertNotEqual(p.re.search('BABA'), None)
        m = p.search(Seq.Seq('BABA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 4)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Good pattern BBAB
        self.assertNotEqual(p.re.search('BBAB'), None)
        m = p.search(Seq.Seq('BBAB'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 4)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'AB') 
        # Bad pattern ABC
        self.assertEqual(p.re.search("ABC"), None)
        self.assertEqual(p.search(Seq.Seq("ABC")), None)
        # Bad pattern CBA
        self.assertEqual(p.re.search("CBA"), None)
        self.assertEqual(p.search(Seq.Seq("CBA")), None)
        # Bad pattern QQQQ
        self.assertEqual(p.re.search("QQQQ"), None)
        self.assertEqual(p.search(Seq.Seq("QQQQ")), None)

    def test_pattern64(self):
        "Testing Prosite pattern '{AB}(3,5)-A(2).'"
        p = Pattern.Prosite(pattern = "{AB}(3,5)-A(2).")
        self.assertEqual(repr(p.re.pattern), "'[^AB]{3,5}A{2}'")
        self.assertEqual(repr(p.grouped_re.pattern), "'([^AB]{3,5})(A{2})'")
        self.assertEqual(p.tostring(), "{AB}(3,5)-A(2).")
        # Good pattern CCCAA
        self.assertNotEqual(p.re.search('CCCAA'), None)
        m = p.search(Seq.Seq('CCCAA'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 0)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 5)
        self.assertEqual(m.start(), 0)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[4]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Good pattern QCQAQCQAAQ
        self.assertNotEqual(p.re.search('QCQAQCQAAQ'), None)
        m = p.search(Seq.Seq('QCQAQCQAAQ'))
        self.assertNotEqual(m, None)
        self.assertEqual(m.start(), 4)
        mapped_pattern = m.mapped_pattern()
        self.assertEqual(len(mapped_pattern), 5)
        self.assertEqual(m.start(), 4)
        mpi = mapped_pattern[0]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[1]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[2]
        self.assertEqual(mpi.ignore, 1) 
        self.assertEqual(mpi.letters, 'AB') 
        mpi = mapped_pattern[3]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        mpi = mapped_pattern[4]
        self.assertEqual(mpi.ignore, 0) 
        self.assertEqual(mpi.letters, 'A') 
        # Bad pattern ABAAA
        self.assertEqual(p.re.search("ABAAA"), None)
        self.assertEqual(p.search(Seq.Seq("ABAAA")), None)

    def test_bad_pattern(self):
        "Testing bad Prosite patterns"
        self.assertRaises(TypeError, Pattern.compile, "A")
        self.assertRaises(TypeError, Pattern.compile, "A<")
        self.assertRaises(TypeError, Pattern.compile, ">A")
        self.assertRaises(TypeError, Pattern.compile, "A<.")
        self.assertRaises(TypeError, Pattern.compile, ">A.")
        self.assertRaises(TypeError, Pattern.compile, "[AB]<.")
        self.assertRaises(TypeError, Pattern.compile, ">[AB].")
        self.assertRaises(TypeError, Pattern.compile, "AB.")
        self.assertRaises(TypeError, Pattern.compile, "A-B<.")
        self.assertRaises(TypeError, Pattern.compile, "A(B).")
        self.assertRaises(TypeError, Pattern.compile, "A.B")
        self.assertRaises(TypeError, Pattern.compile, "A.B.")
        self.assertRaises(TypeError, Pattern.compile, "A-B")
        self.assertRaises(TypeError, Pattern.compile, "[A-B]")
        self.assertRaises(TypeError, Pattern.compile, "{[A-B]}.")
        self.assertRaises(TypeError, Pattern.compile, "[{AB}].")
        self.assertRaises(TypeError, Pattern.compile, "(A).")
        self.assertRaises(TypeError, Pattern.compile, "[A>]-B.")
        self.assertRaises(TypeError, Pattern.compile, "A-[<B].")
        self.assertRaises(TypeError, Pattern.compile, "A(3,).")
        self.assertRaises(TypeError, Pattern.compile, "A(-1).")
        #"A(5,2).",  # too complicated to check via re
        self.assertRaises(TypeError, Pattern.compile, "A(2, 5).")
        self.assertRaises(TypeError, Pattern.compile, "A(B,C).")
        self.assertRaises(TypeError, Pattern.compile, "A-[<G].")
        self.assertRaises(TypeError, Pattern.compile, "[G>]-A.")
        self.assertRaises(TypeError, Pattern.compile, "1.")
        self.assertRaises(TypeError, Pattern.compile, "[1].")
        self.assertRaises(TypeError, Pattern.compile, "[A1B].")
        #"[AA].",   # too complicated to check via re
        self.assertRaises(TypeError, Pattern.compile, "{1}.")
        self.assertRaises(TypeError, Pattern.compile, "[<].")
        self.assertRaises(TypeError, Pattern.compile, "[>].")
        self.assertRaises(TypeError, Pattern.compile, "[A<B].")
        self.assertRaises(TypeError, Pattern.compile, "[A>B].")

    def test_verify_pattern(self):
        "Verify Prosite patterns found in Prosite files"
        filenames = ('ps00107.txt',
                     'ps00159.txt',
                     'ps00165.txt',
                     'ps00432.txt',
                     'ps00488.txt',
                     'ps00546.txt')

        for filename in filenames:
            filename = os.path.join('Prosite', filename)
            input = open(filename)
            pattern = ""
            for line in input:
                if line[:2] != "PA":
                    continue

                pattern += line[5:].strip()
                if pattern.endswith("."):
                    Pattern.prosite_to_re(pattern)
                    Pattern.prosite_to_grouped_re(pattern)
                    p = Pattern.compile(pattern)
                    self.assertEqual(str(p), pattern)
                    pattern = ""

    def test_conversion(self):
        "Test a conversion of a prosite pattern to a regular expression."

        pattern = '[LIV]-G-{P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-[GSTACLIVMFY]-x(5,18)-[LIVMFYWCSTAR]-[AIVP]-[LIVMFAGCKR]-K.'
        regular_expression = Pattern.prosite_to_re( pattern )
        self.assertEqual(regular_expression, '[LIV]G[^P]G[^P][FYWMGSTNH][SGA][^PW][LIVCAT][^PD].[GSTACLIVMFY].{5,18}[LIVMFYWCSTAR][AIVP][LIVMFAGCKR]K')

        pattern = '[IV]-x-D-S-[GAS]-[GASC]-[GAST]-[GA]-T.'
        regular_expression = Pattern.prosite_to_re(pattern)
        self.assertEqual(regular_expression, '[IV].DS[GAS][GASC][GAST][GA]T')

        pattern = 'G-[LIVM]-x(3)-E-[LIV]-T-[LF]-R.'
        regular_expression = Pattern.prosite_to_re(pattern)
        self.assertEqual(regular_expression, 'G[LIVM].{3}E[LIV]T[LF]R')

        pattern = '[DESH]-x(4,5)-[STVG]-x-[AS]-[FYI]-K-[DLIFSA]-[RVMF]-[GA]-[LIVMGA].'
        regular_expression = Pattern.prosite_to_re(pattern)
        self.assertEqual(regular_expression, '[DESH].{4,5}[STVG].[AS][FYI]K[DLIFSA][RVMF][GA][LIVMGA]')

        pattern = 'W-[IV]-[STA]-[RK]-x-[DE]-Y-[DNE]-[DE].'
        regular_expression = Pattern.prosite_to_re(pattern)
        self.assertEqual(regular_expression, 'W[IV][STA][RK].[DE]Y[DNE][DE]')

    def test_verify_pattern( self ):
        "Test verification of a pattern"

        # Good patterns
        pattern = 'W-[IV]-[STA]-[RK]-x-[DE]-Y-[DNE]-[DE].'
        self.assertTrue(Pattern.verify_pattern(pattern))

        pattern = '[LIV]-G-{P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-[GSTACLIVMFY]-x(5,18)-[LIVMFYWCSTAR]-[AIVP]-[LIVMFAGCKR]-K.'
        self.assertTrue(Pattern.verify_pattern(pattern))

        # Bad patterns
        pattern = 'W-[IV]-[STA*]-[RK]-x-[DE]-Y-[DNE]-[DE].'
        self.assertTrue(not Pattern.verify_pattern(pattern))

        pattern = 'W-[IV]-[STA-[RK]-x-[DE]-Y-[DNE]-[DE].'
        self.assertTrue(not Pattern.verify_pattern(pattern))

        pattern = '[LIV]-G-P}-G-{P}-[FYWMGSTNH]-[SGA]-{PW}-[LIVCAT]-{PD}-x-[GSTACLIVMFY]-x(5,18)-[LIVMFYWCSTAR]-[AIVP]-[LIVMFAGCKR]-K.'
        self.assertTrue(not Pattern.verify_pattern(pattern))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
