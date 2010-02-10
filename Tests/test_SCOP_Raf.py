
"""Unit test for Raf"""

import unittest

from Bio.SCOP import Raf




class RafTests(unittest.TestCase):
    rafLine = "101m_ 0.01 38 010301 111011    0  153    0 mm   1 vv   2 ll   3 ss   4 ee   5 gg   6 ee   7 ww   8 qq   9 ll  10 vv  11 ll  12 hh  13 vv  14 ww  15 aa  16 kk  17 vv  18 ee  19 aa  20 dd  21 vv  22 aa  23 gg  24 hh  25 gg  26 qq  27 dd  28 ii  29 ll  30 ii  31 rr  32 ll  33 ff  34 kk  35 ss  36 hh  37 pp  38 ee  39 tt  40 ll  41 ee  42 kk  43 ff  44 dd  45 rr  46 vv  47 kk  48 hh  49 ll  50 kk  51 tt  52 ee  53 aa  54 ee  55 mm  56 kk  57 aa  58 ss  59 ee  60 dd  61 ll  62 kk  63 kk  64 hh  65 gg  66 vv  67 tt  68 vv  69 ll  70 tt  71 aa  72 ll  73 gg  74 aa  75 ii  76 ll  77 kk  78 kk  79 kk  80 gg  81 hh  82 hh  83 ee  84 aa  85 ee  86 ll  87 kk  88 pp  89 ll  90 aa  91 qq  92 ss  93 hh  94 aa  95 tt  96 kk  97 hh  98 kk  99 ii 100 pp 101 ii 102 kk 103 yy 104 ll 105 ee 106 ff 107 ii 108 ss 109 ee 110 aa 111 ii 112 ii 113 hh 114 vv 115 ll 116 hh 117 ss 118 rr 119 hh 120 pp 121 gg 122 nn 123 ff 124 gg 125 aa 126 dd 127 aa 128 qq 129 gg 130 aa 131 mm 132 nn 133 kk 134 aa 135 ll 136 ee 137 ll 138 ff 139 rr 140 kk 141 dd 142 ii 143 aa 144 aa 145 kk 146 yy 147 kk 148 ee 149 ll 150 gg 151 yy 152 qq 153 gg"

    rafLine2 = "101mA 0.01 38 010301 111011    0  153    0 mm   1 vv   2 ll   3 ss   4 ee   5 gg   6Aee   7Aww   8Aqq"

    rafLine3 = "101mB 0.01 38 010301 111011    0  153   90 mm  91 vv  92 ll  939ss  94 ee  95 gg"


    def testParse(self):
        """Can we parse a RAF record?"""
        r = Raf.SeqMap(self.rafLine)

        self.assertEqual(r.pdbid, "101m")
        self.assertEqual(r.pdb_datestamp, "010301")
        self.assertEqual(r.flags, "111011")
      
        i = r.index("143")
        res = r.res[i]
        self.assertEqual(res.chainid, "_")
        self.assertEqual(res.resid, "143")
        self.assertEqual(res.seqres, "A")
        self.assertEqual(res.atom, "A")

        r = Raf.SeqMap(self.rafLine2)   
        res = r.res[r.index("6A", chainid="A")]
        self.assertEqual(res.resid, "6A")
        self.assertEqual(res.atom, "E")

    def testSeqMapAdd(self):
        r2 = Raf.SeqMap(self.rafLine2)
        r3 = Raf.SeqMap(self.rafLine3)

        l = len(r2.res) + len(r3.res)
        r2 += r3
        self.assertEqual(len(r2.res), l)


        r2.extend(r2)
        self.assertEqual(len(r2.res), l*2)

        r4 = r2 + r2
        self.assertEqual(len(r4.res), l*4)

        r4.append(Raf.Res())
        self.assertEqual(len(r4.res), (l*4)+1)
        

    def testSeqMapSlice(self):
        r = Raf.SeqMap(self.rafLine)
        r = r[ r.index("124"): r.index("135")+1]
        self.assertEqual(len(r.res), 12)


    def testSeqMapIndex(self):
        filename = ("./SCOP/raftest.txt")
        
        index = Raf.SeqMapIndex(filename)
        r = index.getSeqMap("103m")
        self.assertEqual(r.pdbid, "103m")
        self.assertEqual(len(r.res), 154)
        self.assertEqual(r.pdb_datestamp, "010301")
        self.assertEqual(r.flags, "111011")

        r = index.getSeqMap("103m 1-10")
        self.assertEqual(r.pdbid, "103m",)
        self.assertEqual(len(r.res), 10)
        self.assertEqual(r.pdb_datestamp, "010301")
        self.assertEqual(r.flags, "111011")

        r = index.getSeqMap("104l A:")
        self.assertEqual(r.pdbid, "104l")

        r = index.getSeqMap("104l A:112-113")
        self.assertEqual(r.pdbid, "104l")
        self.assertEqual(len(r.res), 2)

        r = index.getSeqMap("104l A:112-113,B:146-148")
        self.assertEqual(r.pdbid, "104l")
        self.assertEqual(len(r.res), 5)

        
if __name__=='__main__':
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
