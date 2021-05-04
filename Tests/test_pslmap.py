import os
import random
import unittest
from numpy import array

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, PairwiseAlignment

def map_alignment(alignment1, alignment2):
    path1 = array(alignment1.path)
    path2 = array(alignment2.path)
    path = []
    iEnd, qEnd = path2[0, :]
    for row2 in path2[1:, :]:
        iStart, qStart = iEnd, qEnd
        iEnd, qEnd = row2
        iBlockSize = iEnd - iStart
        qBlockSize = qEnd - qStart
        if qBlockSize == 0:
            continue
        if iBlockSize == 0:
            continue
        assert qBlockSize == iBlockSize
        index = path1[:, 1].searchsorted(iStart, side='right') - 1
        if index >= 0:
            offset = iStart - path1[index, 1]
            row = array([path1[index, 0] + offset, qStart])
        else:
            offset = path1[0, 1] - iStart
            row = array([0, offset])
            qBlockSize -= offset
        path.append(row)
        tStart = row[0]
        for row1 in path1[1:, :]:
            tBlockSize = row1[0] - tStart
            if tBlockSize > qBlockSize:
                row = path[-1] + qBlockSize
                path.append(row)
                break
            else:
                qBlockSize -= tBlockSize
                row =  path[-1] + tBlockSize
                path.append(row)
    target = alignment1.target
    query = alignment2.query
    alignment = PairwiseAlignment(target, query, path, None)
    return alignment

def map_check(alignment1, alignment2):
    psl1 = format(alignment1, "psl")
    handle = open("transcript.psl", "w")
    handle.write(psl1)
    handle.close()
    psl2 = format(alignment2, "psl")
    handle = open("sequence.psl", "w")
    handle.write(psl2)
    handle.close()
    stdout = os.popen("pslMap sequence.psl transcript.psl stdout")
    psl = stdout.read()
    os.remove("transcript.psl")
    os.remove("sequence.psl")
    return psl

class TestOneBlockTargetOneBlockQuery(unittest.TestCase):
    def test_internal(self):
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA")
        chromosome.id = "chromosome"
        transcript = Seq("GGGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("GGCCCCCGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(alignment1.path, ((12, 0), (31, 19)))
        self.assertEqual(str(alignment1), """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
            |||||||||||||||||||         
            GGGGGGGCCCCCGGGGGGA         
""")
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(alignment2.path, ((5, 0), (15, 10)))
        self.assertEqual(str(alignment2), """\
GGGGGGGCCCCCGGGGGGA
     ||||||||||    
     GGCCCCCGGG    
""")
        alignment = map_alignment(alignment1, alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0].tolist(), [17, 0])
        self.assertSequenceEqual(alignment.path[1].tolist(), [27, 10])
        self.assertEqual(str(alignment), """\
AAAAAAAAAAAAGGGGGGGCCCCCGGGGGGAAAAAAAAAA
                 ||||||||||             
                 GGCCCCCGGG             
""")
        psl = format(alignment, "psl")
        self.assertEqual(psl, """\
10	0	0	0	0	0	0	0	+	sequence	10	0	10	chromosome	40	17	27	1	10,	0,	17,
""")

    def test_left_overhang(self):
        chromosome = Seq("GGGCCCCCGGGGGGAAAAAAAAAA")
        chromosome.id = "chromosome"
        transcript = Seq("AGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("GGGGGCCCCCGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(str(alignment1), """\
   GGGCCCCCGGGGGGAAAAAAAAAA
   |||||||||||||||         
AGGGGGCCCCCGGGGGGA         
""")
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(str(alignment2), """\
AGGGGGCCCCCGGGGGGA
 |||||||||||||    
 GGGGGCCCCCGGG    
""")
        alignment = map_alignment(alignment1, alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0].tolist(), [0, 2])
        self.assertSequenceEqual(alignment.path[1].tolist(), [11, 13])
        self.assertEqual(str(alignment), """\
  GGGCCCCCGGGGGGAAAAAAAAAA
  |||||||||||             
GGGGGCCCCCGGG             
""")
        psl = format(alignment, "psl")
        self.assertEqual(psl, """\
11	0	0	0	0	0	0	0	+	sequence	13	2	13	chromosome	24	0	11	1	11,	2,	0,
""")

    def test_right_overhang(self):
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGG")
        chromosome.id = "chromosome"
        transcript = Seq("GGGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("GGCCCCCGGGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        self.assertEqual(len(alignments1), 1)
        alignment1 = alignments1[0]
        self.assertEqual(str(alignment1), """\
AAAAAAAAAAAAGGGGGGGCCCCCGGG    
            |||||||||||||||    
            GGGGGGGCCCCCGGGGGGA
""")
        alignments2 = aligner.align(transcript, sequence)
        self.assertEqual(len(alignments2), 1)
        alignment2 = alignments2[0]
        self.assertEqual(str(alignment2), """\
GGGGGGGCCCCCGGGGGGA
     ||||||||||||  
     GGCCCCCGGGGG  
""")
        alignment = map_alignment(alignment1, alignment2)
        self.assertEqual(len(alignment.path), 2)
        self.assertSequenceEqual(alignment.path[0].tolist(), [17, 0])
        self.assertSequenceEqual(alignment.path[1].tolist(), [27, 10])
        self.assertEqual(str(alignment), """\
AAAAAAAAAAAAGGGGGGGCCCCCGGG  
                 ||||||||||  
                 GGCCCCCGGGGG
""")
        psl = format(alignment, "psl")
        self.assertEqual(psl, """\
10	0	0	0	0	0	0	0	+	sequence	12	0	10	chromosome	27	17	27	1	10,	0,	17,
""")

aligner = PairwiseAligner()
aligner.internal_open_gap_score = -1
aligner.internal_extend_gap_score = -0.0
aligner.match_score = +1
aligner.mismatch_score = -1
aligner.mode = "local"

path1 = ((65, 0), (132, 67), (207, 67), (281, 141))
path2 = ((32, 0), (68, 36))
target1 = Seq(None, length=1000)
query1 = Seq(None, length=141)
target2 = Seq(None, length=141)
query2 = Seq(None, length=36)

target1.id = "chromosome"
query1.id = "transcript"
target2.id = "transcript"
query2.id = "sequence"

alignment1 = PairwiseAlignment(target1, query1, path1, None)
alignment2 = PairwiseAlignment(target2, query2, path2, None)
alignment = map_alignment(alignment1, alignment2)
psl = map_check(alignment1, alignment2)
assert format(alignment, "psl") == psl

if True:
    chromosome = "".join(['ACGT'[random.randint(0,3)] for i in range(1000)])
    nBlocks = random.randint(5, 10)
    nBlocks = 2
    transcript = ""
    position = 0
    for i in range(nBlocks):
        position += random.randint(60,80)
        blockSize = random.randint(60,80)
        transcript += chromosome[position:position+blockSize]
        position += blockSize
    nBlocks = 1
    sequence = ""
    position = 0
    for i in range(nBlocks):
        position += random.randint(20,40)
        blockSize = random.randint(20,40)
        sequence += transcript[position:position+blockSize]
        position += blockSize
    chromosome = Seq(chromosome)
    transcript = Seq(transcript)
    sequence = Seq(sequence)
    chromosome.id = "chromosome"
    transcript.id = "transcript"
    sequence.id = "sequence"
    alignments = aligner.align(chromosome, transcript)
    alignment1 = alignments[0]
    alignments = aligner.align(transcript, sequence)
    alignment2 = alignments[0]
    alignment = map_alignment(alignment1, alignment2)
    psl = map_check(alignment1, alignment2)
    assert format(alignment, "psl") == psl
    print("Randomized test OK")


for i in range(10000):
    test_random()

if False:
        chromosome = Seq("AAAAAAAAAAAAGGGGGGGCCCCCGGG")
        chromosome.id = "chromosome"
        transcript = Seq("GGGGGGGCCCCCGGGGGGA")
        transcript.id = "transcript"
        sequence = Seq("GGCCCCCGGGGG")
        sequence.id = "sequence"
        alignments1 = aligner.align(chromosome, transcript)
        assert len(alignments1) == 1
        alignment1 = alignments1[0]
        assert str(alignment1) == """\
AAAAAAAAAAAAGGGGGGGCCCCCGGG    
            |||||||||||||||    
            GGGGGGGCCCCCGGGGGGA
"""
        alignments2 = aligner.align(transcript, sequence)
        assert len(alignments2) == 1
        alignment2 = alignments2[0]
        assert str(alignment2) == """\
GGGGGGGCCCCCGGGGGGA
     ||||||||||||  
     GGCCCCCGGGGG  
"""
        alignment = map_alignment(alignment1, alignment2)
        assert str(alignment) == """\
AAAAAAAAAAAAGGGGGGGCCCCCGGG  
                 ||||||||||  
                 GGCCCCCGGGGG
"""
        psl = format(alignment, "psl")
        assert psl == """\
10	0	0	0	0	0	0	0	+	sequence	12	0	10	chromosome	27	17	27	1	10,	0,	17,
"""
        print("TEST OK")


# if __name__ == "__main__":
    # runner = unittest.TextTestRunner(verbosity=2)
    # unittest.main(testRunner=runner)
