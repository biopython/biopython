import os
import random
import unittest
from numpy import array

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, PairwiseAlignment


def map_alignment(alignment1, alignment2):
    if len(alignment1.query) != len(alignment2.target):
        raise ValueError("length of alignment1 query sequence (%d) != length of alignment2 target sequence (%d)" % (len(alignment1.query), len(alignment2.target)))
    target = alignment1.target
    query = alignment2.query
    path1 = alignment1.path
    path2 = alignment2.path
    n1 = len(alignment1.query)
    n2 = len(alignment2.query)
    if path1[0][1] < path1[-1][1]:  # mapped to forward strand
        strand1 = "+"
    else:  # mapped to reverse strand
        strand1 = "-"
    if path2[0][1] < path2[-1][1]:  # mapped to forward strand
        strand2 = "+"
    else:  # mapped to reverse strand
        strand2 = "-"
    path1 = array(path1)
    path2 = array(path2)
    if strand1 == "+":
        if strand2 == "-": # mapped to reverse strand
            path2[:, 1] = n2 - path2[:, 1]
    else:  # mapped to reverse strand
        path1[:, 1] = n1 - path1[:, 1]
        path2[:, 0] = n1 - path2[::-1, 0]
        if strand2 == "+":
            path2[:, 1] = n2 - path2[::-1, 1]
        else:  # mapped to reverse strand
            path2[:, 1] = path2[::-1, 1]
    path = []
    iMapBlk = 0
    previous2 = path2[0]
    for row in path2[1:]:
        if row[0] == previous2[0] or row[1] == previous2[1]:
            previous2 = row
            continue
        qStart = previous2[1]
        qEnd = row[1]
        tStart = previous2[0]
        tEnd = row[0]
        previous2 = row
        while True:
            if qStart >= qEnd or tStart >= tEnd:
                break
            blockCount = 0
            previous = path1[0]
            for row in path1[1:]:
                if previous[0] != row[0] and previous[1] != row[1]:
                    blockCount += 1
                previous = row
            for iBlk in range(iMapBlk, blockCount):
                i = 0
                previous = path1[0]
                for row in path1[1:]:
                    if previous[0] != row[0] and previous[1] != row[1]:
                        if i == iBlk:
                            break
                        i += 1
                    previous = row
                mqStart = previous[1]
                mqEnd = row[1]
                if tStart < mqStart:
                    iMapBlk = iBlk
                    if tEnd < mqStart:
                        size = tEnd - tStart
                    else:
                        size = mqStart - tStart
                    mappedBlk_qStart = qStart
                    mappedBlk_qEnd = qStart + size
                    mappedBlk_tStart = 0
                    mappedBlk_tEnd = 0
                    break
                elif tStart < mqEnd:
                    iMapBlk = iBlk
                    mtStart = previous[0]
                    off = tStart - mqStart
                    if tEnd > mqEnd:
                        size = mqEnd - tStart
                    else:
                        size = tEnd - tStart
                    mappedBlk_qStart = qStart
                    mappedBlk_qEnd = qStart + size
                    mappedBlk_tStart = mtStart + off
                    mappedBlk_tEnd = mtStart + off + size
                    break
            else:
                iMapBlk = iBlk
                mappedBlk_qStart = qStart
                mappedBlk_qEnd = qEnd
                mappedBlk_tStart = 0
                mappedBlk_tEnd = 0
            if mappedBlk_qEnd != 0 and mappedBlk_tEnd != 0:
                if path:
                    previous_tStart, previous_qStart = path[-1]
                    tGap = mappedBlk_tStart - previous_tStart
                    qGap = mappedBlk_qStart - previous_qStart
                    if tGap != 0 and qGap != 0:
                        # adding a gap both in target and in query;
                        # add gap to target first:
                        path.append([mappedBlk_tStart, previous_qStart])
                path.append([mappedBlk_tStart, mappedBlk_qStart])
                path.append([mappedBlk_tEnd, mappedBlk_qEnd])
            if mappedBlk_qEnd != 0:
                size = mappedBlk_qEnd - mappedBlk_qStart
            else:
                size = mappedBlk_tEnd - mappedBlk_tStart
            qStart += size
            tStart += size
    if strand1 != strand2:
        path = tuple((c1, n2 - c2) for (c1, c2) in path)
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

def test_random(nBlocks1=1, nBlocks2=1, strand1='+', strand2='+'):
    chromosome = "".join(['ACGT'[random.randint(0,3)] for i in range(1000)])
    nBlocks = nBlocks1
    transcript = ""
    position = 0
    for i in range(nBlocks):
        position += random.randint(60,80)
        blockSize = random.randint(60,80)
        transcript += chromosome[position:position+blockSize]
        position += blockSize
    nBlocks = nBlocks2
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
    if strand1 == '-':
        chromosome = chromosome.reverse_complement()
    if strand2 == '-':
        sequence = sequence.reverse_complement()
    chromosome.id = "chromosome"
    transcript.id = "transcript"
    sequence.id = "sequence"
    alignments1 = aligner.align(chromosome, transcript, strand=strand1)
    alignment1 = alignments1[0]
    alignments2 = aligner.align(transcript, sequence, strand=strand2)
    alignment2 = alignments2[0]
    alignment = map_alignment(alignment1, alignment2)
    psl_check = map_check(alignment1, alignment2)
    psl = format(alignment, "psl")
    assert psl == psl_check
    print("Randomized test %d, %d, %s, %s OK" % (nBlocks1, nBlocks2, strand1, strand2))

def test_random_sequences(strand1='+', strand2='+'):
    chromosome = "".join(['ACGT'[random.randint(0,3)] for i in range(1000)])
    transcript = "".join(['ACGT'[random.randint(0,3)] for i in range(300)])
    sequence = "".join(['ACGT'[random.randint(0,3)] for i in range(100)])
    chromosome = Seq(chromosome)
    transcript = Seq(transcript)
    sequence = Seq(sequence)
    chromosome.id = "chromosome"
    transcript.id = "transcript"
    sequence.id = "sequence"
    alignments = aligner.align(chromosome, transcript, strand=strand1)
    alignment1 = alignments[0]
    alignments = aligner.align(transcript, sequence, strand=strand2)
    alignment2 = alignments[0]
    psl_check = map_check(alignment1, alignment2)
    alignment = map_alignment(alignment1, alignment2)
    psl_check = psl_check.split()
    psl = format(alignment, "psl")
    psl = psl.split()
    assert psl[8:] == psl_check[8:]
    psl1 = format(alignment1, "psl")
    words = psl1.split()
    nBlocks1 = int(words[17])
    psl2 = format(alignment2, "psl")
    words = psl2.split()
    nBlocks2 = int(words[17])
    print("Randomized sequence test %d, %d, %s, %s OK" % (nBlocks1, nBlocks2, strand1, strand2))


for i in range(1000):
    nBlocks1 = random.randint(1,10)
    nBlocks2 = random.randint(1,10)
    test_random(nBlocks1, nBlocks2, '+', '+')
    test_random(nBlocks1, nBlocks2, '+', '-')
    test_random(nBlocks1, nBlocks2, '-', '+')
    test_random(nBlocks1, nBlocks2, '-', '-')
    test_random_sequences('+', '+')
    test_random_sequences('+', '-')
    test_random_sequences('-', '+')
    test_random_sequences('-', '-')

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
