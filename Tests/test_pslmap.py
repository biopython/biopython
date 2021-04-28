import os
import random
from numpy import array
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, PairwiseAlignment

aligner = PairwiseAligner()
aligner.internal_open_gap_score = -1
aligner.internal_extend_gap_score = -0.0
aligner.match_score = +1
aligner.mismatch_score = -1
aligner.mode = "local"

chromosome = "".join(['ACGT'[random.randint(0,3)] for i in range(1000)])
nBlocks = random.randint(5, 10)
nBlocks = 1
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
psl1 = format(alignment1, "psl")
handle = open("transcript.psl", "w")
handle.write(psl1)
handle.close()

alignments = aligner.align(transcript, sequence)
alignment2 = alignments[0]
psl2 = format(alignment2, "psl")
handle = open("sequence.psl", "w")
handle.write(psl2)
handle.close()

stdout = os.popen("pslMap sequence.psl transcript.psl stdout")
psl = stdout.read()

path1 = array(alignment1.path)
path2 = array(alignment2.path)

path = []
iEnd, qEnd = path2[0, :]
for row in path2[1:, :]:
    iStart, qStart = iEnd, qEnd
    iEnd, qEnd = row
    iBlockSize = iEnd - iStart
    qBlockSize = qEnd - qStart
    if qBlockSize == 0:
        continue
    if iBlockSize == 0:
        continue
    assert qBlockSize == iBlockSize
    index = path1[:, 1].searchsorted(iStart, side='right') - 1
    offset = iStart - path1[index, 1]
    row = array([path1[index, 0] + offset, qStart])
    path.append(row)
    tStart = row[0]
    for row in path1[index+1:, :]:
        tBlockSize = row[0] - tStart
        if tBlockSize > qBlockSize:
            row = path[-1] + qBlockSize
            path.append(row)
            break
        else:
            qBlockSize -= tBlockSize
            row =  path[-1] + tBlockSize
            path.append(row)

alignment = PairwiseAlignment(chromosome, sequence, path, None)
assert format(alignment, "psl") == psl
print("OK")
