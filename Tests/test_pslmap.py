import os
import random
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

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
    position += random.randint(20,40)
    blockSize = random.randint(20,40)
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
psl = format(alignment1, "psl")
handle = open("transcript.psl", "w")
handle.write(psl)
handle.close()

alignments = aligner.align(transcript, sequence)
alignment2 = alignments[0]
psl = format(alignment2, "psl")
handle = open("sequence.psl", "w")
handle.write(psl)
handle.close()

stdout = os.popen("pslMap sequence.psl transcript.psl stdout")
psl = stdout.read()

