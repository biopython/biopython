import random
from Bio.Align import PairwiseAligner, PairwiseAligner_Old

def make_wsb_aligner(mode):
    algorithm = "Waterman-Smith-Beyer %s alignment algorithm" % mode
    breaks = [0, 3, 7]
    specificgaps = lambda x, y: (-2 - y) if x in breaks else (-0.5 - y)
    penalty  = -random.randint(0,4)/4.0
    nogaps = lambda x, y: penalty - y
    aligner = PairwiseAligner()
    aligner.match = 1
    aligner.mismatch = -1
    aligner.target_gap_score = nogaps
    aligner.query_gap_score = specificgaps
    aligner.mode = mode
    aligner_old = PairwiseAligner_Old()
    aligner_old.match = 1
    aligner_old.mismatch = -1
    aligner_old.target_gap_score = nogaps
    aligner_old.query_gap_score = specificgaps
    aligner_old.mode = mode
    return aligner, aligner_old

def make_wsb_aligner2(mode):
    def gap_score(i, n):
        if i==3:
            return -8
        if n == 1:
            return -1
        return -10
    algorithm = "Waterman-Smith-Beyer %s alignment algorithm" % mode
    aligner = PairwiseAligner()
    aligner.match = 1
    aligner.mismatch = -10
    aligner.target_gap_score = gap_score
    aligner.query_gap_score = gap_score
    aligner.mode = mode
    aligner_old = PairwiseAligner_Old()
    aligner_old.match = 1
    aligner_old.mismatch = -10
    aligner_old.target_gap_score = gap_score
    aligner_old.query_gap_score = gap_score
    aligner_old.mode = mode
    return aligner, aligner_old

def make_nwsw_aligner(mode):
    aligner = PairwiseAligner()
    aligner_old = PairwiseAligner_Old()
    aligner.mode = mode
    aligner.target_gap_score = -random.randint(1,16)/4.0
    aligner_old.target_gap_score = aligner.target_gap_score
    aligner.target_left_gap_score = -random.randint(1,16)/4.0
    aligner.target_right_gap_score = -random.randint(1,16)/4.0
    aligner.query_gap_score = -random.randint(1,16)/4.0
    aligner_old.query_gap_score = aligner.query_gap_score
    aligner.query_left_gap_score = -random.randint(1,16)/4.0
    aligner.query_right_gap_score = -random.randint(1,16)/4.0
    aligner_old.mode = mode
    aligner_old.target_left_gap_score = aligner.target_left_gap_score
    aligner_old.target_right_gap_score = aligner.target_right_gap_score
    aligner_old.query_left_gap_score = aligner.query_left_gap_score
    aligner_old.query_right_gap_score = aligner.query_right_gap_score
    return aligner, aligner_old

def make_gotoh_aligner(mode):
    aligner = PairwiseAligner()
    aligner.mode = mode
    algorithm = "Gotoh %s alignment algorithm" % mode
    while True:
        aligner.target_open_gap_score = -random.randint(1,16)/4.0
        aligner.target_extend_gap_score = -random.randint(1,16)/4.0
        aligner.target_left_open_gap_score = -random.randint(1,16)/4.0
        aligner.target_left_extend_gap_score = -random.randint(1,16)/4.0
        aligner.target_right_open_gap_score = -random.randint(1,16)/4.0
        aligner.target_right_extend_gap_score = -random.randint(1,16)/4.0
        aligner.query_open_gap_score = -random.randint(1,16)/4.0
        aligner.query_extend_gap_score = -random.randint(1,16)/4.0
        aligner.query_left_open_gap_score = -random.randint(1,16)/4.0
        aligner.query_left_extend_gap_score = -random.randint(1,16)/4.0
        aligner.query_right_open_gap_score = -random.randint(1,16)/4.0
        aligner.query_right_extend_gap_score = -random.randint(1,16)/4.0
        if aligner.algorithm == algorithm:
            break
    aligner_old = PairwiseAligner_Old()
    aligner_old.mode = mode
    aligner_old.target_internal_open_gap_score = aligner.target_internal_open_gap_score
    aligner_old.target_internal_extend_gap_score = aligner.target_internal_extend_gap_score
    aligner_old.target_left_open_gap_score = aligner.target_left_open_gap_score
    aligner_old.target_left_extend_gap_score = aligner.target_left_extend_gap_score
    aligner_old.target_right_open_gap_score = aligner.target_right_open_gap_score
    aligner_old.target_right_extend_gap_score = aligner.target_right_extend_gap_score
    aligner_old.query_internal_open_gap_score = aligner.query_internal_open_gap_score
    aligner_old.query_internal_extend_gap_score = aligner.query_internal_extend_gap_score
    aligner_old.query_left_open_gap_score = aligner.query_left_open_gap_score
    aligner_old.query_left_extend_gap_score = aligner.query_left_extend_gap_score
    aligner_old.query_right_open_gap_score = aligner.query_right_open_gap_score
    aligner_old.query_right_extend_gap_score = aligner.query_right_extend_gap_score
    return aligner, aligner_old

def test_gap_here_only_local_1():
    seq1 = "AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA"
    seq2 = "AABBBAAAACCCCAAAABBBAA"
    breaks = [0, 11, len(seq2)]
    # Very expensive to open a gap in seq1:
    nogaps = lambda x, y: -2000 - y
    # Very expensive to open a gap in seq2 unless it is in one of the allowed positions
    specificgaps = lambda x, y: (-2 - y) if x in breaks else (-2000 - y)
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match = 1
    aligner.mismatch = -1
    aligner.target_gap_score = nogaps
    aligner.query_gap_score = specificgaps
    assert aligner.algorithm == "Waterman-Smith-Beyer local alignment algorithm"
    score = aligner.score(seq1, seq2)
    assert score == 13
    alignments = aligner.align(seq1, seq2)
    assert len(alignments) == 2
    alignment = alignments[0]
    assert alignment.score == 13
    print alignment
    assert str(alignment) == """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
..|||||||||||||.....................
..AABBBAAAACCCCAAAABBBAA............
"""
    alignment = alignments[1]
    assert alignment.score == 13
    assert str(alignment) == """\
AAAABBBAAAACCCCCCCCCCCCCCAAAABBBAAAA
.....................|||||||||||||..
............AABBBAAAACCCCAAAABBBAA..
"""

if True:
    breaks = [0, 3, 7]
    specificgaps = lambda x, y: (-2 - y) if x in breaks else (-0.5 - y)
    penalty  = -random.randint(0,4)/4.0
    nogaps = lambda x, y: penalty - y
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match = 1
    aligner.mismatch = -1
    aligner.target_gap_score = nogaps
    aligner.query_gap_score = specificgaps
    s1 = "ACGT"
    s2 = "GACT"
    print aligner.algorithm
    counter = 0
    while True:
        counter += 1
        try:
            # alignments = aligner.score(s1, s2)
            alignments = aligner.align(s1, s2)
            alignment = next(alignments)
        except MemoryError:
            continue

raise Exception

if False:
    test_gap_here_only_local_1()

if False:
    for counter in range(1000):
        aligner, aligner_old = make_nwsw_aligner('global')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter

if False:
    for counter in range(1000):
        aligner, aligner_old = make_nwsw_aligner('local')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter

if False:
    for counter in range(1000):
        aligner, aligner_old = make_gotoh_aligner('global')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter

if False:
    for counter in range(1000):
        aligner, aligner_old = make_gotoh_aligner('local')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        n1 = 1
        n2 = 1
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter

if False:
    for counter in range(1000):
        aligner, aligner_old = make_wsb_aligner('global')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter

if False:
    for counter in range(1000):
        aligner, aligner_old = make_wsb_aligner2('global')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter

if False:
    for counter in range(1000):
        aligner, aligner_old = make_wsb_aligner('local')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter

if False:
    for counter in range(1000):
        aligner, aligner_old = make_wsb_aligner2('local')
        n1 = random.randint(16,24)
        n2 = random.randint(16,24)
        s1 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n1)])
        s2 = ''.join(['ACGT'[random.randint(0,3)] for i in range(n2)])
        score = aligner.score(s1, s2)
        score_old = aligner_old.score(s1, s2)
        alignments = aligner.align(s1, s2)
        alignments_old = aligner_old.align(s1, s2)
        n = 0
        for alignment, alignment_old in zip(alignments, alignments_old):
            assert alignment == alignment_old
            n += 1
        assert n == len(alignments)
        assert n == len(alignments_old)
        assert alignments.score == score
        assert alignments_old.score == score_old
        assert score == score_old
        if counter % 100 == 0: print counter
