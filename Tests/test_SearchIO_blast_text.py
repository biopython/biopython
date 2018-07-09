# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO BlastIO plain text parsers."""

import os
import unittest

from Bio.SearchIO import parse


# test case files are in the Blast directory
TEST_DIR = 'Blast'
FMT = 'blast-text'


def get_file(filename):
    """Returns the path of a test file."""
    return os.path.join(TEST_DIR, filename)


class BaseBlastCases(unittest.TestCase):

    def check_common_attrs(self, qresults):
        # check common attributes
        for qresult in qresults:
            for hit in qresult:
                self.assertEqual(qresult.id, hit.query_id)
                for hsp in hit:
                    self.assertEqual(hit.id, hsp.hit_id)
                    self.assertEqual(qresult.id, hsp.query_id)


class BlastnCases(BaseBlastCases):

    def test_text_2226_blastn_001(self):
        """Test parsing blastn output (text_2226_blastn_001.txt)"""

        blast_file = get_file('text_2226_blastn_001.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

    def test_text_2226_blastn_002(self):
        """Test parsing blastn output (text_2226_blastn_002.txt)"""

        blast_file = get_file('text_2226_blastn_002.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|356995852:1-490', qresult.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcriptionfactor 1 (Pou5f1), transcript variant 1, mRNA', qresult.description)
        self.assertEqual(490, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(8, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|NM_013633.3|', hit.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcription factor 1 (Pou5f1), transcript variant 1, mRNA', hit.description)
        self.assertEqual(1353, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(490, hsp.aln_span)
        self.assertEqual(0.0, hsp.evalue)
        self.assertEqual(905.0, hsp.bitscore)
        self.assertEqual(490.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(490, hsp.query_end)
        self.assertEqual(490, hsp.hit_end)
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.hit.seq)[:40])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|XR_141831.1|', hit.id)
        self.assertEqual('PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA ref|XR_105837.2| PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA ref|XR_141464.1| PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA ref|XR_141446.1| PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA', hit.description)
        self.assertEqual(570, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(490, hsp.aln_span)
        self.assertEqual(0.0, hsp.evalue)
        self.assertEqual(900.0, hsp.bitscore)
        self.assertEqual(487.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(489, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(80, hsp.hit_start)
        self.assertEqual(490, hsp.query_end)
        self.assertEqual(570, hsp.hit_end)
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.hit.seq)[:40])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.hit.seq)[-40:])

    def test_text_2226_blastn_003(self):
        """Test parsing blastn output (text_2226_blastn_003.txt)"""

        blast_file = get_file('text_2226_blastn_003.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 strand=+repeatMasking=none', qresult.description)
        self.assertEqual(66, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|XM_003267724.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys ATG14 autophagy related 14 homolog (S. cerevisiae) (ATG14), mRNA', hit.description)
        self.assertEqual(4771, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(62, hsp.aln_span)
        self.assertEqual(3e-24, hsp.evalue)
        self.assertEqual(115.0, hsp.bitscore)
        self.assertEqual(62.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(2864, hsp.hit_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual(2926, hsp.hit_end)
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.hit.seq)[:40])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCaaaaaaaaaaaaaaaaaa', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|NM_001040441.1|', hit.id)
        self.assertEqual('Homo sapiens zinc finger and BTB domain containing 8A (ZBTB8A), mRNA', hit.description)
        self.assertEqual(7333, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(62, hsp.aln_span)
        self.assertEqual(3e-24, hsp.evalue)
        self.assertEqual(115.0, hsp.bitscore)
        self.assertEqual(62.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(3676, hsp.hit_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual(3738, hsp.hit_end)
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.hit.seq)[:40])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCaaaaaaaaaaaaaaaaaa', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', str(hsp.hit.seq)[-40:])
        # first qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(53, hsp.aln_span)
        self.assertEqual(3e-19, hsp.evalue)
        self.assertEqual(99.0, hsp.bitscore)
        self.assertEqual(53.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, hsp.ident_num)
        self.assertEqual(53, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(5, hsp.query_start)
        self.assertEqual(2823, hsp.hit_start)
        self.assertEqual(58, hsp.query_end)
        self.assertEqual(2876, hsp.hit_end)
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGT', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGT', str(hsp.hit.seq)[:40])
        self.assertEqual('GCCTGGGCAACAAGAGCGAAACTCCGTCTCaaaaaaaaaa', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('GCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', str(hsp.hit.seq)[-40:])

    def test_text_2226_blastn_004(self):
        """Test parsing blastn output (text_2226_blastn_004.txt)"""

        blast_file = get_file('text_2226_blastn_004.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(3, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|356995852:1-490', qresult.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcriptionfactor 1 (Pou5f1), transcript variant 1, mRNA', qresult.description)
        self.assertEqual(490, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(5, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|356995852|ref|NM_013633.3|', hit.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcription factor 1 (Pou5f1), transcript variant 1, mRNA', hit.description)
        self.assertEqual(1353, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(490, hsp.aln_span)
        self.assertEqual(0.0, hsp.evalue)
        self.assertEqual(905.0, hsp.bitscore)
        self.assertEqual(490.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(490, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(490, hsp.query_end)
        self.assertEqual(490, hsp.hit_end)
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.hit.seq)[:40])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.hit.seq)[-40:])
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|377833530|ref|XR_141831.1|', hit.id)
        self.assertEqual('PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA', hit.description)
        self.assertEqual(570, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(490, hsp.aln_span)
        self.assertEqual(0.0, hsp.evalue)
        self.assertEqual(900.0, hsp.bitscore)
        self.assertEqual(487.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(489, hsp.ident_num)
        self.assertEqual(490, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(80, hsp.hit_start)
        self.assertEqual(490, hsp.query_end)
        self.assertEqual(570, hsp.hit_end)
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GAGGTGAAACCGTCCCTAGGTGAGCCGTCTTTCCACCAGG', str(hsp.hit.seq)[:40])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AGTCCCAGGACATGAAAGCCCTGCAGAAGGAGCTAGAACA', str(hsp.hit.seq)[-40:])

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207307-1207372 5\'pad=0 3\'pad=0 strand=+repeatMasking=none', qresult.description)
        self.assertEqual(66, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('blastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(5, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|94721341|ref|NM_001040441.1|', hit.id)
        self.assertEqual('Homo sapiens zinc finger and BTB domain containing 8A (ZBTB8A), mRNA', hit.description)
        self.assertEqual(7333, hit.seq_len)
        self.assertEqual(2, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(62, hsp.aln_span)
        self.assertEqual(6e-29, hsp.evalue)
        self.assertEqual(115.0, hsp.bitscore)
        self.assertEqual(62.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(3676, hsp.hit_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual(3738, hsp.hit_end)
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.hit.seq)[:40])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCaaaaaaaaaaaaaaaaaa', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', str(hsp.hit.seq)[-40:])
        # third qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(53, hsp.aln_span)
        self.assertEqual(6e-24, hsp.evalue)
        self.assertEqual(99.0, hsp.bitscore)
        self.assertEqual(53.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(53, hsp.ident_num)
        self.assertEqual(53, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(5, hsp.query_start)
        self.assertEqual(2823, hsp.hit_start)
        self.assertEqual(58, hsp.query_end)
        self.assertEqual(2876, hsp.hit_end)
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGT', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('CCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCGT', str(hsp.hit.seq)[:40])
        self.assertEqual('GCCTGGGCAACAAGAGCGAAACTCCGTCTCaaaaaaaaaa', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('GCCTGGGCAACAAGAGCGAAACTCCGTCTCAAAAAAAAAA', str(hsp.hit.seq)[-40:])
        # third qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|332237160|ref|XM_003267724.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys ATG14 autophagy related 14 homolog (S. cerevisiae) (ATG14), mRNA', hit.description)
        self.assertEqual(4771, hit.seq_len)
        self.assertEqual(1, len(hit))
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(62, hsp.aln_span)
        self.assertEqual(6e-29, hsp.evalue)
        self.assertEqual(115.0, hsp.bitscore)
        self.assertEqual(62.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(62, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(4, hsp.query_start)
        self.assertEqual(2864, hsp.hit_start)
        self.assertEqual(66, hsp.query_end)
        self.assertEqual(2926, hsp.hit_end)
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.query.seq)[:40])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GCCATTGCACTCCAGCCTGGGCAACAAGAGCGAAACTCCG', str(hsp.hit.seq)[:40])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCaaaaaaaaaaaaaaaaaa', str(hsp.query.seq)[-40:])
        self.assertEqual('||||||||||||||||||||||||||||||||||||||||', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('AACAAGAGCGAAACTCCGTCTCAAAAAAAAAAAAAAAAAA', str(hsp.hit.seq)[-40:])


class BlastpCases(BaseBlastCases):

    def test_text_2226_blastp_001(self):
        """Test parsing blastp output (text_2226_blastp_001.txt)"""

        blast_file = get_file('text_2226_blastp_001.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual('NCBI Protein Reference Sequences', qresult.target)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

    def test_text_2226_blastp_002(self):
        """Test parsing blastp output (text_2226_blastp_002.txt)"""

        blast_file = get_file('text_2226_blastp_002.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillussubtilis subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual('NCBI Protein Reference Sequences', qresult.target)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|NP_391444.1|', hit.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168] ref|ZP_03593363.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. 168] ref|ZP_03597648.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. NCIB 3610] ref|ZP_03602051.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. JH642] ref|ZP_03606337.1| membrane bound lipoprotein [Bacillus subtilis subsp. subtilis str. SMY] ref|YP_004205398.1| unnamed protein product [Bacillus subtilis BSn5]', hit.description)
        self.assertEqual(102, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(102, hsp.aln_span)
        self.assertEqual(1e-66, hsp.evalue)
        self.assertEqual(205.0, hsp.bitscore)
        self.assertEqual(521.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(102, hsp.ident_num)
        self.assertEqual(102, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(102, hsp.query_end)
        self.assertEqual(102, hsp.hit_end)
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEG', str(hsp.query.seq)[:40])
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEG', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEG', str(hsp.hit.seq)[:40])
        self.assertEqual('DITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', str(hsp.query.seq)[-40:])
        self.assertEqual('DITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('DITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|YP_003922001.1|', hit.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]', hit.description)
        self.assertEqual(100, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(102, hsp.aln_span)
        self.assertEqual(1e-40, hsp.evalue)
        self.assertEqual(139.0, hsp.bitscore)
        self.assertEqual(350.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(81, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(102, hsp.query_end)
        self.assertEqual(100, hsp.hit_end)
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEG', str(hsp.query.seq)[:40])
        self.assertEqual('MKK    LFFILLL+GCGV ++KSQGED      + TKEG', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEG', str(hsp.hit.seq)[:40])
        self.assertEqual('DITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', str(hsp.query.seq)[-40:])
        self.assertEqual('DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('DITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN', str(hsp.hit.seq)[-40:])

    def test_text_2226_blastp_003(self):
        """Test parsing blastp output (text_2226_blastp_003.txt)"""

        blast_file = get_file('text_2226_blastp_003.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual('NCBI Protein Reference Sequences', qresult.target)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|NP_062422.1|', hit.id)
        self.assertEqual('pleckstrin [Mus musculus]', hit.description)
        self.assertEqual(350, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(1e-63, hsp.evalue)
        self.assertEqual(205.0, hsp.bitscore)
        self.assertEqual(522.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(101, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # first qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(100, hsp.aln_span)
        self.assertEqual(0.002, hsp.evalue)
        self.assertEqual(43.5, hsp.bitscore)
        self.assertEqual(101.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(245, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(345, hsp.hit_end)
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGED', str(hsp.hit.seq)[:40])
        self.assertEqual('FGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('VKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|XP_003502426.1|', hit.id)
        self.assertEqual('PREDICTED: pleckstrin-like [Cricetulus griseus]', hit.description)
        self.assertEqual(350, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(2e-63, hsp.evalue)
        self.assertEqual(205.0, hsp.bitscore)
        self.assertEqual(521.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(101, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # first qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(100, hsp.aln_span)
        self.assertEqual(0.001, hsp.evalue)
        self.assertEqual(43.9, hsp.bitscore)
        self.assertEqual(102.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(30, hsp.ident_num)
        self.assertEqual(50, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(245, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(345, hsp.hit_end)
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGED', str(hsp.hit.seq)[:40])
        self.assertEqual('GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('GK+     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('GKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])

    def test_text_2226_blastp_004(self):
        """Test parsing blastp output (text_2226_blastp_004.txt)"""

        blast_file = get_file('text_2226_blastp_004.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(3, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual('minirefseq_prot', qresult.target)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillussubtilis subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual('minirefseq_prot', qresult.target)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(5, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|308175296|ref|YP_003922001.1|', hit.id)
        self.assertEqual('membrane bound lipoprotein [Bacillus amyloliquefaciens DSM 7]', hit.description)
        self.assertEqual(100, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(102, hsp.aln_span)
        self.assertEqual(2e-46, hsp.evalue)
        self.assertEqual(139.0, hsp.bitscore)
        self.assertEqual(350.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(2, hsp.gap_num)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(81, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(102, hsp.query_end)
        self.assertEqual(100, hsp.hit_end)
        self.assertEqual('MKKFIALLFFILLLSGCGVNSQKSQGEDVSPDSNIETKEG', str(hsp.query.seq)[:40])
        self.assertEqual('MKK    LFFILLL+GCGV ++KSQGED      + TKEG', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('MKKIFGCLFFILLLAGCGVTNEKSQGEDAG--EKLVTKEG', str(hsp.hit.seq)[:40])
        self.assertEqual('DITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERAN', str(hsp.query.seq)[-40:])
        self.assertEqual('DITEES  D+   N+G+KVT+ Y+KN +GQL+LKDIE AN', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('DITEESADDVKNLNNGEKVTVKYQKNSKGQLVLKDIEPAN', str(hsp.hit.seq)[-40:])
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|375363999|ref|YP_005132038.1|', hit.id)
        self.assertEqual('lytA gene product [Bacillus amyloliquefaciens subsp. plantarum CAU B946]', hit.description)
        self.assertEqual(105, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(105, hsp.aln_span)
        self.assertEqual(7e-27, hsp.evalue)
        self.assertEqual(89.0, hsp.bitscore)
        self.assertEqual(219.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(48, hsp.ident_num)
        self.assertEqual(69, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(101, hsp.query_end)
        self.assertEqual(104, hsp.hit_end)
        self.assertEqual('MKKFIALLFFILL----LSGCGVNSQKSQGEDVSPDSNIE', str(hsp.query.seq)[:40])
        self.assertEqual('MKK IA  F ILL    L+ CG   Q  +G   S ++  +', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('MKKTIAASFLILLFSVVLAACGTAEQSKKGSG-SSENQAQ', str(hsp.hit.seq)[:40])
        self.assertEqual('LDITEESTSDLDKFNSGDKVTITYEKNDEGQLLLKDIERA', str(hsp.query.seq)[-40:])
        self.assertEqual(' + +++ +  L+KF+  DKV+ITY  ND+GQ  +K+IE+A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FEFSDDFSDVLNKFSENDKVSITYFTNDKGQKEIKEIEKA', str(hsp.hit.seq)[-40:])

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual('minirefseq_prot', qresult.target)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(5, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|11464971|ref|NP_062422.1|', hit.id)
        self.assertEqual('pleckstrin [Mus musculus]', hit.description)
        self.assertEqual(350, hit.seq_len)
        self.assertEqual(2, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(2e-69, hsp.evalue)
        self.assertEqual(205.0, hsp.bitscore)
        self.assertEqual(522.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(101, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # third qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(100, hsp.aln_span)
        self.assertEqual(3e-09, hsp.evalue)
        self.assertEqual(43.5, hsp.bitscore)
        self.assertEqual(101.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(29, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(245, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(345, hsp.hit_end)
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGED', str(hsp.hit.seq)[:40])
        self.assertEqual('FGK--RMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('  K     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('VKKSDEENLFEIITADEVHYYLQAATSKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])
        # third qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|354480464|ref|XP_003502426.1|', hit.id)
        self.assertEqual('PREDICTED: pleckstrin-like [Cricetulus griseus]', hit.description)
        self.assertEqual(350, hit.seq_len)
        self.assertEqual(2, len(hit))
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(3e-69, hsp.evalue)
        self.assertEqual(205.0, hsp.bitscore)
        self.assertEqual(521.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(3, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(101, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # third qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(100, hsp.aln_span)
        self.assertEqual(2e-09, hsp.evalue)
        self.assertEqual(43.9, hsp.bitscore)
        self.assertEqual(102.0, hsp.bitscore_raw)
        self.assertEqual(0, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(30, hsp.ident_num)
        self.assertEqual(50, hsp.pos_num)
        self.assertEqual(0, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(245, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(345, hsp.hit_end)
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGED', str(hsp.hit.seq)[:40])
        self.assertEqual('GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('GK+     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('GKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])


class BlastxCases(BaseBlastCases):

    def test_text_2226_blastx_001(self):
        """Test parsing blastx output (text_2226_blastx_001.txt)"""

        blast_file = get_file('text_2226_blastx_001.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual('NCBI Protein Reference Sequences', qresult.target)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

    def test_text_2226_blastx_002(self):
        """Test parsing blastx output (text_2226_blastx_002.txt)"""

        blast_file = get_file('text_2226_blastx_002.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|356995852:1-490', qresult.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcriptionfactor 1 (Pou5f1), transcript variant 1, mRNA', qresult.description)
        self.assertEqual(490, qresult.seq_len)
        self.assertEqual('NCBI Protein Reference Sequences', qresult.target)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|NP_038661.2|', hit.id)
        self.assertEqual('POU domain, class 5, transcription factor 1 isoform 1 [Mus musculus]', hit.description)
        self.assertEqual(352, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(140, hsp.aln_span)
        self.assertEqual(4e-57, hsp.evalue)
        self.assertEqual(192.0, hsp.bitscore)
        self.assertEqual(487.0, hsp.bitscore_raw)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(140, hsp.ident_num)
        self.assertEqual(140, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(68, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(488, hsp.query_end)
        self.assertEqual(140, hsp.hit_end)
        self.assertEqual('MAGHLAsdfafspppgggdgsagLEPGWVDPRTWLSFQgp', str(hsp.query.seq)[:40])
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGP', str(hsp.hit.seq)[:40])
        self.assertEqual('NSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', str(hsp.query.seq)[-40:])
        self.assertEqual('NSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('NSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|NP_001009178.1|', hit.id)
        self.assertEqual('POU class 5 homeobox 1 [Rattus norvegicus]', hit.description)
        self.assertEqual(352, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(140, hsp.aln_span)
        self.assertEqual(3e-52, hsp.evalue)
        self.assertEqual(179.0, hsp.bitscore)
        self.assertEqual(454.0, hsp.bitscore_raw)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(133, hsp.ident_num)
        self.assertEqual(135, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(68, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(488, hsp.query_end)
        self.assertEqual(140, hsp.hit_end)
        self.assertEqual('MAGHLAsdfafspppgggdgsagLEPGWVDPRTWLSFQgp', str(hsp.query.seq)[:40])
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('MAGHLASDFAFSPPPGGGDGSAGLEPGWVDPRTWLSFQGP', str(hsp.hit.seq)[:40])
        self.assertEqual('NSEGTSSEPCADRPNAVKLEKVEPTPEESQDMKALQKELE', str(hsp.query.seq)[-40:])
        self.assertEqual('NSEG SS PC  RP+AVKLEKVEP+PEESQDMKALQKELE', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('NSEGASSGPCTARPSAVKLEKVEPSPEESQDMKALQKELE', str(hsp.hit.seq)[-40:])

    def test_text_2226_blastx_003(self):
        """Test parsing blastx output (text_2226_blastx_003.txt)"""

        blast_file = get_file('text_2226_blastx_003.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 strand=+repeatMasking=none', qresult.description)
        self.assertEqual(485, qresult.seq_len)
        self.assertEqual('NCBI Protein Reference Sequences', qresult.target)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|XP_003278367.1|', hit.id)
        self.assertEqual('PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]', hit.description)
        self.assertEqual(132, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(95, hsp.aln_span)
        self.assertEqual(2e-32, hsp.evalue)
        self.assertEqual(121.0, hsp.bitscore)
        self.assertEqual(304.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(15, hsp.query_start)
        self.assertEqual(24, hsp.hit_start)
        self.assertEqual(300, hsp.query_end)
        self.assertEqual(119, hsp.hit_end)
        self.assertEqual('LRRSFALVAQAGVQWLDLGppqpppPGFK*FSCLSHPSSW', str(hsp.query.seq)[:40])
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSW', str(hsp.hit.seq)[:40])
        self.assertEqual('VETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', str(hsp.query.seq)[-40:])
        self.assertEqual('VE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('VEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', str(hsp.hit.seq)[-40:])
        # first qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(72, hsp.aln_span)
        self.assertEqual(2e-06, hsp.evalue)
        self.assertEqual(51.6, hsp.bitscore)
        self.assertEqual(122.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(243, hsp.query_start)
        self.assertEqual(31, hsp.hit_start)
        self.assertEqual(459, hsp.query_end)
        self.assertEqual(98, hsp.hit_end)
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPP', str(hsp.query.seq)[:40])
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPP', str(hsp.hit.seq)[:40])
        self.assertEqual('WDCRCPPPHPANffffffffFLRRSFALVAQAGVQWLDLG', str(hsp.query.seq)[-40:])
        self.assertEqual('W+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('WEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|NP_001243358.1|', hit.id)
        self.assertEqual('PDZ and LIM domain protein 5 isoform i [Homo sapiens]', hit.description)
        self.assertEqual(136, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(88, hsp.aln_span)
        self.assertEqual(1e-29, hsp.evalue)
        self.assertEqual(114.0, hsp.bitscore)
        self.assertEqual(286.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(63, hsp.ident_num)
        self.assertEqual(69, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(15, hsp.query_start)
        self.assertEqual(29, hsp.hit_start)
        self.assertEqual(279, hsp.query_end)
        self.assertEqual(117, hsp.hit_end)
        self.assertEqual('VAQAGVQWLDLGppqpppPGFK*FSCLSHPSSWDYRHMPP', str(hsp.query.seq)[:40])
        self.assertEqual('++ AGVQW +LG PQPP P FK FSCLS PSSWDYRH+PP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('ISSAGVQWRNLGSPQPPSPEFKRFSCLSLPSSWDYRHVPP', str(hsp.hit.seq)[:40])
        self.assertEqual('VETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', str(hsp.query.seq)[-40:])
        self.assertEqual('VET F +VGQAGLE P SG+LP  ASQS  ITGVSH A P', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('VETKFPYVGQAGLELPTSGDLPTSASQSAKITGVSHRAWP', str(hsp.hit.seq)[-40:])
        # first qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(69, hsp.aln_span)
        self.assertEqual(1e-06, hsp.evalue)
        self.assertEqual(52.4, hsp.bitscore)
        self.assertEqual(124.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(33, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(258, hsp.query_start)
        self.assertEqual(27, hsp.hit_start)
        self.assertEqual(465, hsp.query_end)
        self.assertEqual(91, hsp.hit_end)
        self.assertEqual('VSVGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCP', str(hsp.query.seq)[:40])
        self.assertEqual('+++  A VQ  +L S QPP+PEFK FS LSL SSWD R  ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('LTISSAGVQWRNLGSPQPPSPEFKRFSCLSLPSSWDYRHV', str(hsp.hit.seq)[:40])
        self.assertEqual('SLQSSWDCRCPPPHPANffffffffFLRRSFALVAQAGVQ', str(hsp.query.seq)[-40:])
        self.assertEqual('SL SSWD R  PP  AN     F F +   F  V QAG++', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('SLPSSWDYRHVPPRLAN-----FVFLVETKFPYVGQAGLE', str(hsp.hit.seq)[-40:])

    def test_text_2226_blastx_004(self):
        """Test parsing blastx output (text_2226_blastx_004.txt)"""

        blast_file = get_file('text_2226_blastx_004.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(2, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual('minirefseq_prot', qresult.target)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 strand=+repeatMasking=none', qresult.description)
        self.assertEqual(485, qresult.seq_len)
        self.assertEqual('minirefseq_prot', qresult.target)
        self.assertEqual('blastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(5, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|332258565|ref|XP_003278367.1|', hit.id)
        self.assertEqual('PREDICTED: UPF0764 protein C16orf89-like [Nomascus leucogenys]', hit.description)
        self.assertEqual(132, hit.seq_len)
        self.assertEqual(2, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(95, hsp.aln_span)
        self.assertEqual(3e-38, hsp.evalue)
        self.assertEqual(121.0, hsp.bitscore)
        self.assertEqual(304.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(69, hsp.ident_num)
        self.assertEqual(74, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(15, hsp.query_start)
        self.assertEqual(24, hsp.hit_start)
        self.assertEqual(300, hsp.query_end)
        self.assertEqual(119, hsp.hit_end)
        self.assertEqual('LRRSFALVAQAGVQWLDLGppqpppPGFK*FSCLSHPSSW', str(hsp.query.seq)[:40])
        self.assertEqual('LRRSFALVAQ  VQW +LG PQPPPPGFK FSCLS  SSW', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('LRRSFALVAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSW', str(hsp.hit.seq)[:40])
        self.assertEqual('VETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', str(hsp.query.seq)[-40:])
        self.assertEqual('VE GF HVGQAGLE   SG+ P   SQS GI GVSH AQP', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('VEMGFLHVGQAGLELVTSGDPPTLTSQSAGIIGVSHCAQP', str(hsp.hit.seq)[-40:])
        # second qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(72, hsp.aln_span)
        self.assertEqual(3e-12, hsp.evalue)
        self.assertEqual(51.6, hsp.bitscore)
        self.assertEqual(122.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(34, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(243, hsp.query_start)
        self.assertEqual(31, hsp.hit_start)
        self.assertEqual(459, hsp.query_end)
        self.assertEqual(98, hsp.hit_end)
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPP', str(hsp.query.seq)[:40])
        self.assertEqual('V   RVQ ++L S QPP P FK FS LSL SSW+ R  PP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('VAQTRVQWYNLGSPQPPPPGFKRFSCLSLLSSWEYRHVPP', str(hsp.hit.seq)[:40])
        self.assertEqual('WDCRCPPPHPANffffffffFLRRSFALVAQAGVQWLDLG', str(hsp.query.seq)[-40:])
        self.assertEqual('W+ R  PPH AN     F F +   F  V QAG++ +  G', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('WEYRHVPPHLAN-----FLFLVEMGFLHVGQAGLELVTSG', str(hsp.hit.seq)[-40:])
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|374093214|ref|NP_001243358.1|', hit.id)
        self.assertEqual('PDZ and LIM domain protein 5 isoform i [Homo sapiens]', hit.description)
        self.assertEqual(136, hit.seq_len)
        self.assertEqual(2, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(88, hsp.aln_span)
        self.assertEqual(2e-35, hsp.evalue)
        self.assertEqual(114.0, hsp.bitscore)
        self.assertEqual(286.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(63, hsp.ident_num)
        self.assertEqual(69, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(15, hsp.query_start)
        self.assertEqual(29, hsp.hit_start)
        self.assertEqual(279, hsp.query_end)
        self.assertEqual(117, hsp.hit_end)
        self.assertEqual('VAQAGVQWLDLGppqpppPGFK*FSCLSHPSSWDYRHMPP', str(hsp.query.seq)[:40])
        self.assertEqual('++ AGVQW +LG PQPP P FK FSCLS PSSWDYRH+PP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('ISSAGVQWRNLGSPQPPSPEFKRFSCLSLPSSWDYRHVPP', str(hsp.hit.seq)[:40])
        self.assertEqual('VETGFYHVGQAGLEPPISGNLPAWASQSVGITGVSHHAQP', str(hsp.query.seq)[-40:])
        self.assertEqual('VET F +VGQAGLE P SG+LP  ASQS  ITGVSH A P', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('VETKFPYVGQAGLELPTSGDLPTSASQSAKITGVSHRAWP', str(hsp.hit.seq)[-40:])
        # second qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(69, hsp.aln_span)
        self.assertEqual(2e-12, hsp.evalue)
        self.assertEqual(52.4, hsp.bitscore)
        self.assertEqual(124.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(5, hsp.gap_num)
        self.assertEqual(33, hsp.ident_num)
        self.assertEqual(41, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(258, hsp.query_start)
        self.assertEqual(27, hsp.hit_start)
        self.assertEqual(465, hsp.query_end)
        self.assertEqual(91, hsp.hit_end)
        self.assertEqual('VSVGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCP', str(hsp.query.seq)[:40])
        self.assertEqual('+++  A VQ  +L S QPP+PEFK FS LSL SSWD R  ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('LTISSAGVQWRNLGSPQPPSPEFKRFSCLSLPSSWDYRHV', str(hsp.hit.seq)[:40])
        self.assertEqual('SLQSSWDCRCPPPHPANffffffffFLRRSFALVAQAGVQ', str(hsp.query.seq)[-40:])
        self.assertEqual('SL SSWD R  PP  AN     F F +   F  V QAG++', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('SLPSSWDYRHVPPRLAN-----FVFLVETKFPYVGQAGLE', str(hsp.hit.seq)[-40:])


class TblastnCases(BaseBlastCases):

    def test_text_2226_tblastn_001(self):
        """Test parsing tblastn output (text_2226_tblastn_001.txt)"""

        blast_file = get_file('text_2226_tblastn_001.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

    def test_text_2226_tblastn_002(self):
        """Test parsing tblastn output (text_2226_tblastn_002.txt)"""

        blast_file = get_file('text_2226_tblastn_002.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillussubtilis subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(4, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|XM_001425911.1|', hit.id)
        self.assertEqual('Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA', hit.description)
        self.assertEqual(4632, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(43, hsp.aln_span)
        self.assertEqual(0.74, hsp.evalue)
        self.assertEqual(34.7, hsp.bitscore)
        self.assertEqual(78.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(30, hsp.query_start)
        self.assertEqual(1743, hsp.hit_start)
        self.assertEqual(73, hsp.query_end)
        self.assertEqual(1872, hsp.hit_end)
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTS', str(hsp.query.seq)[:40])
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*K', str(hsp.hit.seq)[:40])
        self.assertEqual('NIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', str(hsp.query.seq)[-40:])
        self.assertEqual('   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('ATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|XM_003382561.1|', hit.id)
        self.assertEqual('PREDICTED: Amphimedon queenslandica CWF19-like protein 1-like (LOC100635130), mRNA', hit.description)
        self.assertEqual(1811, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(74, hsp.aln_span)
        self.assertEqual(6.4, hsp.evalue)
        self.assertEqual(32.0, hsp.bitscore)
        self.assertEqual(71.0, hsp.bitscore_raw)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(11, hsp.gap_num)
        self.assertEqual(19, hsp.ident_num)
        self.assertEqual(36, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(28, hsp.query_start)
        self.assertEqual(1105, hsp.hit_start)
        self.assertEqual(94, hsp.query_end)
        self.assertEqual(1318, hsp.hit_end)
        self.assertEqual('VSPDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEES', str(hsp.query.seq)[:40])
        self.assertEqual('+  DS +   +G   GL D H + + + + P S+D  +E ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('IGNDSYLALSKG---GLVDEHVLILPIGHYPSSIDAPQEV', str(hsp.hit.seq)[:40])
        self.assertEqual('DITEESTSDLDK--------FNSGDKVTITYEKNDEGQLL', str(hsp.query.seq)[-40:])
        self.assertEqual('D  +E   ++DK        F+S ++  + +E+N   Q L', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('DAPQEVIEEIDKFKVALRKYFSSKNQTCVMFERNFRSQHL', str(hsp.hit.seq)[-40:])

    def test_text_2226_tblastn_003(self):
        """Test parsing tblastn output (text_2226_tblastn_003.txt)"""

        blast_file = get_file('text_2226_tblastn_003.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|XM_003502378.1|', hit.id)
        self.assertEqual('PREDICTED: Cricetulus griseus pleckstrin-like (LOC100773128), mRNA', hit.description)
        self.assertEqual(1119, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(1e-63, hsp.evalue)
        self.assertEqual(205.0, hsp.bitscore)
        self.assertEqual(521.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(98, hsp.ident_num)
        self.assertEqual(98, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(75, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(369, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # first qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(100, hsp.aln_span)
        self.assertEqual(0.0005, hsp.evalue)
        self.assertEqual(43.9, hsp.bitscore)
        self.assertEqual(102.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(30, hsp.ident_num)
        self.assertEqual(50, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(801, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(1101, hsp.hit_end)
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGED', str(hsp.hit.seq)[:40])
        self.assertEqual('GKRM---FVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('GK+     + +I T  +  ++ QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('GKKSDDENLFEIITADEVHYYLQAAAPKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|XM_003360601.2|', hit.id)
        self.assertEqual('PREDICTED: Sus scrofa pleckstrin-like (LOC100626968), mRNA', hit.description)
        self.assertEqual(772, hit.seq_len)
        self.assertEqual(2, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(1e-62, hsp.evalue)
        self.assertEqual(199.0, hsp.bitscore)
        self.assertEqual(506.0, hsp.bitscore_raw)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(94, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(388, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # first qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(71, hsp.aln_span)
        self.assertEqual(2.8, hsp.evalue)
        self.assertEqual(32.7, hsp.bitscore)
        self.assertEqual(73.0, hsp.bitscore_raw)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(29, hsp.query_start)
        self.assertEqual(541, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(754, hsp.hit_end)
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK-', str(hsp.query.seq)[:40])
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ + ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWER', str(hsp.hit.seq)[:40])
        self.assertEqual('GKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('GKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])

    def test_text_2226_tblastn_004(self):
        """Test parsing tblastn output (text_2226_tblastn_004.txt)"""

        blast_file = get_file('text_2226_tblastn_004.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(3, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(32, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|16080617|ref|NP_391444.1|', qresult.id)
        self.assertEqual('membrane bound lipoprotein [Bacillussubtilis subsp. subtilis str. 168]', qresult.description)
        self.assertEqual(102, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(3, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|145479850|ref|XM_001425911.1|', hit.id)
        self.assertEqual('Paramecium tetraurelia hypothetical protein (GSPATT00004923001) partial mRNA', hit.description)
        self.assertEqual(4632, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(43, hsp.aln_span)
        self.assertEqual(1e-05, hsp.evalue)
        self.assertEqual(34.7, hsp.bitscore)
        self.assertEqual(78.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(15, hsp.ident_num)
        self.assertEqual(26, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(30, hsp.query_start)
        self.assertEqual(1743, hsp.hit_start)
        self.assertEqual(73, hsp.query_end)
        self.assertEqual(1872, hsp.hit_end)
        self.assertEqual('PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTS', str(hsp.query.seq)[:40])
        self.assertEqual('P +   TK+GT +GL   HTI   + +  +SL++ E++  ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*K', str(hsp.hit.seq)[:40])
        self.assertEqual('NIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD', str(hsp.query.seq)[-40:])
        self.assertEqual('   TK+GT +GL   HTI   + +  +SL++ E++  D+D', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('ATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID', str(hsp.hit.seq)[-40:])
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|72012412|ref|XM_777959.1|', hit.id)
        self.assertEqual('PREDICTED: Strongylocentrotus purpuratus hypothetical LOC577746 (LOC577746), mRNA', hit.description)
        self.assertEqual(1593, hit.seq_len)
        self.assertEqual(1, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(59, hsp.aln_span)
        self.assertEqual(0.0001, hsp.evalue)
        self.assertEqual(31.6, hsp.bitscore)
        self.assertEqual(70.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(8, hsp.gap_num)
        self.assertEqual(20, hsp.ident_num)
        self.assertEqual(29, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(43, hsp.query_start)
        self.assertEqual(1056, hsp.hit_start)
        self.assertEqual(94, hsp.query_end)
        self.assertEqual(1233, hsp.hit_end)
        self.assertEqual('GLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSG-----', str(hsp.query.seq)[:40])
        self.assertEqual('GL   HT+ + V +    LD+TEE  ++LD+F S      ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('GLVPDHTLILPVGHYQSMLDLTEEVQTELDQFKSALRKYY', str(hsp.hit.seq)[:40])
        self.assertEqual('DITEESTSDLDKFNSG--------DKVTITYEKNDEGQLL', str(hsp.query.seq)[-40:])
        self.assertEqual('D+TEE  ++LD+F S          K  + YE+N   Q L', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('DLTEEVQTELDQFKSALRKYYLSKGKTCVIYERNFRTQHL', str(hsp.hit.seq)[-40:])

        # test third qresult
        qresult = qresults[2]
        self.assertEqual('gi|11464971:4-101', qresult.id)
        self.assertEqual('pleckstrin [Mus musculus]', qresult.description)
        self.assertEqual(98, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('tblastn', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(5, len(qresult))
        # third qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|350596019|ref|XM_003360601.2|', hit.id)
        self.assertEqual('PREDICTED: Sus scrofa pleckstrin-like (LOC100626968), mRNA', hit.description)
        self.assertEqual(772, hit.seq_len)
        self.assertEqual(2, len(hit))
        # third qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(2e-67, hsp.evalue)
        self.assertEqual(199.0, hsp.bitscore)
        self.assertEqual(506.0, hsp.bitscore_raw)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(94, hsp.ident_num)
        self.assertEqual(96, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(94, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(388, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVKKGS+FNTWKPMWV+LLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFV KITTTKQQDHFFQAAFLEERD WVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # third qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(71, hsp.aln_span)
        self.assertEqual(4e-05, hsp.evalue)
        self.assertEqual(32.7, hsp.bitscore)
        self.assertEqual(73.0, hsp.bitscore_raw)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(4, hsp.gap_num)
        self.assertEqual(21, hsp.ident_num)
        self.assertEqual(33, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(29, hsp.query_start)
        self.assertEqual(541, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(754, hsp.hit_end)
        self.assertEqual('IEFYKKKSDNSPKGMIPLKGSTLTS-PCQDFGKRMFVLK-', str(hsp.query.seq)[:40])
        self.assertEqual('+ +Y       P G I L+G  +TS      GK  F+ + ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWER', str(hsp.hit.seq)[:40])
        self.assertEqual('GKRMFVLK---ITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('GK  F+ +     T  +  +F QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('GKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])
        # third qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|301779869|ref|XM_002925302.1|', hit.id)
        self.assertEqual('PREDICTED: Ailuropoda melanoleuca pleckstrin-like (LOC100466932), mRNA', hit.description)
        self.assertEqual(1144, hit.seq_len)
        self.assertEqual(2, len(hit))
        # third qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(98, hsp.aln_span)
        self.assertEqual(2e-67, hsp.evalue)
        self.assertEqual(202.0, hsp.bitscore)
        self.assertEqual(515.0, hsp.bitscore_raw)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(96, hsp.ident_num)
        self.assertEqual(97, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(77, hsp.hit_start)
        self.assertEqual(98, hsp.query_end)
        self.assertEqual(371, hsp.hit_end)
        self.assertEqual('KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('KRIREGYLVK+GSVFNTWKPMWVVLLEDGIEFYKKKSDNS', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNS', str(hsp.hit.seq)[:40])
        self.assertEqual('FGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.query.seq)[-40:])
        self.assertEqual('FGKRMFV KITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('FGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK', str(hsp.hit.seq)[-40:])
        # third qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(100, hsp.aln_span)
        self.assertEqual(3e-09, hsp.evalue)
        self.assertEqual(45.1, hsp.bitscore)
        self.assertEqual(105.0, hsp.bitscore_raw)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(0, hsp.hit_frame)
        self.assertEqual(6, hsp.gap_num)
        self.assertEqual(30, hsp.ident_num)
        self.assertEqual(48, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(0, hsp.hit_strand)
        self.assertEqual(2, hsp.query_start)
        self.assertEqual(803, hsp.hit_start)
        self.assertEqual(96, hsp.query_end)
        self.assertEqual(1103, hsp.hit_end)
        self.assertEqual('IREGYLVKKGSVFNTWKPMWVVLLEDG--IEFYKKKSDNS', str(hsp.query.seq)[:40])
        self.assertEqual('I++G L+K+G     WK    +L ED   + +Y       ', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGED', str(hsp.hit.seq)[:40])
        self.assertEqual('QDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA', str(hsp.query.seq)[-40:])
        self.assertEqual('    +   + +I T  +  +F QAA  +ER  W++ I+ A', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('VRKSEEENLFEIITADEVHYFLQAATPKERTEWIKAIQVA', str(hsp.hit.seq)[-40:])


class TblastxCases(BaseBlastCases):

    def test_text_2226_tblastx_001(self):
        """Test parsing tblastx output (text_2226_tblastx_001.txt)"""

        blast_file = get_file('text_2226_tblastx_001.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

    def test_text_2226_tblastx_002(self):
        """Test parsing tblastx output (text_2226_tblastx_002.txt)"""

        blast_file = get_file('text_2226_tblastx_002.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('gi|356995852:1-490', qresult.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcriptionfactor 1 (Pou5f1), transcript variant 1, mRNA', qresult.description)
        self.assertEqual(490, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|NM_013633.3|', hit.id)
        self.assertEqual('Mus musculus POU domain, class 5, transcription factor 1 (Pou5f1), transcript variant 1, mRNA', hit.description)
        self.assertEqual(1353, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(163, hsp.aln_span)
        self.assertEqual(2e-115, hsp.evalue)
        self.assertEqual(418.0, hsp.bitscore)
        self.assertEqual(908.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(163, hsp.ident_num)
        self.assertEqual(163, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(0, hsp.hit_start)
        self.assertEqual(489, hsp.query_end)
        self.assertEqual(489, hsp.hit_end)
        self.assertEqual('EVKPSLGEPSFHQAPGSGCPPSPWLDTWLQTSPSHPHQVG', str(hsp.query.seq)[:40])
        self.assertEqual('EVKPSLGEPSFHQAPGSGCPPSPWLDTWLQTSPSHPHQVG', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('EVKPSLGEPSFHQAPGSGCPPSPWLDTWLQTSPSHPHQVG', str(hsp.hit.seq)[:40])
        self.assertEqual('TQREPPLSPVPTAPMP*SWRRWNQLPRSPRT*KPCRRS*N', str(hsp.query.seq)[-40:])
        self.assertEqual('TQREPPLSPVPTAPMP*SWRRWNQLPRSPRT*KPCRRS*N', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('TQREPPLSPVPTAPMP*SWRRWNQLPRSPRT*KPCRRS*N', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|XR_141831.1|', hit.id)
        self.assertEqual('PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA ref|XR_105837.2| PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA ref|XR_141464.1| PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA ref|XR_141446.1| PREDICTED: Mus musculus predicted gene, 19553 (Gm19553), miscRNA', hit.description)
        self.assertEqual(570, hit.seq_len)
        self.assertEqual(1, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(163, hsp.aln_span)
        self.assertEqual(3e-114, hsp.evalue)
        self.assertEqual(415.0, hsp.bitscore)
        self.assertEqual(900.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(162, hsp.ident_num)
        self.assertEqual(162, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(-1, hsp.hit_strand)
        self.assertEqual(0, hsp.query_start)
        self.assertEqual(81, hsp.hit_start)
        self.assertEqual(489, hsp.query_end)
        self.assertEqual(570, hsp.hit_end)
        self.assertEqual('EVKPSLGEPSFHQAPGSGCPPSPWLDTWLQTSPSHPHQVG', str(hsp.query.seq)[:40])
        self.assertEqual('EVKPSLGEPSFHQAPGSGCPPSPWLDTWLQTSPSHPHQVG', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('EVKPSLGEPSFHQAPGSGCPPSPWLDTWLQTSPSHPHQVG', str(hsp.hit.seq)[:40])
        self.assertEqual('TQREPPLSPVPTAPMP*SWRRWNQLPRSPRT*KPCRRS*N', str(hsp.query.seq)[-40:])
        self.assertEqual('TQREPPLSPVPTAPMP*SWRRWNQL RSPRT*KPCRRS*N', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('TQREPPLSPVPTAPMP*SWRRWNQLQRSPRT*KPCRRS*N', str(hsp.hit.seq)[-40:])

    def test_text_2226_tblastx_003(self):
        """Test parsing tblastx output (text_2226_tblastx_003.txt)"""

        blast_file = get_file('text_2226_tblastx_003.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('hg19_dna', qresult.id)
        self.assertEqual('range=chr1:1207057-1207541 5\'pad=0 3\'pad=0 strand=+repeatMasking=none', qresult.description)
        self.assertEqual(485, qresult.seq_len)
        self.assertEqual('NCBI Transcript Reference Sequences', qresult.target)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(10, len(qresult))
        # first qresult, first hit
        hit = qresult[0]
        self.assertEqual('ref|NM_002985.2|', hit.id)
        self.assertEqual('Homo sapiens chemokine (C-C motif) ligand 5 (CCL5), mRNA', hit.description)
        self.assertEqual(1237, hit.seq_len)
        self.assertEqual(3, len(hit))
        # first qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(107, hsp.aln_span)
        self.assertEqual(4e-49, hsp.evalue)
        self.assertEqual(118.0, hsp.bitscore)
        self.assertEqual(252.0, hsp.bitscore_raw)
        self.assertEqual(-3, hsp.query_frame)
        self.assertEqual(-1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(68, hsp.ident_num)
        self.assertEqual(72, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(-1, hsp.hit_strand)
        self.assertEqual(138, hsp.query_start)
        self.assertEqual(622, hsp.hit_start)
        self.assertEqual(459, hsp.query_end)
        self.assertEqual(943, hsp.hit_end)
        self.assertEqual('VGPARVQ*HDLSSLQPPAPEFK*FSHLSLQSSWDCRCPPP', str(hsp.query.seq)[:40])
        self.assertEqual('V  A V+ H+LSSLQPP P FK FS LSL SSWD R  PP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('VTQAGVKWHNLSSLQPPPPGFKQFSCLSLPSSWDYRRGPP', str(hsp.hit.seq)[:40])
        self.assertEqual('WLDLGppqpppPGFK*FSCLSHPSSWDYRHMPPCLINFVF', str(hsp.query.seq)[-40:])
        self.assertEqual('W DLG  Q PPPGF  FSCLS PSSWDYR   P   NF++', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('WRDLGSLQAPPPGFTPFSCLSLPSSWDYRRPLPRPANFLY', str(hsp.hit.seq)[-40:])
        # first qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(44, hsp.aln_span)
        self.assertEqual(4e-49, hsp.evalue)
        self.assertEqual(100.0, hsp.bitscore)
        self.assertEqual(214.0, hsp.bitscore_raw)
        self.assertEqual(-2, hsp.query_frame)
        self.assertEqual(-2, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(37, hsp.ident_num)
        self.assertEqual(38, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(-1, hsp.hit_strand)
        self.assertEqual(16, hsp.query_start)
        self.assertEqual(498, hsp.hit_start)
        self.assertEqual(148, hsp.query_end)
        self.assertEqual(630, hsp.hit_end)
        self.assertEqual('FCIFSRDGVLPCWSGWSRTPDLR*SACLGLPKCWDYRCEP', str(hsp.query.seq)[:40])
        self.assertEqual('FCIFSRDGV  CW GWSRTPDL+*S  LGLPKCWDYR EP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('FCIFSRDGVSSCWPGWSRTPDLK*STHLGLPKCWDYRREP', str(hsp.hit.seq)[:40])
        self.assertEqual('SRDGVLPCWSGWSRTPDLR*SACLGLPKCWDYRCEPPRPA', str(hsp.query.seq)[-40:])
        self.assertEqual('SRDGV  CW GWSRTPDL+*S  LGLPKCWDYR EPPRPA', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('SRDGVSSCWPGWSRTPDLK*STHLGLPKCWDYRREPPRPA', str(hsp.hit.seq)[-40:])
        # first qresult, second hit
        hit = qresult[1]
        self.assertEqual('ref|XM_003255417.1|', hit.id)
        self.assertEqual('PREDICTED: Nomascus leucogenys 5\'-nucleotidase, cytosolic II, transcript variant 2 (NT5C2), mRNA', hit.description)
        self.assertEqual(3285, hit.seq_len)
        self.assertEqual(3, len(hit))
        # first qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(94, hsp.aln_span)
        self.assertEqual(9e-49, hsp.evalue)
        self.assertEqual(197.0, hsp.bitscore)
        self.assertEqual(425.0, hsp.bitscore_raw)
        self.assertEqual(-2, hsp.query_frame)
        self.assertEqual(-2, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(78, hsp.ident_num)
        self.assertEqual(79, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(-1, hsp.hit_strand)
        self.assertEqual(16, hsp.query_start)
        self.assertEqual(2744, hsp.hit_start)
        self.assertEqual(298, hsp.query_end)
        self.assertEqual(3026, hsp.hit_end)
        self.assertEqual('ETEFRSCCPGWSAMA*SWPTTASTSWIQVILLPQSPE*LG', str(hsp.query.seq)[:40])
        self.assertEqual('E EFRSCCPGWSAMA SW    S SW+QVIL PQ PE*LG', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('EMEFRSCCPGWSAMAQSWLIATSVSWVQVILWPQPPE*LG', str(hsp.hit.seq)[:40])
        self.assertEqual('SRDGVLPCWSGWSRTPDLR*SACLGLPKCWDYRCEPPRPA', str(hsp.query.seq)[-40:])
        self.assertEqual('SRDGV PCWSGWSRTPDLR*SACLGLPKCWDYR EPP PA', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('SRDGVSPCWSGWSRTPDLR*SACLGLPKCWDYRREPPCPA', str(hsp.hit.seq)[-40:])
        # first qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(94, hsp.aln_span)
        self.assertEqual(4e-43, hsp.evalue)
        self.assertEqual(178.0, hsp.bitscore)
        self.assertEqual(384.0, hsp.bitscore_raw)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(77, hsp.ident_num)
        self.assertEqual(83, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(17, hsp.query_start)
        self.assertEqual(2745, hsp.hit_start)
        self.assertEqual(299, hsp.query_end)
        self.assertEqual(3027, hsp.hit_end)
        self.assertEqual('AGRGGSHL*SQHFGRPRQADYLRSGVRDQPDQHGKTPSLL', str(hsp.query.seq)[:40])
        self.assertEqual('AG GGS L*SQHFGRPRQAD+LRSGVRDQPDQHG+TPSLL', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('AGHGGSRL*SQHFGRPRQADHLRSGVRDQPDQHGETPSLL', str(hsp.hit.seq)[:40])
        self.assertEqual('PSYSGD*GRRIT*IQEVEAVVGQDQAIALQPGQQERNSVS', str(hsp.query.seq)[-40:])
        self.assertEqual('PSYSG *G+RIT* QE E  + QD AIALQPGQQERNS+S', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('PSYSGG*GQRIT*TQETEVAMSQDCAIALQPGQQERNSIS', str(hsp.hit.seq)[-40:])

    def test_text_2226_tblastx_004(self):
        """Test parsing tblastx output (text_2226_tblastx_004.txt)"""

        blast_file = get_file('text_2226_tblastx_004.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(2, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('random_s00', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(128, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(0, len(qresult))

        # test second qresult
        qresult = qresults[1]
        self.assertEqual('gi|296147483:1-350', qresult.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA,complete cds', qresult.description)
        self.assertEqual(350, qresult.seq_len)
        self.assertEqual('minirefseq_mrna', qresult.target)
        self.assertEqual('tblastx', qresult.program)
        self.assertEqual('2.2.26+', qresult.version)
        self.assertEqual(5, len(qresult))
        # second qresult, first hit
        hit = qresult[0]
        self.assertEqual('gi|296147483|ref|NM_001183135.1|', hit.id)
        self.assertEqual('Saccharomyces cerevisiae S288c Mon2p (MON2) mRNA, complete cds', hit.description)
        self.assertEqual(4911, hit.seq_len)
        self.assertEqual(8, len(hit))
        # second qresult, first hit, first hsp
        hsp = qresult[0].hsps[0]
        self.assertEqual(116, hsp.aln_span)
        self.assertEqual(2e-81, hsp.evalue)
        self.assertEqual(289.0, hsp.bitscore)
        self.assertEqual(626.0, hsp.bitscore_raw)
        self.assertEqual(2, hsp.query_frame)
        self.assertEqual(2, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(1, hsp.hit_start)
        self.assertEqual(349, hsp.query_end)
        self.assertEqual(349, hsp.hit_end)
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINP', str(hsp.query.seq)[:40])
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('WP*TLEGLTPCKGNLKQNCVLYLPNRKEEIQPFAMLVINP', str(hsp.hit.seq)[:40])
        self.assertEqual('WQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', str(hsp.query.seq)[-40:])
        self.assertEqual('WQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('WQCNAYRDCQPFHLFLEAGCLKFWMPSLRLLISRWRFN*K', str(hsp.hit.seq)[-40:])
        # second qresult, first hit, second hsp
        hsp = qresult[0].hsps[1]
        self.assertEqual(116, hsp.aln_span)
        self.assertEqual(5e-78, hsp.evalue)
        self.assertEqual(278.0, hsp.bitscore)
        self.assertEqual(602.0, hsp.bitscore_raw)
        self.assertEqual(-2, hsp.query_frame)
        self.assertEqual(-3, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(116, hsp.ident_num)
        self.assertEqual(116, hsp.pos_num)
        self.assertEqual(-1, hsp.query_strand)
        self.assertEqual(-1, hsp.hit_strand)
        self.assertEqual(1, hsp.query_start)
        self.assertEqual(1, hsp.hit_start)
        self.assertEqual(349, hsp.query_end)
        self.assertEqual(349, hsp.hit_end)
        self.assertEqual('LLIESPSRDE*PQ*RHPKFQTAGFEE*MERLTVPVGIALP', str(hsp.query.seq)[:40])
        self.assertEqual('LLIESPSRDE*PQ*RHPKFQTAGFEE*MERLTVPVGIALP', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('LLIESPSRDE*PQ*RHPKFQTAGFEE*MERLTVPVGIALP', str(hsp.hit.seq)[:40])
        self.assertEqual('WIYH*HGEWLNFFFSIRKIKNAILLQVAFAWSQTLQCSWP', str(hsp.query.seq)[-40:])
        self.assertEqual('WIYH*HGEWLNFFFSIRKIKNAILLQVAFAWSQTLQCSWP', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('WIYH*HGEWLNFFFSIRKIKNAILLQVAFAWSQTLQCSWP', str(hsp.hit.seq)[-40:])
        # second qresult, second hit
        hit = qresult[1]
        self.assertEqual('gi|365982352|ref|XM_003667962.1|', hit.id)
        self.assertEqual('Naumovozyma dairenensis CBS 421 hypothetical protein (NDAI0A06120), mRNA', hit.description)
        self.assertEqual(4932, hit.seq_len)
        self.assertEqual(10, len(hit))
        # second qresult, second hit, first hsp
        hsp = qresult[1].hsps[0]
        self.assertEqual(85, hsp.aln_span)
        self.assertEqual(5e-42, hsp.evalue)
        self.assertEqual(152.0, hsp.bitscore)
        self.assertEqual(327.0, hsp.bitscore_raw)
        self.assertEqual(1, hsp.query_frame)
        self.assertEqual(1, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(62, hsp.ident_num)
        self.assertEqual(73, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(93, hsp.query_start)
        self.assertEqual(87, hsp.hit_start)
        self.assertEqual(348, hsp.query_end)
        self.assertEqual(342, hsp.hit_end)
        self.assertEqual('TIRHASDKSIEILKRVHSFEELERHPDFALPFVLACQSRN', str(hsp.query.seq)[:40])
        self.assertEqual('TI+HASDKSI+ILK + + EEL RHPDF  P VLAC SRN', hsp.aln_annotation['similarity'][:40])
        self.assertEqual('TIKHASDKSIDILKTIQNIEELVRHPDFVTPLVLACSSRN', str(hsp.hit.seq)[:40])
        self.assertEqual('LAMQCLQGLSTVPSIPRSRLSEILDAFIEATHLAMEIQLK', str(hsp.query.seq)[-40:])
        self.assertEqual('+AMQCLQGL++VPSIP SR+ E+LD FIEAT LAMEIQLK', hsp.aln_annotation['similarity'][-40:])
        self.assertEqual('IAMQCLQGLASVPSIPESRIPEVLDGFIEATQLAMEIQLK', str(hsp.hit.seq)[-40:])
        # second qresult, second hit, second hsp
        hsp = qresult[1].hsps[1]
        self.assertEqual(14, hsp.aln_span)
        self.assertEqual(5e-42, hsp.evalue)
        self.assertEqual(26.3, hsp.bitscore)
        self.assertEqual(51.0, hsp.bitscore_raw)
        self.assertEqual(3, hsp.query_frame)
        self.assertEqual(3, hsp.hit_frame)
        self.assertEqual(0, hsp.gap_num)
        self.assertEqual(11, hsp.ident_num)
        self.assertEqual(11, hsp.pos_num)
        self.assertEqual(1, hsp.query_strand)
        self.assertEqual(1, hsp.hit_strand)
        self.assertEqual(68, hsp.query_start)
        self.assertEqual(62, hsp.hit_start)
        self.assertEqual(110, hsp.query_end)
        self.assertEqual(104, hsp.hit_end)
        self.assertEqual('FRIEKKKFNHSPC*', str(hsp.query.seq))
        self.assertEqual('FRI KKKFNH  C*', hsp.aln_annotation['similarity'])
        self.assertEqual('FRI*KKKFNH*TC*', str(hsp.hit.seq))

    def test_text_2230_blastp_001(self):
        """Test parsing blastp output (text_2230_blastp_001.txt)"""

        blast_file = get_file('text_2230_blastp_001.txt')
        qresults = list(parse(blast_file, FMT))
        self.assertEqual(1, len(qresults))
        self.check_common_attrs(qresults)

        # test first qresult
        qresult = qresults[0]
        self.assertEqual('TR11080zzzc0_g2_i2_0', qresult.id)
        self.assertEqual('', qresult.description)
        self.assertEqual(1691, qresult.seq_len)
        self.assertEqual('subject.fasta', qresult.target)
        self.assertEqual('blastp', qresult.program)
        self.assertEqual('2.2.30+', qresult.version)
        self.assertEqual(1, len(qresult))
        self.assertEqual(3, len(qresult.hsps))
        hsp = qresult.hsps[0]
        self.assertIn("PTSP" + ("-" * 79) + "AYSP", hsp.query)
        hsp = qresult.hsps[1]
        self.assertTrue(hsp.query.seq.startswith("AYSPTSPAYSPTSPAYSPTSPAYSPTSPAYS----------PTSPAYSPTSPAYSPTSPA"))
        hsp = qresult.hsps[2]
        self.assertTrue(hsp.query.seq.startswith("YSPTSPAYSPTSPAYSPTSPAYSPTSPAYS----------PTSPAYSPTSPAYSPTSPAY"))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
