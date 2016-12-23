# Copyright 2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Simple 'print-and-compare' unit test for fasta-m10 parser.

Created to check for any regressions from my new implementation of the
parser.
"""
from __future__ import print_function

import os
from Bio import AlignIO
from Bio._py3k import StringIO
from Bio.AlignIO import FastaIO
import unittest

# test_files is a list of tuples containing:
# - string:  file format
# - integer: number of sequences per alignment
# - integer: number of alignments
# - string:  relative filename
#
# Most of the input files are also used by test_SeqIO.py,
# and by other additional tests as noted below.
test_files = [
    ("fasta-m10", 2, 4, 'Fasta/output001.m10'),
    ("fasta-m10", 2, 6, 'Fasta/output002.m10'),
    ("fasta-m10", 2, 3, 'Fasta/output003.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output004.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output005.m10'),
    ("fasta-m10", 2, 1, 'Fasta/output006.m10'),
    ("fasta-m10", 2, 9, 'Fasta/output007.m10'),
    ("fasta-m10", 2, 12, 'Fasta/output008.m10'),
    ]

# Main tests...
for (t_format, t_per, t_count, t_filename) in test_files:
    assert t_format == "fasta-m10" and t_per == 2

    print("Testing reading %s format file %s with %i alignments"
          % (t_format, t_filename, t_count))
    assert os.path.isfile(t_filename), t_filename

    # Try as an iterator using handle
    with open(t_filename, "r") as handle:
        alignments = list(AlignIO.parse(handle, format=t_format))
    assert len(alignments) == t_count, \
         "Found %i alignments but expected %i" % (len(alignments), t_count)
    for alignment in alignments:
        assert len(alignment) == t_per, \
            "Expected %i records per alignment, got %i" \
            % (t_per, len(alignment))

    # Print the alignment
    for i, alignment in enumerate(alignments):
        print("=" * 78)
        print("Alignment %i, with %i sequences of length %i"
              % (i,
                 len(alignment),
                 alignment.get_alignment_length()))
        for k in sorted(alignment._annotations):
            print(" - %s: %r" % (k, alignment._annotations[k]))
        assert alignment[0].name == "query"
        assert alignment[1].name == "match"
        # Show each sequence row horizontally
        for record in alignment:
            print("-" * 78)
            print(record.id)
            print(record.description)
            print(repr(record.seq))
            assert not record.features
            assert not record.letter_annotations
            for k in sorted(record.annotations):
                print(" - %s: %r" % (k, record.annotations[k]))
    print("=" * 78)
print("Finished tested reading files")

class TestSelf(unittest.TestCase):
    def test_example(self):
        simple_example = \
    """# /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
    FASTA searches a protein or DNA sequence data bank
     version 34.26 January 12, 2007
    Please cite:
      W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448
        
    Query library NC_002127.faa vs NC_009649.faa library
    searching NC_009649.faa library
      1>>>gi|10955263|ref|NP_052604.1| plasmid mobilization [Escherichia coli O157:H7 s 107 aa - 107 aa
     vs  NC_009649.faa library
        
      45119 residues in   180 sequences
      Expectation_n fit: rho(ln(x))= 6.9146+/-0.0249; mu= -5.7948+/- 1.273
     mean_var=53.6859+/-13.609, 0's: 0 Z-trim: 1  B-trim: 9 in 1/25
     Lambda= 0.175043
         
    FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
     join: 36, opt: 24, open/ext: -10/-2, width:  16
     Scan time:  0.000
    The best scores are:                                      opt bits E(180)
    gi|152973457|ref|YP_001338508.1| ATPase with chape ( 931)   71 24.9    0.58
    gi|152973588|ref|YP_001338639.1| F pilus assembly  ( 459)   63 23.1    0.99
        
    >>>gi|10955263|ref|NP_052604.1|, 107 aa vs NC_009649.faa library
    ; pg_name: /opt/fasta/fasta34
    ; pg_ver: 34.26
    ; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
    ; pg_name: FASTA
    ; pg_ver: 3.5 Sept 2006
    ; pg_matrix: BL50 (15:-5)
    ; pg_open-ext: -10 -2
    ; pg_ktup: 2
    ; pg_optcut: 24
    ; pg_cgap: 36
    ; mp_extrap: 60000 180
    ; mp_stats:  Expectation_n fit: rho(ln(x))= 6.9146+/-0.0249; mu= -5.7948+/- 1.273  mean_var=53.6859+/-13.609, 0's: 0 Z-trim: 1  B-trim: 9 in 1/25  Lambda= 0.175043
    ; mp_KS: -0.0000 (N=0) at 8159228
    >>gi|152973457|ref|YP_001338508.1| ATPase with chaperone activity, ATP-binding subunit [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
    ; fa_frame: f
    ; fa_initn:  65
    ; fa_init1:  43
    ; fa_opt:  71
    ; fa_z-score: 90.3
    ; fa_bits: 24.9
    ; fa_expect:   0.58
    ; sw_score: 71
    ; sw_ident: 0.250
    ; sw_sim: 0.574
    ; sw_overlap: 108
    >gi|10955263| ..
    ; sq_len: 107
    ; sq_offset: 1
    ; sq_type: p
    ; al_start: 5
    ; al_stop: 103
    ; al_display_start: 1
    --------------------------MTKRSGSNT-RRRAISRPVRLTAE
    ED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----
    RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSRLMAD
    >gi|152973457|ref|YP_001338508.1| ..
    ; sq_len: 931
    ; sq_type: p
    ; al_start: 96
    ; al_stop: 195
    ; al_display_start: 66
    SDFFRIGDDATPVAADTDDVVDASFGEPAAAGSGAPRRRGSGLASRISEQ
    SEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKI
    KVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSRAFVASNELGHSY
    VGPEHFLIGLAEEGEGLAANLLRRYGLTPQ
    >>gi|152973588|ref|YP_001338639.1| F pilus assembly protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
    ; fa_frame: f
    ; fa_initn:  33
    ; fa_init1:  33
    ; fa_opt:  63
    ; fa_z-score: 86.1
    ; fa_bits: 23.1
    ; fa_expect:   0.99
    ; sw_score: 63
    ; sw_ident: 0.266
    ; sw_sim: 0.656
    ; sw_overlap: 64
    >gi|10955263| ..
    ; sq_len: 107
    ; sq_offset: 1
    ; sq_type: p
    ; al_start: 32
    ; al_stop: 94
    ; al_display_start: 2
    TKRSGSNTRRRAISRPVRLTAEEDQEIRKRAAECGKTVSGFLRAAALGKK
    VNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAITEYHRAL
    LSRLMAD
    >gi|152973588|ref|YP_001338639.1| ..
    ; sq_len: 459
    ; sq_type: p
    ; al_start: 191
    ; al_stop: 248
    ; al_display_start: 161
    VGGLFPRTQVAQQKVCQDIAGESNIFSDWAASRQGCTVGG--KMDSVQDK
    ASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLTGTLIFG
    ENSEITPLPARTTDQDLIRAMMEGGTAKIYHCNDSDKCLKVVADATVTIT
    SNKALKSQISALLSSIQNKAVADEKLTDQE
      2>>>gi|10955264|ref|NP_052605.1| hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa - 126 aa
     vs  NC_009649.faa library
          
      45119 residues in   180 sequences
      Expectation_n fit: rho(ln(x))= 7.1374+/-0.0246; mu= -7.6540+/- 1.313
    mean_var=51.1189+/-13.171, 0's: 0 Z-trim: 1  B-trim: 8 in 1/25
    Lambda= 0.179384

    FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
     join: 36, opt: 24, open/ext: -10/-2, width:  16
     Scan time:  0.000
    The best scores are:                                      opt bits E(180)
    gi|152973462|ref|YP_001338513.1| hypothetical prot ( 101)   58 22.9    0.29
        
    >>>gi|10955264|ref|NP_052605.1|, 126 aa vs NC_009649.faa library
    ; pg_name: /opt/fasta/fasta34
    ; pg_ver: 34.26
    ; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
    ; pg_name: FASTA
    ; pg_ver: 3.5 Sept 2006
    ; pg_matrix: BL50 (15:-5)
    ; pg_open-ext: -10 -2
    ; pg_ktup: 2
    ; pg_optcut: 24
    ; pg_cgap: 36
    ; mp_extrap: 60000 180
    ; mp_stats:  Expectation_n fit: rho(ln(x))= 7.1374+/-0.0246; mu= -7.6540+/- 1.313  mean_var=51.1189+/-13.171, 0's: 0 Z-trim: 1  B-trim: 8 in 1/25  Lambda= 0.179384
    ; mp_KS: -0.0000 (N=0) at 8159228
    >>gi|152973462|ref|YP_001338513.1| hypothetical protein KPN_pKPN3p05904 [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
    ; fa_frame: f
    ; fa_initn:  50
    ; fa_init1:  50
    ; fa_opt:  58
    ; fa_z-score: 95.8
    ; fa_bits: 22.9
    ; fa_expect:   0.29
    ; sw_score: 58
    ; sw_ident: 0.289
    ; sw_sim: 0.632
    ; sw_overlap: 38
    >gi|10955264| ..
    ; sq_len: 126
    ; sq_offset: 1
    ; sq_type: p
    ; al_start: 1
    ; al_stop: 38
    ; al_display_start: 1
    ------------------------------MKKDKKYQIEAIKNKDKTLF
    IVYATDIYSPSEFFSKIESDLKKKKSKGDVFFDLIIPNGGKKDRYVYTSF
    NGEKFSSYTLNKVTKTDEYN
    >gi|152973462|ref|YP_001338513.1| ..
    ; sq_len: 101
    ; sq_type: p
    ; al_start: 44
    ; al_stop: 81
    ; al_display_start: 14
    DALLGEIQRLRKQVHQLQLERDILTKANELIKKDLGVSFLKLKNREKTLI
    VDALKKKYPVAELLSVLQLARSCYFYQNVCTISMRKYA
      3>>>gi|10955265|ref|NP_052606.1| hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa - 346 aa
     vs  NC_009649.faa library
          
      45119 residues in   180 sequences
      Expectation_n fit: rho(ln(x))= 6.0276+/-0.0276; mu= 3.0670+/- 1.461
     mean_var=37.1634+/- 8.980, 0's: 0 Z-trim: 1  B-trim: 14 in 1/25
     Lambda= 0.210386

    FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
    join: 37, opt: 25, open/ext: -10/-2, width:  16
    Scan time:  0.020
    The best scores are:                                      opt bits E(180)
    gi|152973545|ref|YP_001338596.1| putative plasmid  ( 242)   70 27.5   0.082
    >>>gi|10955265|ref|NP_052606.1|, 346 aa vs NC_009649.faa library
    ; pg_name: /opt/fasta/fasta34
    ; pg_ver: 34.26
    ; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
    ; pg_name: FASTA
    ; pg_ver: 3.5 Sept 2006
    ; pg_matrix: BL50 (15:-5)
    ; pg_open-ext: -10 -2
    ; pg_ktup: 2
    ; pg_optcut: 25
    ; pg_cgap: 37
    ; mp_extrap: 60000 180
    ; mp_stats:  Expectation_n fit: rho(ln(x))= 6.0276+/-0.0276; mu= 3.0670+/- 1.461  mean_var=37.1634+/- 8.980, 0's: 0 Z-trim: 1  B-trim: 14 in 1/25  Lambda= 0.210386
    ; mp_KS: -0.0000 (N=0) at 8159228
    >>gi|152973545|ref|YP_001338596.1| putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
    ; fa_frame: f
    ; fa_initn:  52
    ; fa_init1:  52
    ; fa_opt:  70
    ; fa_z-score: 105.5
    ; fa_bits: 27.5
    ; fa_expect:  0.082
    ; sw_score: 70
    ; sw_ident: 0.279
    ; sw_sim: 0.651
    ; sw_overlap: 43
    >gi|10955265| ..
    ; sq_len: 346
    ; sq_offset: 1
    ; sq_type: p
    ; al_start: 197
    ; al_stop: 238
    ; al_display_start: 167
    DFMCSILNMKEIVEQKNKEFNVDIKKETIESELHSKLPKSIDKIHEDIKK
    QLSC-SLIMKKIDVEMEDYSTYCFSALRAIEGFIYQILNDVCNPSSSKNL
    GEYFTENKPKYIIREIHQET
    >gi|152973545|ref|YP_001338596.1| ..
    ; sq_len: 242
    ; sq_type: p
    ; al_start: 52
    ; al_stop: 94
    ; al_display_start: 22
    IMTVEEARQRGARLPSMPHVRTFLRLLTGCSRINSDVARRIPGIHRDPKD
    RLSSLKQVEEALDMLISSHGEYCPLPLTMDVQAENFPEVLHTRTVRRLKR
    QDFAFTRKMRREARQVEQSW
    >>><<<

    579 residues in 3 query   sequences
    45119 residues in 180 library sequences
    Scomplib [34.26]
     start: Tue May 20 16:38:45 2008 done: Tue May 20 16:38:45 2008
     Total Scan time:  0.020 Total Display time:  0.010
        
    Function used was FASTA [version 34.26 January 12, 2007]
    """

    
        alignments = list(FastaIO.FastaM10Iterator(StringIO(simple_example)))
        assert len(alignments) == 4, len(alignments)
        assert len(alignments[0]) == 2
        for a in alignments:
            self.assertEqual((len(a), a.get_alignment_length()))
            for r in a:
                print("%s %s %i" % (r.seq, r.id, r.annotations["original_length"]))
                # print(a.annotations)
        print("Done")
            
        path = "../../Tests/Fasta/"
        files = sorted(f for f in os.listdir(path) if os.path.splitext(f)[-1] == ".m10")
        for filename in files:
            if os.path.splitext(filename)[-1] == ".m10":
                print("")
                print(filename)
                print("=" * len(filename))
                for i, a in enumerate(FastaIO.FastaM10Iterator(open(os.path.join(path, filename)))):
                    print("#%i, %s" % (i + 1, a))
                    for r in a:
                        if "-" in r.seq:
                            self.assertEqual(r.seq.alphabet.gap_char, "-")
                        else:
                            assert not hasattr(r.seq.alphabet, "gap_char")
