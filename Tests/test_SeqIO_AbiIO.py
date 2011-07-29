# Copyright 2011 by Wibowo Arindrarto (w.arindrarto@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import unittest

from os.path import join, basename

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio._py3k import _bytes_to_string, _as_bytes

test_data = {
'data_empty': { 
              'path': ['Abi', 'empty.ab1'],
              'seq': 'NNNNN',
              'qual': [0, 0, 0, 0, 0],
              'sample': '226041_C-ME-19_pCAGseqF',
              'sample_well': 'C10',
              'machine_model': '3730',
              'run_start': '2009-12-12 09:56:53',
              'run_finish': '2009-12-12 11:44:49',
              },       
'data_3730': { 
             'path': ['Abi', '3730.ab1'],
             'seq': 
'GGGCGAGCKYYAYATTTTGGCAAGAATTGAGCTCTATGGCCACAACCATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGGGCGGCAGCGGCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAARACCCGCGCCGAGGTGAARTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAAGGGGCAYCGCACCTTTC',
             'qual': 
[20, 3, 4, 4, 4, 6, 4, 4, 0, 0, 0, 6, 0, 10, 20, 26, 22, 17, 21, 31, 29, 32, 28, 18, 23, 17, 19, 35, 36, 50, 39, 50, 50, 50, 50, 50, 25, 35, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 35, 39, 33, 20, 35, 31, 50, 50, 50, 50, 50, 50, 50, 50, 50, 31, 50, 35, 31, 23, 28, 31, 21, 43, 39, 35, 24, 30, 26, 35, 31, 50, 50, 50, 50, 50, 50, 50, 50, 50, 39, 31, 24, 39, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 31, 31, 43, 43, 50, 50, 50, 50, 50, 31, 31, 31, 31, 50, 50, 50, 50, 50, 50, 50, 50, 31, 31, 35, 50, 50, 50, 50, 31, 36, 55, 55, 55, 55, 36, 55, 55, 55, 55, 55, 36, 55, 55, 55, 55, 55, 36, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 40, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 36, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 40, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 43, 43, 50, 43, 43, 50, 43, 43, 50, 43, 43, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 43, 43, 50, 43, 43, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 28, 28, 35, 28, 28, 35, 28, 28, 35, 28, 28, 35, 28, 28, 35, 28, 21, 28, 35, 28, 28, 35, 35, 35, 35, 35, 37, 38, 21, 28, 35, 28, 28, 35, 35, 35, 35, 35, 35, 35, 36, 36, 21, 39, 35, 35, 35, 39, 35, 37, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 28, 28, 35, 35, 28, 28, 35, 35, 35, 36, 36, 22, 39, 35, 35, 35, 35, 35, 35, 37, 38, 28, 35, 21, 36, 36, 37, 35, 35, 20, 39, 39, 35, 35, 35, 35, 37, 38, 28, 35, 37, 34, 35, 24, 24, 27, 25, 20, 24, 37, 35, 27, 21, 20, 21, 27, 17, 20, 24, 32, 26, 20, 12, 20, 10, 20, 24, 25, 23, 20, 32, 24, 24, 23, 20, 24, 23, 18, 34, 34, 34, 22, 26, 24, 24, 18, 22, 22, 23, 25, 20, 12, 20, 24, 23, 24, 23, 22, 20, 20, 0, 20, 24, 23, 20, 8, 10, 4, 20, 20, 3, 7, 19, 20, 4, 4, 7, 7, 0, 7, 11, 18, 8, 3, 23, 23, 20, 11, 4, 20, 18, 12, 20, 20, 20, 4, 20, 4, 2, 3, 21, 21, 21, 21, 10, 15, 14, 15, 19, 2, 4, 3, 6, 11, 3, 4, 6, 21, 16, 20, 11, 1, 4, 12, 0, 15, 8, 1, 3, 3, 12, 1, 11, 20, 4],
             'sample': '226032_C-ME-18_pCAGseqF',
             'sample_well': 'B9',
             'dye': 'Z-BigDyeV3',
             'polymer': 'POP7                            ',
             'machine_model': '3730',
             'run_start': '2009-12-12 09:56:53',
             'run_finish': '2009-12-12 11:44:49',
            },
'data_3100': { 
             'path': ['Abi', '3100.ab1'],
             'seq': 
'CAAGATTGCATTCATGATCTACGATTACTAGCGATTCCAGCTTCATATAGTCGAGTTGCAGACTACAATCCGAACTGAGAACAACTTTATGGGATTTGCTTGACCTCGCGGTTTCGCTGCCCTTTGTATTGTCCATTGTAGCACGTGTGTAGCCCAAATCATAAGGGGCATGATGATTTGACGTCATCCCCACCTTCCTCCGGTTTGTCACCGGCAGTCAACTTAGAGTGCCCAACTTAATGATGGCAACTAAGCTTAAGGGTTGCGCTCGTTGCGGGACTTAACCCAACATCTCACGACACGAGCTGACGACAACCATGCACCACCTGTCACTCTGTCCCCCGAAGGGGAAAACTCTATCTCTAGAGGAGTCAGAGGATGTCAAGATTTGGTAAGGTTCTTCGCGTTGCTTCGAATTAAACCACATGCTCCACCGCTTGTGCGGGTCCCCGTCAATTCCTTTGAGTTTCAACCTTGCGGTCGTACTCCCCAGGCGGAGTGCTTAATGCGTTAGCTGCAGCACTAAGGGGCGGAAACCCCCTAACACTTAGCACTCATCGTTTACGGCGTGGACTACCAGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGCACATCAGCGTCAGTTACAGACCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCATATCTCTGCGCATTTCACCGCTACACATGGAAATTCCACTTTCCTCTTCTGCACTCAAGTTTTCCAGATTTCGATGAACCTTCAACGATGGAGGCCCGTGGCTTTTCACCATTCAAGGAACCTTTTA',
             'qual': 
[5, 3, 4, 4, 4, 5, 9, 4, 4, 4, 5, 4, 4, 4, 4, 4, 6, 13, 23, 20, 15, 17, 10, 9, 13, 23, 18, 14, 19, 34, 15, 24, 34, 18, 30, 35, 41, 34, 22, 39, 53, 58, 35, 46, 24, 42, 46, 46, 54, 42, 53, 58, 52, 41, 22, 62, 52, 62, 52, 62, 41, 38, 16, 42, 62, 56, 62, 62, 62, 62, 62, 62, 54, 56, 39, 62, 49, 35, 19, 44, 44, 42, 62, 62, 62, 62, 54, 62, 62, 62, 62, 49, 49, 43, 39, 62, 62, 62, 62, 62, 62, 62, 59, 43, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 49, 59, 62, 62, 59, 59, 62, 62, 62, 62, 62, 59, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 59, 62, 62, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 59, 59, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 50, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 54, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 59, 59, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 49, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 42, 31, 31, 26, 18, 24, 33, 56, 56, 62, 56, 62, 47, 39, 62, 62, 62, 62, 62, 62, 49, 47, 62, 62, 62, 49, 62, 62, 62, 62, 49, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 58, 58, 58, 58, 58, 58, 58, 58, 58, 54, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 54, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 53, 58, 58, 58, 51, 49, 58, 58, 58, 58, 50, 58, 58, 58, 58, 58, 58, 51, 32, 58, 58, 42, 58, 51, 58, 49, 39, 58, 58, 51, 58, 53, 58, 50, 58, 58, 58, 53, 54, 53, 58, 54, 53, 58, 53, 58, 32, 58, 44, 46, 50, 58, 54, 58, 58, 49, 41, 58, 58, 58, 54, 54, 54, 51, 52, 53, 51, 54, 54, 54, 38, 52, 53, 47, 54, 54, 53, 54, 54, 54, 52, 54, 42, 47, 35, 42, 54, 45, 54, 53, 54, 23, 34, 24, 51, 52, 33, 46, 44, 33, 49, 51, 45, 53, 32, 52, 53, 52, 51, 36, 52, 53, 52, 52, 37, 53, 52, 40, 27, 24, 46, 16, 33, 45, 51, 43, 51, 51, 46, 46, 45, 42, 31, 43, 41, 43, 35, 43, 43, 42, 31, 46, 46, 32, 46, 35, 42, 46, 42, 31, 31, 43, 42, 35, 21, 42, 35, 28, 46, 47, 37, 41, 29, 39, 31, 32, 39, 31, 41, 43, 34, 43, 36, 27, 31, 27, 35, 36, 24, 44, 33, 31, 32, 38, 35, 31, 32, 29, 31, 31, 31, 38, 38, 36, 31, 34, 30, 21, 28, 32, 32, 47, 31, 33, 30, 23, 22, 25, 25, 26, 30, 30, 21, 28, 31, 21, 14, 5, 3, 16, 9, 18, 27, 8, 19, 4, 5, 13, 47, 27, 14, 21, 9, 17, 21, 17, 27, 14, 23, 23, 14, 15, 28, 15, 12, 8, 5, 3, 17, 28, 16, 13, 17, 11, 2, 2, 3, 29, 7, 3, 1, 3, 3, 8, 7, 7, 24, 12, 3, 4, 15, 3, 4, 7, 5, 1, 2, 2, 8, 24, 6, 4, 11, 23, 9, 6, 4, 3, 5, 3, 5, 17, 17, 22, 11, 4, 5, 8, 16, 10, 6, 11, 6, 11, 5, 5, 27, 18, 7, 15, 16, 30, 21, 14, 7],
             'sample': '16S_S2_1387R',
             'sample_well': 'H3',
             'dye': 'Z-BigDyeV3',
             'polymer': 'POP7                            ',
             'machine_model': '3100',
             'run_start': '2010-01-27 09:52:45',
             'run_finish': '2010-01-27 10:41:07',
            },
'data_310': { 
            'path': ['Abi', '310.ab1'],
            'seq': 
'TGATNTTNACNNTTTTGAANCANTGAGTTAATAGCAATNCTTTACNAATAAGAATATACACTTTCTGCTTAGGGATGATAATTGGCAGGCAAGTGAATCCCTGAGCGTGNATTTGATAATGACCTAAATAATGGATGGGGTTTTAATTCCCAGACCTTCCCCTTTTTAANNGGNGGATTANTGGGGGNNNAACNNGGGGGGCCCTTNCCNAAGGGGGAAAAAATTTNAAACCCCCCNAGGNNGGGNAAAAAAAAATTTCCAAATTNCCGGGGTNNCCCCCAANTTTTTNCCGCNGGGAAAANNNNCCCCCCCNGGGNCCCCCCCCNNAAAAAAAAAAAAAAAAACCCCCCCCCCNTTGGGGNGGTNTNCNCCCCCNNANAANNGGGGGNNAAAAAAAAAGGCCCCCCCCAAAAAAAACCCNCNTTCTNNCNNNNNGNNCNGNNCCCCCNNCCNTNTNGGGGGGGGGGGNGGAAAAAAAACCCCTTTNTGNNNANANNAACCCNCTCNTNTTTTTTTTTTTANGNNNNCNNNNCAAAAAAAAANCNCCCCCNNCNNNCNNNCNCCCCNNNNTNAAAANANNAANNNNTTTTTTTNGGGGGGGTGNGCGNCCCNNANCNNNNNNNNGCGNGGNCNCCNNCCCNCNANAAANNNTNTTTTTTTTTTTTTTTNTNNTCNNCCCNNNCCCCNNCCCCCCCCCCCCCNCCNCNNNNNGGGGNNNCGGNNCNNNNNNNCCNTNCTNNANATNCCNTTNNNNNNNNGNNNNNNNNACNNNNNTNNTNNNCNNNNNNNNNNNNNNCNNNNNNCNNCCCNNCANNNNNNNCNNNNNNNNNNNNNNNNNNNNNTCNCTNCNCNCCCCNCCCNNNNNNNG',
            'qual': 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            'sample': 'D11F',
            'sample_well': 'C5',
            'dye': None,
            'polymer': None,
            'machine_model': '310 ',
            'run_start': '2009-02-19 01:19:30',
            'run_finish': '2009-02-19 04:04:15',
           },
}

test_data_fake = {
'fake': {
              'path': ['Abi', 'fake.ab1'],
             }
}

def open_files(test_array):
    for trace in test_array:
        test_array[trace]['handle'] = open(join(*test_array[trace]['path']), 'rb')

def open_files_wrong_mode(test_array):
    for trace in test_array:
        test_array[trace]['handle'] = open(join(*test_array[trace]['path']))

def close_files(test_array):
    for trace in test_array:
        test_array[trace]['handle'].close()

class TestAbi(unittest.TestCase):
    
    def setUp(self):
        open_files(test_data)

    def tearDown(self):
        close_files(test_data)

    def test_file_type(self):
        """Test if filetype is ABIF."""
        for trace in test_data:
            self.assertEqual(test_data[trace]['handle'].read(4), _as_bytes('ABIF'))

    def test_seqrecord(self):
        """Test if the extracted seqrecords data are equal to expected values."""
        for trace in test_data:
            record = SeqIO.read(test_data[trace]['handle'], 'abi')
            self.assertEqual(basename(test_data[trace]["path"][-1]).replace('.ab1',''), record.name)
            self.assertEqual(test_data[trace]['seq'], str(record.seq))
            self.assertEqual(test_data[trace]['qual'], record.letter_annotations['phred_quality'])
            self.assertEqual(test_data[trace]['sample'], record.id)
            self.assertEqual(test_data[trace]['sample_well'], record.annotations['sample_well'])
            self.assertEqual(test_data[trace]['machine_model'], record.annotations['machine_model'])
            self.assertEqual(test_data[trace]['run_start'], record.annotations['run_start'])
            self.assertEqual(test_data[trace]['run_finish'], record.annotations['run_finish'])

    def test_trim(self):
        """Test if trim works."""
        for trace in test_data:
            record = SeqIO.read(test_data[trace]['handle'], 'abi-trim')
            if trace != 'data_empty':
                self.assertNotEqual(str(record.seq), test_data[trace]['seq'])
                self.assertTrue(str(record.seq) in test_data[trace]['seq'])
            else:
                self.assertEqual(str(record.seq), test_data[trace]['seq'])


class TestAbiWrongMode(unittest.TestCase):
    
    def test_file_mode(self):
        """Test if exception is raised if file is not opened in 'rb' mode."""
        open_files_wrong_mode(test_data)
        for trace in test_data:
            self.assertRaises(ValueError, SeqIO.read, test_data[trace]['handle'], 'abi')
        close_files(test_data)


class TestAbiFake(unittest.TestCase):

    def test_file_type(self):
        """Test if error is raised if filetype is not ABIF."""
        open_files(test_data_fake)
        for trace in test_data_fake:
            self.assertRaises(IOError, SeqIO.read, test_data_fake[trace]['handle'], 'abi')
        close_files(test_data_fake)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
