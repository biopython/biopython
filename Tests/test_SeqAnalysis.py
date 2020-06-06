# Copyright 2019-2020 by moozeq, KacDom, erpeg and s0lczi. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for SeqAnalysis module."""

import os
import shutil
import unittest

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqAnalysis import database, downloader, analyzer
from Bio.SeqRecord import SeqRecord


def print_success(msg: str):
    print(f'\033[92m[+] {msg}\033[0m')


class TestSeqDatabase(unittest.TestCase):
    def setUp(self) -> None:
        # Declaration of variables needed for running tests.
        self.test_root_dir = 'SeqAnalysis'
        self.filename = f"{self.test_root_dir}/database_tests.txt"  # Additional file made only for reading function test.
        self.dirs = [f'{self.test_root_dir}/Downloads', f'{self.test_root_dir}/Failed', f'{self.test_root_dir}/Database']
        self.test_dir = "test"
        self.test_tuple = (f'{self.test_root_dir}/Downloads/test/A0A1A2.fasta', 'https://www.uniprot.org/uniprot/A0A1A2.fasta')
        self.test_obs = []
        self.test_http = []

    def test_read_list_of_ids(self):
        """
        Check if the function reads the file correctly.
        """
        entries = downloader.read_list_of_ids(self.filename)
        assert entries == ['A0A1A2', 'B0B1B2', 'C0C1C2', 'D0D1D2', 'E0E1E2', 'F0F1F2', 'G0G1G2', 'H0H1H2',
                           'I0I1I2', 'J0J1J2'], \
            'Function execution failure. Check if "test.txt" is in the same directory.'

    def test_mk_dirs(self):
        """
        Check if the function create directories correctly
        """
        for single_dir in self.dirs:
            downloader.mk_dirs(single_dir)
        for single_dir in self.dirs:
            assert os.path.exists(
                single_dir) is True, \
                "Function execution failure. Check if all test directories were deleted before running test."

    def test_mk_subdirs(self):
        """
        Check if the function create directories correctly
        """
        downloader.mk_subdirs(self.test_dir, self.test_root_dir)
        assert os.path.exists(
            f"{self.test_root_dir}/Downloads/test") is True, \
            "Function execution failure. Check if all test directories were deleted before running test."
        assert os.path.exists(
            f"{self.test_root_dir}/Failed/test") is True, \
            "Function execution failure. Check if all test directories were deleted before running test."

    def test_paths_and_urls(self):
        """
        Check if the function create proper tuple of path and url.
        """
        urls = downloader.paths_and_urls(["A0A1A2"], self.test_dir, self.test_root_dir)
        assert urls == [(f'{self.test_root_dir}/Downloads/test/A0A1A2.fasta',
                         'https://www.uniprot.org/uniprot/A0A1A2.fasta')], \
            "Function execution failure."

    def test_fetch_url(self):
        """
        Check if the function returns file path, as well as ,downloads and saves file correctly.
        """
        fetched_url = downloader.fetch_url(self.test_tuple, self.test_obs, self.test_http)
        assert fetched_url == f'{self.test_root_dir}/Downloads/test/A0A1A2.fasta', \
            "Function execution failure."
        assert os.path.exists(f"{self.test_root_dir}/Downloads/test/A0A1A2.fasta") is True, \
            "Function execution failure."

    def test_download_ncbi(self):
        """
        Check if the function downloads and saves NCBI file correctly using Entrez.
        """
        downloader.download_ncbi("YP_025292.1", self.test_dir, self.test_obs, self.test_http, self.test_root_dir)
        assert os.path.exists(f"{self.test_root_dir}/Downloads/test/YP_025292.1.fasta") is True, \
            "Function execution failure."

    def test_downloader_class(self):
        """
        Check if the downloader class downloads and saves files correctly.
        """
        downloader(["A0A1A2", "YP_025292.1"], self.test_dir, self.test_root_dir)
        assert os.path.exists(
            f"{self.test_root_dir}/Downloads/test/A0A1A2.fasta") is True, \
            "Class execution failure. Something wrong with one of the files."
        assert os.path.exists(
            f"{self.test_root_dir}/Downloads/test/YP_025292.1.fasta") is True, \
            "Class execution failure. Something wrong with one of the files."

    def test_get_and_database_class(self):
        """
        Check if the class database creates dict. Also checks if get() function works properly.
        """
        test_database = database(["A0A1A2", "YP_025292.1"], self.test_dir, self.test_root_dir)
        assert type(test_database.database) == dict, \
            "Class execution failure."

        record = list(test_database.get().keys())[0:1]
        same_record = test_database.get(["A0A1A2", "YP_025292.1"])
        assert record != same_record, \
            "Function execution failure. Records are not the same."

    def rm_dirs(self):
        """
        Quick cleanup after certain tests and at the end of testing.
        """
        shutil.rmtree(f"{self.test_root_dir}/Database")
        shutil.rmtree(f"{self.test_root_dir}/Downloads")
        shutil.rmtree(f"{self.test_root_dir}/Failed")


class TestSeqAnalyzer(unittest.TestCase):
    def setUp(self) -> None:
        # Declaring round value - values are rounded to round_value decimal places
        self.round_value = 3

        # Creating test dictionary for testing methods from analyzer class
        self.test_data_dict_of_str = {
            # invalid aa in sequence, should be discarded
            'id1': 'MGKKRAVGAQALACFECERCKSDNEQYXCTNDHVLTMWTPYKDG',
            # to short aa sequence, should be discarded
            'id2': 'MSVKPSKK',
            # to short sequence and invalid aa in sequence, should be discarded
            'id3': 'MSVPKKB',
            # medium length sequence
            'id4': 'MRALAYFGKGNIRFTNHLCKEPHIVAPDELVIDIEWCGICGTDLHEYTDGPIFFPEDGHTHEISHNPLPQAMGHEMAGTVLEVGPGVKNLKVGDKVVVEPTGTCRDYICSYLGLCGAGVQSGGFAERVVMNESHCYKVPDFVPLDVAALIQPLAVCWH',
            # short sequence
            'id5': 'GVKNLKVGDKVVVEPTGVMNESHCYKVPDFVPLDVA',
            # very long sequence
            'id6': 'MGNQSLVVLTESKGEYENCKEPHIVAPDELVIDIMSVQIFGDQVTEERAENARLSAFVGAIAVGDLVKSTLGPKGMDKLLQSASSNTMSLQLLNPKAESLRRDAALKVNVTSAEGLQSVLETNLGPKGTLKMLVDGAGNIKLTKDGKILLTEMQIQSPTAVMIARAAAAQDEITGDGTTTVVCLVGELLKQAYRFIQEGVHPRIITDGFEIARTKALEFLDEYKIEDKVVNDGNVDREFLLQVARSSLSTKVTPELTEVLTPIVTDAVLSVASKDTLNLDLYMVEIMQMQHLSPKDTTFIKGLVLDHGARHPDMPMRVENAHVLILNVSLEYEKTEVNSGFFYSSAEQRDKLAASERKFVDEKLKKIIDLKNEVCGLDNKQGFVIINQKGIDPMSLDVLAKHNILALRRAKRRNMERLQLVTGGEAQNSVDDLSPSILGYSGLVYQETIGEEKFTYVTENKDPKSCTILIKGSTNYALNQTKDAVRDGLRAVANVIKDKCVVPGAGAFFIAASKHLKSSNYSKLGVKGKTKTGIEAFSEALLVVPKTLVKNSGFDPLDVLALCEDELEDAKESEEKRYVGVDLKIGDSCDPTIDGVWDSYRVIRNAINGATGIASNLLLCDELLRAGCMVTNDGATILKSIPLDNPAAKVLVNISKVQDDEVGDGTTSVTVLSAELLREAEKLIDQSKIHPQTIIEGYRLASAAALDALTKAAVDNSHDKTMFREDLIHIAKTTLSSKILSQDKDHFAELATNAILRLKGSTNLEHIQIIKILGGKLSDSFLDEGFILAKKFGNNQPKRIENAKILIANTTLDTDKVKIFGTKFKVDSTAKLAQLEKAEREKMKNKIAKISKFGINTFINRQLIYDYPEQLFTDLGINSIEHADFEGVERLALVTGGEVVSTFDEPSKCKLGECDVIEEIMLGEQPFLKFSGCKAGEACTIVLRGATDQTLDEAERSLHDALSVLSQTTKETRTVLGGGCAEMVMSKAVDTEAQNIDGKKSLAVEAFARALRQLPTILADNAGFDSSELVSKLRSSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKRAVVSSASEAAEVLLRVDNIIRARPRTANRQHMEWCGICGTDLHEYTDGPIFFPEDGHTHEISHNPLPQAMGHEMAGTVLEVGPGVKNLKVGDKVVVEPTGTCRDYICSYLGLCGAGVQSGGFAERVVMNESHCYKVPDFVPLDVAALIQPLAVETELPVKKSSRDNNIGESLTATAFTQSEDEMVDSNQKWQNPNYFKYAWQEYLFIFTCMISQLLNQAGTTQTLSIMNILSDSFGSEGNSKWLMASFPLVSGSFILISGRLGDIYGLKKMLLVGYVLVIIWSLICGITKYSGSDTFFIISRAFQGLGIAFVLPNVLGIIGNIYVGGTFRKNIVISFVGAMAPIGATLGCLFAGLIGTEDPKQWPWAFYAYSIAAFINFVLSIYAIPSTIPTNIHHFSMDWIGSVLGVIGLILLNFVWNQAPISGWNQAYIIVILIISVIFLVVFIIYEIRFAKTPLLPRAVIKDRHMIQIMLALFFGWGSFGIFTFYYFQFQLNIRQYTALWAGGTYFMFLIWGIIAALLVGFTIKNVSPSVFLFFSMVAFNVGSIMASVTPV',
            # empty sequence, should be discarded
            'id7': ''
        }
        # Creating test dictionary for testing method from analyzer class. Contains only proper protein sequences
        self.test_data_dict_of_str_only_proper_sequence = {
            'id4': self.test_data_dict_of_str['id4'],
            'id5': self.test_data_dict_of_str['id5'],
            'id6': self.test_data_dict_of_str['id6']
        }

        # Creating test dictionary for testing analyzer class
        seq1 = {'id1': SeqRecord(Seq(self.test_data_dict_of_str['id1'], IUPAC.protein))}
        seq2 = {'id2': SeqRecord(Seq(self.test_data_dict_of_str['id2'], IUPAC.protein))}
        seq3 = {'id3': SeqRecord(Seq(self.test_data_dict_of_str['id3'], IUPAC.protein))}
        seq4 = {'id4': SeqRecord(Seq(self.test_data_dict_of_str['id4'], IUPAC.protein))}
        seq5 = {'id5': SeqRecord(Seq(self.test_data_dict_of_str['id5'], IUPAC.protein))}
        seq6 = {'id6': SeqRecord(Seq(self.test_data_dict_of_str['id6'], IUPAC.protein))}
        seq7 = {'id7': SeqRecord(Seq(self.test_data_dict_of_str['id7'], IUPAC.protein))}

        # making dict from created SeqRecord objects
        self.test_data_dict_of_seqrecords = {**seq1, **seq2, **seq3, **seq4, **seq5, **seq6, **seq7}

    def test_seqanalyzer_class(self):
        """
        Function that tests the overall functioning of analyzer class and analyzer.results() function

        :return: returns nothing if class passes tests, else function returns text message
        """
        results = analyzer(self.test_data_dict_of_seqrecords)
        assert results.results() == (
        {'median_length': 158.0, 'average_length': 605.667, 'standard_deviation_length': 721.085,
         'id': ['id4', 'id5', 'id6'], 'sequence_length': [158, 36, 1623],
         'molecular_weight': [17245.603, 3884.478, 177765.884], 'theoretical_pi': [5.197, 5.518, 5.556],
         'instability_index': [28.522, 9.625, 33.646], 'flexibility': [0.994, 1.004, 0.998],
         'gravy': [0.042, 0.064, 0.074], 'aromaticity': [0.082, 0.056, 0.076], 'c_term_freq': {'M': 1.0, 'G': 0.5},
         'n_term_freq': {'A': 0.333, 'H': 0.333, 'V': 0.333},
         'secondary_structure_fraction': [{'Helix': 31.646, 'Turn': 23.418, 'Sheet': 23.418},
                                          {'Helix': 36.111, 'Turn': 25.0, 'Sheet': 16.667},
                                          {'Helix': 33.826, 'Turn': 22.058, 'Sheet': 26.802}],
         'extinction_coefficient': [{'Cysteines reduced': 18450, 'Cysteines residues': 18950},
                                    {'Cysteines reduced': 1490, 'Cysteines residues': 1490},
                                    {'Cysteines reduced': 133620, 'Cysteines residues': 134995}], 'amino_acid_comp': [
            {'A': 10, 'C': 8, 'D': 9, 'E': 11, 'F': 6, 'G': 17, 'H': 9, 'I': 9, 'K': 6, 'L': 12, 'M': 4, 'N': 5,
             'P': 11, 'Q': 3, 'R': 4, 'S': 4, 'T': 7, 'V': 16, 'W': 2, 'Y': 5},
            {'A': 1, 'C': 1, 'D': 3, 'E': 2, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4, 'L': 2, 'M': 1, 'N': 2, 'P': 3,
             'Q': 0, 'R': 0, 'S': 1, 'T': 1, 'V': 9, 'W': 0, 'Y': 1},
            {'A': 129, 'C': 23, 'D': 91, 'E': 98, 'F': 72, 'G': 122, 'H': 25, 'I': 131, 'K': 105, 'L': 171, 'M': 37,
             'N': 72, 'P': 51, 'Q': 55, 'R': 56, 'S': 113, 'T': 97, 'V': 123, 'W': 14, 'Y': 38}],
         'amino_acid_percent': [
             {'A': 6.329, 'C': 5.063, 'D': 5.696, 'E': 6.962, 'F': 3.797, 'G': 10.759, 'H': 5.696, 'I': 5.696,
              'K': 3.797, 'L': 7.595, 'M': 2.532, 'N': 3.165, 'P': 6.962, 'Q': 1.899, 'R': 2.532, 'S': 2.532, 'T': 4.43,
              'V': 10.127, 'W': 1.266, 'Y': 3.165},
             {'A': 2.778, 'C': 2.778, 'D': 8.333, 'E': 5.556, 'F': 2.778, 'G': 8.333, 'H': 2.778, 'I': 0.0, 'K': 11.111,
              'L': 5.556, 'M': 2.778, 'N': 5.556, 'P': 8.333, 'Q': 0.0, 'R': 0.0, 'S': 2.778, 'T': 2.778, 'V': 25.0,
              'W': 0.0, 'Y': 2.778},
             {'A': 7.948, 'C': 1.417, 'D': 5.607, 'E': 6.038, 'F': 4.436, 'G': 7.517, 'H': 1.54, 'I': 8.071, 'K': 6.47,
              'L': 10.536, 'M': 2.28, 'N': 4.436, 'P': 3.142, 'Q': 3.389, 'R': 3.45, 'S': 6.962, 'T': 5.977, 'V': 7.579,
              'W': 0.863, 'Y': 2.341}]}, {
            'id4': {'id': ['id4'], 'sequence_length': [158], 'molecular_weight': [17245.603], 'theoretical_pi': [5.197],
                    'instability_index': [28.522], 'flexibility': [0.994], 'gravy': [0.042], 'aromaticity': [0.082],
                    'secondary_structure_fraction': [{'Helix': 31.646, 'Turn': 23.418, 'Sheet': 23.418}],
                    'extinction_coefficient': [{'Cysteines reduced': 18450, 'Cysteines residues': 18950}],
                    'amino_acid_comp': [
                        {'A': 10, 'C': 8, 'D': 9, 'E': 11, 'F': 6, 'G': 17, 'H': 9, 'I': 9, 'K': 6, 'L': 12, 'M': 4,
                         'N': 5, 'P': 11, 'Q': 3, 'R': 4, 'S': 4, 'T': 7, 'V': 16, 'W': 2, 'Y': 5}],
                    'amino_acid_percent': [
                        {'A': 6.329, 'C': 5.063, 'D': 5.696, 'E': 6.962, 'F': 3.797, 'G': 10.759, 'H': 5.696,
                         'I': 5.696, 'K': 3.797, 'L': 7.595, 'M': 2.532, 'N': 3.165, 'P': 6.962, 'Q': 1.899, 'R': 2.532,
                         'S': 2.532, 'T': 4.43, 'V': 10.127, 'W': 1.266, 'Y': 3.165}], 'c_term_freq': {'H': 1},
                    'n_term_freq': {'M': 1}},
            'id5': {'id': ['id5'], 'sequence_length': [36], 'molecular_weight': [3884.478], 'theoretical_pi': [5.518],
                    'instability_index': [9.625], 'flexibility': [1.004], 'gravy': [0.064], 'aromaticity': [0.056],
                    'secondary_structure_fraction': [{'Helix': 36.111, 'Turn': 25.0, 'Sheet': 16.667}],
                    'extinction_coefficient': [{'Cysteines reduced': 1490, 'Cysteines residues': 1490}],
                    'amino_acid_comp': [
                        {'A': 1, 'C': 1, 'D': 3, 'E': 2, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4, 'L': 2, 'M': 1, 'N': 2,
                         'P': 3, 'Q': 0, 'R': 0, 'S': 1, 'T': 1, 'V': 9, 'W': 0, 'Y': 1}], 'amino_acid_percent': [
                    {'A': 2.778, 'C': 2.778, 'D': 8.333, 'E': 5.556, 'F': 2.778, 'G': 8.333, 'H': 2.778, 'I': 0.0,
                     'K': 11.111, 'L': 5.556, 'M': 2.778, 'N': 5.556, 'P': 8.333, 'Q': 0.0, 'R': 0.0, 'S': 2.778,
                     'T': 2.778, 'V': 25.0, 'W': 0.0, 'Y': 2.778}], 'c_term_freq': {'A': 1}, 'n_term_freq': {'G': 1}},
            'id6': {'id': ['id6'], 'sequence_length': [1623], 'molecular_weight': [177765.884],
                    'theoretical_pi': [5.556], 'instability_index': [33.646], 'flexibility': [0.998], 'gravy': [0.074],
                    'aromaticity': [0.076],
                    'secondary_structure_fraction': [{'Helix': 33.826, 'Turn': 22.058, 'Sheet': 26.802}],
                    'extinction_coefficient': [{'Cysteines reduced': 133620, 'Cysteines residues': 134995}],
                    'amino_acid_comp': [
                        {'A': 129, 'C': 23, 'D': 91, 'E': 98, 'F': 72, 'G': 122, 'H': 25, 'I': 131, 'K': 105, 'L': 171,
                         'M': 37, 'N': 72, 'P': 51, 'Q': 55, 'R': 56, 'S': 113, 'T': 97, 'V': 123, 'W': 14, 'Y': 38}],
                    'amino_acid_percent': [
                        {'A': 7.948, 'C': 1.417, 'D': 5.607, 'E': 6.038, 'F': 4.436, 'G': 7.517, 'H': 1.54, 'I': 8.071,
                         'K': 6.47, 'L': 10.536, 'M': 2.28, 'N': 4.436, 'P': 3.142, 'Q': 3.389, 'R': 3.45, 'S': 6.962,
                         'T': 5.977, 'V': 7.579, 'W': 0.863, 'Y': 2.341}], 'c_term_freq': {'V': 1},
                    'n_term_freq': {'M': 1}}}), \
            "class didn't pass the test"

    def test_seqanalyzer_median_standard_deviation_average_length(self):
        """
        Function that tests if the median_standard_deviation_average_length() from analyzer class works properly

        :return: returns nothing if class passes tests, else function returns text message
        """
        results = analyzer.median_standard_deviation_average_length(self.test_data_dict_of_str, self.round_value)
        assert results == (555.45, 268.0, 36.0), \
            "function median_standard_deviation_average_length() didn't pass the test"

    def test_seqanalyzer_main_calculations(self):
        """
        Function that tests if the main_calculations() from analyzer class works properly

        :return: returns nothing if class passes tests, else function returns text message
        """

        results = analyzer.main_calculations(self.test_data_dict_of_str_only_proper_sequence, self.round_value)
        assert results == ([17245.603, 3884.478, 177765.884], [5.197, 5.518, 5.556], [
            {'A': 10, 'C': 8, 'D': 9, 'E': 11, 'F': 6, 'G': 17, 'H': 9, 'I': 9, 'K': 6, 'L': 12, 'M': 4, 'N': 5,
             'P': 11, 'Q': 3, 'R': 4, 'S': 4, 'T': 7, 'V': 16, 'W': 2, 'Y': 5},
            {'A': 1, 'C': 1, 'D': 3, 'E': 2, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4, 'L': 2, 'M': 1, 'N': 2, 'P': 3,
             'Q': 0, 'R': 0, 'S': 1, 'T': 1, 'V': 9, 'W': 0, 'Y': 1},
            {'A': 129, 'C': 23, 'D': 91, 'E': 98, 'F': 72, 'G': 122, 'H': 25, 'I': 131, 'K': 105, 'L': 171, 'M': 37,
             'N': 72, 'P': 51, 'Q': 55, 'R': 56, 'S': 113, 'T': 97, 'V': 123, 'W': 14, 'Y': 38}], [0.082, 0.056, 0.076],
                           [{'A': 6.329, 'C': 5.063, 'D': 5.696, 'E': 6.962, 'F': 3.797, 'G': 10.759, 'H': 5.696,
                             'I': 5.696, 'K': 3.797, 'L': 7.595, 'M': 2.532, 'N': 3.165, 'P': 6.962, 'Q': 1.899,
                             'R': 2.532, 'S': 2.532, 'T': 4.43, 'V': 10.127, 'W': 1.266, 'Y': 3.165},
                            {'A': 2.778, 'C': 2.778, 'D': 8.333, 'E': 5.556, 'F': 2.778, 'G': 8.333, 'H': 2.778,
                             'I': 0.0, 'K': 11.111, 'L': 5.556, 'M': 2.778, 'N': 5.556, 'P': 8.333, 'Q': 0.0, 'R': 0.0,
                             'S': 2.778, 'T': 2.778, 'V': 25.0, 'W': 0.0, 'Y': 2.778},
                            {'A': 7.948, 'C': 1.417, 'D': 5.607, 'E': 6.038, 'F': 4.436, 'G': 7.517, 'H': 1.54,
                             'I': 8.071, 'K': 6.47, 'L': 10.536, 'M': 2.28, 'N': 4.436, 'P': 3.142, 'Q': 3.389,
                             'R': 3.45, 'S': 6.962, 'T': 5.977, 'V': 7.579, 'W': 0.863, 'Y': 2.341}],
                           [28.522, 9.625, 33.646], [0.994, 1.004, 0.998],
                           [{'Helix': 31.646, 'Turn': 23.418, 'Sheet': 23.418},
                            {'Helix': 36.111, 'Turn': 25.0, 'Sheet': 16.667},
                            {'Helix': 33.826, 'Turn': 22.058, 'Sheet': 26.802}], [0.042, 0.064, 0.074],
                           [{'Cysteines reduced': 18450, 'Cysteines residues': 18950},
                            {'Cysteines reduced': 1490, 'Cysteines residues': 1490},
                            {'Cysteines reduced': 133620, 'Cysteines residues': 134995}]), \
            "function main_calculations() didn't pass the test"

    def test_seqanalyzer_terminal_aa_counter(self):
        """
        Function that tests if the terminal_aa_counter() from analyzer class works properly

        :return: returns nothing if class passes tests, else function returns text message
        """
        results = analyzer.terminal_aa_counter(self.test_data_dict_of_str_only_proper_sequence, self.round_value)
        assert results == ({'A': 0.333, 'H': 0.333, 'V': 0.333}, {'M': 1.0, 'G': 0.5}), \
            "function test_analyzer_terminal_aa_counter() didn't pass the test"


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
