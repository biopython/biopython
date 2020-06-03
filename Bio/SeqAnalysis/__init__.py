# Copyright 2019-2020 by moozeq, KacDom, erpeg and s0lczi. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

r"""Proteins sequences analysis and visualization.

Bio.SeqAnalysis provides python components for storing, analyzing and visualizing proteins sequences.
Components can be used separately or through API as fully working, automated pipeline.

Components
==========
There are 3 main components which basically are 3 classes:
    **SeqDatabase** -- downloading, storing and providing protein sequences
    **SeqAnalyzer** -- complex analysis of protein sequences
    **SeqVisualizer** -- protein sequences analysis report visualization

SeqDatabase
-----------
Component responsible for downloading, storing and providing sequences.
Proteins sequences are then concurrently downloaded, which is much faster than
downloading them one by one. After that, they are converted to proper format and
stored on disk. Report about which sequences were properly downloaded is printed.
Then user can use SeqDatabase object to retrieve and manage those sequences.
Sequences are available in SeqIO format.

:Input:

Creating SeqDatabase object is done by providing list of UniProt/NCBI IDs
which sequences should be downloaded and converted to :class:`SeqRecord <Bio.SeqRecord.SeqRecord>` format:

>>> from Bio import SeqAnalysis
>>> seq_db = SeqAnalysis.database(
...     [
...         'YP_025292.1',
...         '1JOY',
...         '8OO85XXX'
...     ]
... )
[*] Downloading sequences...
[+] Properly downloaded and stored sequences (2/3):
    YP_025292.1
    1JOY
[-] Unable to download and store sequences (1/3):
    8OO85XXX


:Output:

After creating database, retrieving sequences in SeqRecord format is available in two ways:
    - providing sequence ID or list of IDs
    - calling `get()` without parameters and use index

Retrieving single sequence:

>>> record = seq_db.get('YP_025292.1')
>>> same_record = seq_db.get()[0]
>>> record == same_record
True
>>> print(record)
ID: YP_025292.1
Name: HokC
Description: toxic membrane protein
Number of features: 0
Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())


Retrieving multiple sequences:

>>> records = seq_db.get(['YP_025292.1', '1JOY'])
>>> same_records = seq_db.get()[0:2]
>>> records == same_records
True
>>> print(records[1])
ID: 1JOY
Name: EnvZ
Description: Homodimeric domain of EnvZ from E. coli
Number of features: 1
Per letter annotation for: secondary_structure
Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR', IUPACProtein())


SeqAnalyzer
-----------
Component responsible for complex analysis of received protein sequences. Analyzed are:
    - sequences lengths
    - molecular weights
    - theoretical isoelectric point
    - amino acid composition
    - atomic composition
    - extinction coefficient
    - instability index
    - average hydrophobicity and hydrophilicity (GRAVY)
    - estimated half-life
    - amino acids frequencies on C/N terminals

:Results as dict example:

    .. code-block:: python

        data = {
            'id': ['AsDF', 'asd'],
            'sequence_length': [12, 14],
            'molecular_weight': [44, 123],
            'theoretical_pl': [1.4, 24.5],
            'amino_acid_comp': [{'A': 12, 'G': 24}, {'L': 124, 'A': 52}],
            'atomic_comp': [{'C': 42, 'H': 335, 'N': 55}, {'C': 55, 'H': 523, 'N': 13, 'O': 53}],
            'extinction_coefficient': [12.4, 43.5],
            'instability_index': [1.2, 0.4],
            'gravy': [23.5, 53.5],
            'est_half_life': [244.5, 3333.3],
            'c_term_freq': [{'A': 0.2, 'G': 0.8}, {'A': 0.3, 'L': 0.7}],
            'n_term_freq': [{'A': 0.8, 'G': 0.2}, {'A': 0.7, 'L': 0.3}],
        }


:Input:

Object of SeqAnalyzer class is created by providing list of SeqRecord objects
which will be analyzed. They can be passed from SeqDatabase:

>>> from Bio import SeqAnalysis
>>> seq_db = SeqAnalysis.database(
...     [
...         'YP_025292.1',
...         '1JOY'
...     ]
... )
>>> records = seq_db.get()
>>> seq_analyzer = SeqAnalysis.analyzer(records)


:Output:

Results are available in pandas `DataFrame <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html>`
format or as simple dict. Analysis are per protein sequence but additionally contains medians,
averages and deviations for all sequences.

>>> results = seq_analyzer.results()
>>> print(results)
            id  sequence_length  ...           c_term_freq           n_term_freq
0  YP_025292.1               12  ...  {'A': 0.2, 'G': 0.8}  {'A': 0.8, 'G': 0.2}
1         1JOY               14  ...  {'A': 0.3, 'L': 0.7}  {'A': 0.7, 'L': 0.3}

[2 rows x 12 columns]


SeqVisualizer
-------------
Component responsible for results from analysis visualization. Available visualizations:
    - histograms
    - box-plots
    - scatter plot

All plots are available also in animated/interactive forms.


:Input:

Object of SeqVisualizer class is created by providing DataFrame object obtained from
:func:`SeqAnalyzer.results() <Bio.SeqAnalysis.SeqAnalyzer.results>` method, it can be
combined with receiving sequences from SeqDatabase object:

>>> from Bio import SeqAnalysis
>>> seq_db = SeqAnalysis.database(
...     [
...         'YP_025292.1',
...         '1JOY'
...     ]
... )
>>> records = seq_db.get()
>>> seq_analyzer = SeqAnalysis.analyzer(records)
>>> results = seq_analyzer.results()
>>> seq_visualizer = SeqAnalysis.visualizer(results)


:Output:

Specific visualization can be obtained by calling proper method:

>>> seq_visualizer.histogram()
>>> seq_visualizer.box_plot()
>>> seq_visualizer.scatter_plot()


User can provide additional parameters if want to save visualization:
>>> seq_visualizer.histogram(output='seqs_hist.html')

Run animation:
>>> seq_visualizer.box_plot(animation=True)

See plot in interactive mode:
>>> seq_visualizer.scatter_plot(interactive=True)
"""

import os
import shutil
import sys

from Bio.Seq import Seq
from Bio.SeqAnalysis.SeqDatabase import SeqDatabase, SeqDownloader
from Bio.SeqAnalysis.SeqAnalyzer import SeqAnalyzer
from Bio.SeqAnalysis.SeqVisualizer import SeqVisualizer
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def print_success(msg: str):
    print(f'\033[92m[+] {msg}\033[0m')


def test_database():
    # Declaration of variables needed for running tests.
    filename = "./test.txt"  # Additional file made only for reading function test.
    dirs = ['Downloads', 'Failed', 'Database']
    test_dir = "test"
    test_tuple = ('./Downloads/test/A0A1A2.fasta', 'https://www.uniprot.org/uniprot/A0A1A2.fasta')
    test_obs = []
    test_http = []

    def blockPrint():
        """
        Prevents some functions from printing their outputs.
        """
        sys.stdout = open(os.devnull, 'w')

    def test_read_list_of_ids():
        """
        Check if the function reads the file correctly.
        """
        entries = SeqDownloader.read_list_of_ids(filename)
        assert entries == ['A0A1A2', 'B0B1B2', 'C0C1C2', 'D0D1D2', 'E0E1E2', 'F0F1F2', 'G0G1G2', 'H0H1H2',
                           'I0I1I2', 'J0J1J2'], \
            'Function execution failure. Check if "test.txt" is in the same directory.'

    def test_mk_dirs():
        """
        Check if the function create directories correctly
        """
        for single_dir in dirs:
            SeqDownloader.mk_dirs(single_dir)
        for single_dir in dirs:
            assert os.path.exists(
                single_dir) is True, \
                "Function execution failure. Check if all test directories were deleted before running test."

    def test_mk_subdirs():
        """
        Check if the function create directories correctly
        """
        SeqDownloader.mk_subdirs(test_dir)
        assert os.path.exists(
            "./Downloads/test") is True, \
            "Function execution failure. Check if all test directories were deleted before running test."
        assert os.path.exists(
            "./Failed/test") is True, \
            "Function execution failure. Check if all test directories were deleted before running test."

    def test_paths_and_urls():
        """
        Check if the function create proper tuple of path and url.
        """
        urls = SeqDownloader.paths_and_urls(["A0A1A2"], test_dir)
        assert urls == [('./Downloads/test/A0A1A2.fasta',
                         'https://www.uniprot.org/uniprot/A0A1A2.fasta')], \
            "Function execution failure."

    def test_fetch_url():
        """
        Check if the function returns file path, as well as ,downloads and saves file correctly.
        """
        fetched_url = SeqDownloader.fetch_url(test_tuple, test_obs, test_http)
        assert fetched_url == './Downloads/test/A0A1A2.fasta', \
            "Function execution failure."
        assert os.path.exists("./Downloads/test/A0A1A2.fasta") is True, \
            "Function execution failure."

    def test_download_ncbi():
        """
        Check if the function downloads and saves NCBI file correctly using Entrez.
        """
        SeqDownloader.download_ncbi("YP_025292.1", test_dir, test_obs, test_http)
        assert os.path.exists("./Downloads/test/YP_025292.1.fasta") is True, \
            "Function execution failure."

    def test_SeqDownloader_class():
        """
        Check if the SeqDownloader class downloads and saves files correctly.
        """
        SeqDownloader(["A0A1A2", "YP_025292.1"], test_dir)
        assert os.path.exists(
            "./Downloads/test/A0A1A2.fasta") is True, \
            "Class execution failure. Something wrong with one of the files."
        assert os.path.exists(
            "./Downloads/test/YP_025292.1.fasta") is True, \
            "Class execution failure. Something wrong with one of the files."

    def test_get_and_SeqDatabase_class():
        """
        Check if the class SeqDatabase creates dict. Also checks if get() function works properly.
        """
        test_database = SeqDatabase(["A0A1A2", "YP_025292.1"], test_dir)
        assert type(test_database.database) == dict, \
            "Class execution failure."

        record = list(test_database.get().keys())[0:1]
        same_record = test_database.get(["A0A1A2", "YP_025292.1"])
        assert record != same_record, \
            "Function execution failure. Records are not the same."

    def rm_dirs():
        """
        Quick cleanup after certain tests and at the end of testing.
        """
        shutil.rmtree("./Database")
        shutil.rmtree("./Downloads")
        shutil.rmtree("./Failed")

    mod = 'SeqDatabase'
    test_read_list_of_ids()
    print_success(f'[{mod}] Passed: reading list test')
    test_mk_dirs()
    test_mk_subdirs()
    print_success(f'[{mod}] Passed: making dirs/subdirs test')
    test_paths_and_urls()
    test_fetch_url()
    test_download_ncbi()
    print_success(f'[{mod}] Passed: fetch and download sequences test')
    rm_dirs()
    test_SeqDownloader_class()
    print_success(f'[{mod}] Passed: SeqDownloader class test')
    rm_dirs()
    test_get_and_SeqDatabase_class()
    print_success(f'[{mod}] Passed: SeqDatabase class test')
    rm_dirs()


def test_analyzer():
    # Declaring round value - values are rounded to round_value decimal places
    round_value = 3

    # Creating test dictionary for testing methods from SeqAnalyzer class
    test_data_dict_of_str = {
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
    # Creating test dictionary for testing method from SeqAnalyzer class. Contains only proper protein sequences
    test_data_dict_of_str_only_proper_sequence = {
        'id4': test_data_dict_of_str['id4'],
        'id5': test_data_dict_of_str['id5'],
        'id6': test_data_dict_of_str['id6']
    }

    # Creating test dictionary for testing SeqAnalyzer class
    seq1 = {'id1': SeqRecord(Seq(test_data_dict_of_str['id1'], IUPAC.protein))}
    seq2 = {'id2': SeqRecord(Seq(test_data_dict_of_str['id2'], IUPAC.protein))}
    seq3 = {'id3': SeqRecord(Seq(test_data_dict_of_str['id3'], IUPAC.protein))}
    seq4 = {'id4': SeqRecord(Seq(test_data_dict_of_str['id4'], IUPAC.protein))}
    seq5 = {'id5': SeqRecord(Seq(test_data_dict_of_str['id5'], IUPAC.protein))}
    seq6 = {'id6': SeqRecord(Seq(test_data_dict_of_str['id6'], IUPAC.protein))}
    seq7 = {'id7': SeqRecord(Seq(test_data_dict_of_str['id7'], IUPAC.protein))}

    # making dict from created SeqRecord objects
    test_data_dict_of_seqrecords = {**seq1, **seq2, **seq3, **seq4, **seq5, **seq6, **seq7}

    def test_seqanalyzer_class():
        """
        Function that tests the overall functioning of SeqAnalyzer class and SeqAnalyzer.results() function

        :return: returns nothing if class passes tests, else function returns text message
        """
        results = SeqAnalyzer(test_data_dict_of_seqrecords)
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

    def test_seqanalyzer_median_standard_deviation_average_length():
        """
        Function that tests if the median_standard_deviation_average_length() from SeqAnalyzer class works properly

        :return: returns nothing if class passes tests, else function returns text message
        """
        results = SeqAnalyzer.median_standard_deviation_average_length(test_data_dict_of_str, round_value)
        assert results == (555.45, 268.0, 36.0), \
            "function median_standard_deviation_average_length() didn't pass the test"

    def test_seqanalyzer_main_calculations():
        """
        Function that tests if the main_calculations() from SeqAnalyzer class works properly

        :return: returns nothing if class passes tests, else function returns text message
        """

        results = SeqAnalyzer.main_calculations(test_data_dict_of_str_only_proper_sequence, round_value)
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

    def test_seqanalyzer_terminal_aa_counter():
        """
        Function that tests if the terminal_aa_counter() from SeqAnalyzer class works properly

        :return: returns nothing if class passes tests, else function returns text message
        """
        results = SeqAnalyzer.terminal_aa_counter(test_data_dict_of_str_only_proper_sequence, round_value)
        assert results == ({'A': 0.333, 'H': 0.333, 'V': 0.333}, {'M': 1.0, 'G': 0.5}), \
            "function test_SeqAnalyzer_terminal_aa_counter() didn't pass the test"

    """
    If the script is run directly it performs multiple test of SeqAnalyzer class and its methods. If the SeqAnalyzer
    class passes tests the functions return nothing. If any one the tests fails the functions return messages about
    which part of the SeqAnalyzer class failed to produce correct results
    """
    mod = 'SeqAnalyzer'
    test_seqanalyzer_class()
    print_success(f'[{mod}] Passed: SeqAnalyzer class test')
    test_seqanalyzer_median_standard_deviation_average_length()
    print_success(f'[{mod}] Passed: median/std/avg lengths test')
    test_seqanalyzer_main_calculations()
    print_success(f'[{mod}] Passed: amino acids analysis test')
    test_seqanalyzer_terminal_aa_counter()
    print_success(f'[{mod}] Passed: terminals amino acids counter test')
