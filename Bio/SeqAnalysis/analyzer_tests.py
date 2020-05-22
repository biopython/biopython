from Bio.SeqAnalysis.SeqAnalyzer import SeqAnalyzer
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

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
    assert results.results() == ({'median_length': 158.0, 'average_length': 605.667, 'standard_deviation_length': 721.085, 'id': ['id4', 'id5', 'id6'], 'sequence_length': [158, 36, 1623], 'molecular_weight': [17245.603, 3884.478, 177765.884], 'theoretical_pi': [5.197, 5.518, 5.556], 'instability_index': [28.522, 9.625, 33.646], 'flexibility': [0.994, 1.004, 0.998], 'gravy': [0.042, 0.064, 0.074], 'aromaticity': [0.082, 0.056, 0.076], 'c_term_freq': {'M': 1.0, 'G': 0.5}, 'n_term_freq': {'A': 0.333, 'H': 0.333, 'V': 0.333}, 'secondary_structure_fractions': [{'Helix': 31.646, 'Turn': 23.418, 'Sheet': 23.418}, {'Helix': 36.111, 'Turn': 25.0, 'Sheet': 16.667}, {'Helix': 33.826, 'Turn': 22.058, 'Sheet': 26.802}], 'extinction_coefficient': [{'Cysteines_reduced': 18450, 'Cysteines_residues': 18950}, {'Cysteines_reduced': 1490, 'Cysteines_residues': 1490}, {'Cysteines_reduced': 133620, 'Cysteines_residues': 134995}], 'amino_acid_comp': [{'A': 10, 'C': 8, 'D': 9, 'E': 11, 'F': 6, 'G': 17, 'H': 9, 'I': 9, 'K': 6, 'L': 12, 'M': 4, 'N': 5, 'P': 11, 'Q': 3, 'R': 4, 'S': 4, 'T': 7, 'V': 16, 'W': 2, 'Y': 5}, {'A': 1, 'C': 1, 'D': 3, 'E': 2, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4, 'L': 2, 'M': 1, 'N': 2, 'P': 3, 'Q': 0, 'R': 0, 'S': 1, 'T': 1, 'V': 9, 'W': 0, 'Y': 1}, {'A': 129, 'C': 23, 'D': 91, 'E': 98, 'F': 72, 'G': 122, 'H': 25, 'I': 131, 'K': 105, 'L': 171, 'M': 37, 'N': 72, 'P': 51, 'Q': 55, 'R': 56, 'S': 113, 'T': 97, 'V': 123, 'W': 14, 'Y': 38}], 'amino_acid_percent': [{'A': 6.329, 'C': 5.063, 'D': 5.696, 'E': 6.962, 'F': 3.797, 'G': 10.759, 'H': 5.696, 'I': 5.696, 'K': 3.797, 'L': 7.595, 'M': 2.532, 'N': 3.165, 'P': 6.962, 'Q': 1.899, 'R': 2.532, 'S': 2.532, 'T': 4.43, 'V': 10.127, 'W': 1.266, 'Y': 3.165}, {'A': 2.778, 'C': 2.778, 'D': 8.333, 'E': 5.556, 'F': 2.778, 'G': 8.333, 'H': 2.778, 'I': 0.0, 'K': 11.111, 'L': 5.556, 'M': 2.778, 'N': 5.556, 'P': 8.333, 'Q': 0.0, 'R': 0.0, 'S': 2.778, 'T': 2.778, 'V': 25.0, 'W': 0.0, 'Y': 2.778}, {'A': 7.948, 'C': 1.417, 'D': 5.607, 'E': 6.038, 'F': 4.436, 'G': 7.517, 'H': 1.54, 'I': 8.071, 'K': 6.47, 'L': 10.536, 'M': 2.28, 'N': 4.436, 'P': 3.142, 'Q': 3.389, 'R': 3.45, 'S': 6.962, 'T': 5.977, 'V': 7.579, 'W': 0.863, 'Y': 2.341}]}, {'id4': {'id': ['id4'], 'sequence_length': [158], 'molecular_weight': [17245.603], 'theoretical_pi': [5.197], 'instability_index': [28.522], 'flexibility': [0.994], 'gravy': [0.042], 'aromaticity': [0.082], 'secondary_structure_fraction': [{'Helix': 31.646, 'Turn': 23.418, 'Sheet': 23.418}], 'extinction_coefficient': [{'Cysteines_reduced': 18450, 'Cysteines_residues': 18950}], 'amino_acid_comp': [{'A': 10, 'C': 8, 'D': 9, 'E': 11, 'F': 6, 'G': 17, 'H': 9, 'I': 9, 'K': 6, 'L': 12, 'M': 4, 'N': 5, 'P': 11, 'Q': 3, 'R': 4, 'S': 4, 'T': 7, 'V': 16, 'W': 2, 'Y': 5}], 'amino_acid_percent': [{'A': 6.329, 'C': 5.063, 'D': 5.696, 'E': 6.962, 'F': 3.797, 'G': 10.759, 'H': 5.696, 'I': 5.696, 'K': 3.797, 'L': 7.595, 'M': 2.532, 'N': 3.165, 'P': 6.962, 'Q': 1.899, 'R': 2.532, 'S': 2.532, 'T': 4.43, 'V': 10.127, 'W': 1.266, 'Y': 3.165}], 'c_term_freq': {'H': 1}, 'n_term_freq': {'M': 1}}, 'id5': {'id': ['id5'], 'sequence_length': [36], 'molecular_weight': [3884.478], 'theoretical_pi': [5.518], 'instability_index': [9.625], 'flexibility': [1.004], 'gravy': [0.064], 'aromaticity': [0.056], 'secondary_structure_fraction': [{'Helix': 36.111, 'Turn': 25.0, 'Sheet': 16.667}], 'extinction_coefficient': [{'Cysteines_reduced': 1490, 'Cysteines_residues': 1490}], 'amino_acid_comp': [{'A': 1, 'C': 1, 'D': 3, 'E': 2, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4, 'L': 2, 'M': 1, 'N': 2, 'P': 3, 'Q': 0, 'R': 0, 'S': 1, 'T': 1, 'V': 9, 'W': 0, 'Y': 1}], 'amino_acid_percent': [{'A': 2.778, 'C': 2.778, 'D': 8.333, 'E': 5.556, 'F': 2.778, 'G': 8.333, 'H': 2.778, 'I': 0.0, 'K': 11.111, 'L': 5.556, 'M': 2.778, 'N': 5.556, 'P': 8.333, 'Q': 0.0, 'R': 0.0, 'S': 2.778, 'T': 2.778, 'V': 25.0, 'W': 0.0, 'Y': 2.778}], 'c_term_freq': {'A': 1}, 'n_term_freq': {'G': 1}}, 'id6': {'id': ['id6'], 'sequence_length': [1623], 'molecular_weight': [177765.884], 'theoretical_pi': [5.556], 'instability_index': [33.646], 'flexibility': [0.998], 'gravy': [0.074], 'aromaticity': [0.076], 'secondary_structure_fraction': [{'Helix': 33.826, 'Turn': 22.058, 'Sheet': 26.802}], 'extinction_coefficient': [{'Cysteines_reduced': 133620, 'Cysteines_residues': 134995}], 'amino_acid_comp': [{'A': 129, 'C': 23, 'D': 91, 'E': 98, 'F': 72, 'G': 122, 'H': 25, 'I': 131, 'K': 105, 'L': 171, 'M': 37, 'N': 72, 'P': 51, 'Q': 55, 'R': 56, 'S': 113, 'T': 97, 'V': 123, 'W': 14, 'Y': 38}], 'amino_acid_percent': [{'A': 7.948, 'C': 1.417, 'D': 5.607, 'E': 6.038, 'F': 4.436, 'G': 7.517, 'H': 1.54, 'I': 8.071, 'K': 6.47, 'L': 10.536, 'M': 2.28, 'N': 4.436, 'P': 3.142, 'Q': 3.389, 'R': 3.45, 'S': 6.962, 'T': 5.977, 'V': 7.579, 'W': 0.863, 'Y': 2.341}], 'c_term_freq': {'V': 1}, 'n_term_freq': {'M': 1}}}),\
    'class didn\nt pass the test'


def test_seqanalyzer_median_standard_deviation_average_length():
    """
    Function that tests if the median_standard_deviation_average_length() from SeqAnalyzer class works properly

    :return: returns nothing if class passes tests, else function returns text message
    """
    results = SeqAnalyzer.median_standard_deviation_average_length(test_data_dict_of_str)
    assert results == (555.45, 268.0, 36.0), \
        'function median_standard_deviation_average_length() didn\'t pass the test'


def test_seqanalyzer_main_calculations():
    """
    Function that tests if the main_calculations() from SeqAnalyzer class works properly

    :return: returns nothing if class passes tests, else function returns text message
    """

    results = SeqAnalyzer.main_calculations(test_data_dict_of_str_only_proper_sequence)
    assert results == ([17245.603, 3884.478, 177765.884], [5.197, 5.518, 5.556], [{'A': 10, 'C': 8, 'D': 9, 'E': 11, 'F': 6, 'G': 17, 'H': 9, 'I': 9, 'K': 6, 'L': 12, 'M': 4, 'N': 5, 'P': 11, 'Q': 3, 'R': 4, 'S': 4, 'T': 7, 'V': 16, 'W': 2, 'Y': 5}, {'A': 1, 'C': 1, 'D': 3, 'E': 2, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4, 'L': 2, 'M': 1, 'N': 2, 'P': 3, 'Q': 0, 'R': 0, 'S': 1, 'T': 1, 'V': 9, 'W': 0, 'Y': 1}, {'A': 129, 'C': 23, 'D': 91, 'E': 98, 'F': 72, 'G': 122, 'H': 25, 'I': 131, 'K': 105, 'L': 171, 'M': 37, 'N': 72, 'P': 51, 'Q': 55, 'R': 56, 'S': 113, 'T': 97, 'V': 123, 'W': 14, 'Y': 38}], [0.082, 0.056, 0.076], [{'A': 6.329, 'C': 5.063, 'D': 5.696, 'E': 6.962, 'F': 3.797, 'G': 10.759, 'H': 5.696, 'I': 5.696, 'K': 3.797, 'L': 7.595, 'M': 2.532, 'N': 3.165, 'P': 6.962, 'Q': 1.899, 'R': 2.532, 'S': 2.532, 'T': 4.43, 'V': 10.127, 'W': 1.266, 'Y': 3.165}, {'A': 2.778, 'C': 2.778, 'D': 8.333, 'E': 5.556, 'F': 2.778, 'G': 8.333, 'H': 2.778, 'I': 0.0, 'K': 11.111, 'L': 5.556, 'M': 2.778, 'N': 5.556, 'P': 8.333, 'Q': 0.0, 'R': 0.0, 'S': 2.778, 'T': 2.778, 'V': 25.0, 'W': 0.0, 'Y': 2.778}, {'A': 7.948, 'C': 1.417, 'D': 5.607, 'E': 6.038, 'F': 4.436, 'G': 7.517, 'H': 1.54, 'I': 8.071, 'K': 6.47, 'L': 10.536, 'M': 2.28, 'N': 4.436, 'P': 3.142, 'Q': 3.389, 'R': 3.45, 'S': 6.962, 'T': 5.977, 'V': 7.579, 'W': 0.863, 'Y': 2.341}], [28.522, 9.625, 33.646], [0.994, 1.004, 0.998], [{'Helix': 31.646, 'Turn': 23.418, 'Sheet': 23.418}, {'Helix': 36.111, 'Turn': 25.0, 'Sheet': 16.667}, {'Helix': 33.826, 'Turn': 22.058, 'Sheet': 26.802}], [0.042, 0.064, 0.074], [{'Cysteines_reduced': 18450, 'Cysteines_residues': 18950}, {'Cysteines_reduced': 1490, 'Cysteines_residues': 1490}, {'Cysteines_reduced': 133620, 'Cysteines_residues': 134995}]),\
        'function main_calculations() didn\'t pass the test'


def test_seqanalyzer_terminal_aa_counter():
    """
    Function that tests if the terminal_aa_counter() from SeqAnalyzer class works properly

    :return: returns nothing if class passes tests, else function returns text message
    """
    results = SeqAnalyzer.terminal_aa_counter(test_data_dict_of_str_only_proper_sequence)
    assert results == ({'A': 0.333, 'H': 0.333, 'V': 0.333}, {'M': 1.0, 'G': 0.5}),\
        'function test_SeqAnalyzer_terminal_aa_counter() didn\'t pass the test'


if __name__ == '__main__':
    """
    If the script is run directly it performs multiple test of SeqAnalyzer class and its methods. If the SeqAnalyzer
    class passes tests the functions return nothing. If any one the tests fails the functions return messages about
    which part of the SeqAnalyzer class failed to produce correct results
    """
    test_seqanalyzer_class()
    test_seqanalyzer_median_standard_deviation_average_length()
    test_seqanalyzer_main_calculations()
    test_seqanalyzer_terminal_aa_counter()
