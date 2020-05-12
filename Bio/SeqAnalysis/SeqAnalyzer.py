"""
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
"""

from Bio import SeqAnalysis

seq_db = SeqAnalysis.database(['YP_025292.1', '1JOY'])
records = seq_db.get()

# making a dict out of records, where id is the key and its sequence its value
dict = {}
for _ in records:
    dict[_[0]] = (_[1].seq)[0]



class analyzer():
    def __init__(self):






seq_analyzer = SeqAnalysis.analyzer(records)

