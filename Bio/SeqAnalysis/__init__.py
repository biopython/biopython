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

from pandas import DataFrame
from typing import List
from Bio.SeqRecord import SeqRecord
from Bio.SeqAnalysis.SeqDatabase import SeqDatabase
from Bio.SeqAnalysis.SeqAnalyzer import SeqAnalyzer
from Bio.SeqAnalysis.SeqVisualizer import SeqVisualizer


def database(sequences_ids: List[str]) -> SeqDatabase:
    return SeqDatabase(sequences_ids)

#kacdom
def analyzer(sequences_records: List[SeqRecord]) -> SeqAnalyzer:
    return SeqAnalyzer(sequences_records)


def visualizer(sequences_analysis_report: DataFrame) -> SeqVisualizer:
    return SeqVisualizer(sequences_analysis_report)
