# SeqAnalysis

## Description

Proteins sequences analysis and visualization

**Bio.SeqAnalysis** provides python components for storing, analyzing and visualizing proteins sequences.
Components can be used separately or through API as fully working, automated pipeline.

SeqAnalysis module was made for proteomic analysis performed on protein sequences downloaded from NCBI. Due to the
ability of efficient downloading great amount of protein sequences SeqAnalysis can be used to analyze/compare proteins,
from single proteins, through protein families, to whole proteomes of organisms. Thanks to the complex visualization
capacities of SeqAnalysis module of calculated results, users can swiftly asses said results, moreover graphical
representation of results makes them easier to read for users less familiar with proteomic analysis. SeqAnalysis module
can be also used as a efficient tool for downloading huge chunks of protein sequences from NCBI and saving them in a
txt file, so that those sequences can be used for further processing. Therefore SeqAnalysis module can be used as a fast
sequences downloading tool with a built in proteomic pre-processing analysis.

The SeqAnalysis:
- Calculates **molecular weight** of protein sequences
- Calculates **isoelectric point** for protein sequences
- Counts **standard aa** in protein sequences
- Calculates **aromaticity** for protein sequences according to Lobry 1994 method
- Calculates **amino acid percentage** content based aa content in protein sequences
- Calculates **instability indexes** for protein sequences according to Guruprasad et al 1990
- Calculates **flexibility** of protein sequences according to Vihinen, 1994
- Calculates **percentage of helix, turn and sheet** in protein sequences
- Calculates **GRAVY** of protein sequences according to Kyte and Doolittle
- Calculates **molar extinction coefficients** for protein sequences
- Calculates **median length** of sequences
- Calculates **average length** of sequences
- Calculates **standard deviation** of sequence length
- Calculates **'N' and 'C' terminus** amino acid occurrences frequencies

## Components
There are 3 main components which basically are 3 classes:
- **SeqDatabase**: downloading, storing and providing protein sequences
- **SeqAnalyzer**: complex analysis of protein sequences
- **SeqVisualizer**: protein sequences analysis report visualization

## Usage

All you need is file with sequences Uniprot IDs which you'd like to analyze. You can use **SeqAnalysis** as a single
module or (when it'll come to official biopython) as biopython module.

### Build

To use as single module, you need to extract folder with this module from this biopython branch and install requirements
(just copy and paste to terminal script below):

```bash
pip3 install virtualenv

# extracting module
git clone --single-branch --branch sequences_analysis https://github.com/moozeq/biopython.git
mv biopython/Bio/SeqAnalysis .

# extracting example file
mv biopython/Tests/SeqAnalysis/sequences.txt .

# install requirements
python3 -m venv venv
source venv/bin/activate
pip3 install biopython matplotlib plotly pandas numpy
```

### Run

To run whole analysis you only need to specify filename of file with sequences IDs (below example file *sequences.txt*):

```bash
python3 -m SeqAnalysis sequences.txt
```

IDs in file should be separated. Available delimiters: `';' ':' ',' '\n' '\t'`.

## Development

For development purposes biopython needs to be built and installed. Some files are required to be compiled and copied to
proper folders.

### Build

You need to specify name for `BRANCH_FOR_DEVELOPMENT` first:

```bash
pip3 install virtualenv

git clone https://github.com/moozeq/biopython.git
cd biopython
git checkout -b ${BRANCH_FOR_DEVELOPMENT} origin/sequences_analysis
python3 -m venv venv
source venv/bin/activate
python3 setup.py build
python3 setup.py install
cp -r venv/lib/*/site-packages/biopython*/Bio/Align/* Bio/Align
```

### Tests

Tests are in biopython `Tests` directory and are performed during `python3 setup.py test`.
When wants to run tests only for `SeqAnalysis` module, you need to run them from `unittest`:

```bash
python3 -m unittest Tests/test_SeqAnalysis.py
```

When running from Pycharm, set working directory to biopython, go to
`test_SeqAnalysis.py` file and run tests.

### Run

Running as single module with input from example file, analysis output should appear in `seq_analysis` directory:

```bash
source venv/bin/activate
python3 -m Bio.SeqAnalysis Tests/SeqAnalysis/sequences.txt -o
mv seq_analysis ..
```