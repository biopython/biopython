# SeqAnalysis

Proteins sequences analysis and visualization

**Bio.SeqAnalysis** provides python components for storing, analyzing and visualizing proteins sequences.
Components can be used separately or through API as fully working, automated pipeline.

## Components
There are 3 main components which basically are 3 classes:
- **SeqDatabase**: downloading, storing and providing protein sequences
- **SeqAnalyzer**: complex analysis of protein sequences
- **SeqVisualizer**: protein sequences analysis report visualization

## Usage

All you need is file with sequences Uniprot IDs which you'd like to analyze. You can use **SeqAnalysis** as a single
module or (when it'll come to official biopython) as biopython module.

### Build

To use as single module, you need to extract folder with this module from this biopython branch and install requirements:

```bash
pip3 install virtualenv

# extracting module
git clone --single-branch --branch sequences_analysis https://github.com/moozeq/biopython.git
mv biopython/Bio/SeqAnalysis .

# install requirements
python3 -m venv venv
source venv/bin/activate
pip3 install biopython matplotlib plotly pandas numpy
```

### Run

You need to specify name for `FILE_WITH_SEQUENCES` first (example file: `Tests/SeqAnalysis/sequences.txt`):

```bash
source venv/bin/activate
python3 -m SeqAnalysis ${FILE_WITH_SEQUENCES} -o
```

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
cp venv/lib/*/site-packages/biopython*/Bio/Align/* Bio/Align
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