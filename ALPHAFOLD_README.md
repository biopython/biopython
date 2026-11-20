# AlphaFoldList Module

## Overview

`Bio.PDB.AlphaFoldList` provides a PDBList-like interface for accessing and downloading protein structures from the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/).

This module addresses Issue #3990 and enables users to:
- Download individual AlphaFold predictions by UniProt accession
- Bulk download complete organism proteomes
- Access confidence scores (pLDDT) and predicted aligned error (PAE) data
- Query the AlphaFold API for structure metadata
- Integrate seamlessly with existing Bio.PDB parsers

## Quick Start

```python
from Bio.PDB.AlphaFoldList import AlphaFoldList

# Initialize
afl = AlphaFoldList(local_dir="alphafold_structures")

# Download a single structure
structure_file = afl.retrieve_alphafold_structure("P00520", file_format="pdb")

# Download with confidence data
structure_file = afl.retrieve_alphafold_structure(
    "P00520", include_pae=True, include_confidence=True
)
```

## Features

### 1. Single Structure Download

Download individual protein structures by UniProt accession:

```python
afl = AlphaFoldList()

# PDB format
pdb_file = afl.retrieve_alphafold_structure("P00520", file_format="pdb")

# mmCIF format (default)
cif_file = afl.retrieve_alphafold_structure("P00520", file_format="mmCif")

# Specific version
old_version = afl.retrieve_alphafold_structure("P00520", version=4)
```

### 2. Bulk Downloads

Download multiple structures in parallel:

```python
afl = AlphaFoldList()

# List of UniProt IDs
proteins = ["P00520", "P00533", "P04629", "P06239"]

# Download with threading
afl.download_structures_by_accessions(proteins, file_format="mmCif", max_num_threads=4)
```

### 3. Proteome Downloads

Download complete organism proteomes from the FTP server:

```python
afl = AlphaFoldList()

# List available proteomes
proteomes = afl.get_available_proteomes()
print(f"Found {len(proteomes)} proteomes")

# Download E. coli proteome (~456 MB)
afl.download_proteome("ECOLI", extract=True)

# Download human proteome (~4.8 GB - large!)
# afl.download_proteome("HUMAN", extract=True)
```

### 4. Quality Metrics

Download confidence scores and predicted aligned error:

```python
afl = AlphaFoldList()

structure = afl.retrieve_alphafold_structure(
    "P00520",
    include_pae=True,  # Predicted Aligned Error matrix
    include_confidence=True,  # Per-residue pLDDT scores
)
```

These JSON files can be used to:
- Assess prediction reliability
- Identify confident vs. uncertain regions
- Evaluate domain-domain interactions

### 5. API Queries

Query structure metadata before downloading:

```python
afl = AlphaFoldList()

metadata = afl.get_structure_metadata("P00520")
if metadata:
    entry = metadata[0]
    print(f"Protein: {entry['uniprotDescription']}")
    print(f"Organism: {entry['organismScientificName']}")
    print(f"Latest version: {entry['latestVersion']}")
    print(f"Confidence: {entry['globalMetricValue']}")
```

## Integration with Bio.PDB

AlphaFoldList works seamlessly with BioPython's structure parsers:

```python
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.AlphaFoldList import AlphaFoldList

# Download structure
afl = AlphaFoldList()
pdb_file = afl.retrieve_alphafold_structure("P00520", file_format="pdb")

# Parse with Bio.PDB
parser = PDBParser()
structure = parser.get_structure("ABL1", pdb_file)

# Analyze
for model in structure:
    for chain in model:
        print(f"Chain {chain.id}: {len(chain)} residues")
```

## Configuration Options

### Directory Organization

Two modes available:

1. **Hierarchical tree** (PDB-style, default):
```python
afl = AlphaFoldList(local_dir="structures")
afl.flat_tree = False  # default
# Creates: structures/p0/AF-P00520-F1-model_v6.cif
```

2. **Flat directory**:
```python
afl = AlphaFoldList(local_dir="structures")
afl.flat_tree = True
# Creates: structures/AF-P00520-F1-model_v6.cif
```

3. **Custom directory** per structure:
```python
afl.retrieve_alphafold_structure("P00520", pdir="custom/path/")
```

### Server Configuration

```python
afl = AlphaFoldList(
    api_server="https://alphafold.ebi.ac.uk",
    ftp_server="https://ftp.ebi.ac.uk/pub/databases/alphafold",
    local_dir="/path/to/structures",
    verbose=True,  # Enable progress messages
)
```

## Command-Line Usage

AlphaFoldList can be used directly from the command line:

```bash
# Download single structure
python -m Bio.PDB.AlphaFoldList P00520

# Download in PDB format with extras
python -m Bio.PDB.AlphaFoldList P00520 -pdb -pae -conf

# List available proteomes
python -m Bio.PDB.AlphaFoldList proteomes

# Download organism proteome
python -m Bio.PDB.AlphaFoldList proteome HUMAN /path/to/save

# Specific version
python -m Bio.PDB.AlphaFoldList P00520 -v 4

# Custom directory
python -m Bio.PDB.AlphaFoldList P00520 -d /custom/path
```

## API Reference

### Class: AlphaFoldList

#### `__init__(api_server, ftp_server, local_dir, verbose)`

Initialize the AlphaFoldList instance.

**Parameters:**
- `api_server` (str): AlphaFold API URL (default: https://alphafold.ebi.ac.uk)
- `ftp_server` (str): FTP server URL (default: https://ftp.ebi.ac.uk/pub/databases/alphafold)
- `local_dir` (str): Local directory for downloads (default: current directory)
- `verbose` (bool): Enable verbose output (default: True)

#### `get_structure_metadata(uniprot_id)`

Retrieve metadata for a structure from the AlphaFold API.

**Parameters:**
- `uniprot_id` (str): UniProt accession (e.g., "P00520")

**Returns:** List of dict containing structure metadata

#### `retrieve_alphafold_structure(uniprot_id, pdir, file_format, version, overwrite, include_pae, include_confidence)`

Download a single AlphaFold structure.

**Parameters:**
- `uniprot_id` (str): UniProt accession
- `pdir` (str, optional): Custom directory
- `file_format` (str): "mmCif" (default) or "pdb"
- `version` (int, optional): Specific version number
- `overwrite` (bool): Overwrite existing files (default: False)
- `include_pae` (bool): Download PAE JSON (default: False)
- `include_confidence` (bool): Download confidence JSON (default: False)

**Returns:** Path to downloaded file (str) or None if failed

#### `download_structures_by_accessions(uniprot_ids, pdir, file_format, overwrite, max_num_threads)`

Bulk download multiple structures in parallel.

**Parameters:**
- `uniprot_ids` (list): List of UniProt accessions
- `pdir` (str, optional): Custom directory
- `file_format` (str): File format
- `overwrite` (bool): Overwrite existing
- `max_num_threads` (int, optional): Max parallel downloads

#### `get_available_proteomes()`

List all available organism proteomes from FTP.

**Returns:** List of dict with proteome information

#### `download_proteome(organism_code, pdir, extract, overwrite)`

Download complete organism proteome.

**Parameters:**
- `organism_code` (str): Organism code (e.g., "HUMAN", "MOUSE")
- `pdir` (str, optional): Download directory
- `extract` (bool): Extract tar archive (default: True)
- `overwrite` (bool): Overwrite existing (default: False)

**Returns:** Path to downloaded/extracted files

## File Formats

### Structure Files

- **PDB**: Traditional PDB format (human-readable, limited to 99,999 atoms)
- **mmCIF**: Modern macromolecular CIF format (recommended for large structures)

### Quality Metrics Files

- **PAE JSON** (`*-predicted_aligned_error_v*.json`): 
  - Matrix of predicted aligned errors between residue pairs
  - Useful for assessing domain interactions and overall fold confidence
  
- **Confidence JSON** (`*-confidence_v*.json`):
  - Per-residue pLDDT (predicted local distance difference test) scores
  - Scale: 0-100, higher = more confident
  - >90: very high confidence
  - 70-90: confident
  - 50-70: low confidence
  - <50: very low confidence

## Known Organisms/Proteomes

Common organisms available for bulk download:

| Organism | Code | Size | Description |
|----------|------|------|-------------|
| Human | HUMAN | 4.8 GB | Homo sapiens |
| Mouse | MOUSE | 3.5 GB | Mus musculus |
| E. coli | ECOLI | 456 MB | Escherichia coli |
| Yeast | YEAST | 1.0 GB | Saccharomyces cerevisiae |
| Fruit fly | DROME | 2.2 GB | Drosophila melanogaster |
| C. elegans | CAEEL | 2.6 GB | Caenorhabditis elegans |
| Zebrafish | DANRE | 4.6 GB | Danio rerio |
| Arabidopsis | ARATH | 3.6 GB | Arabidopsis thaliana |

...and 40+ more organisms. Use `get_available_proteomes()` for the complete list.

## Comparison with PDBList

`AlphaFoldList` follows the same design patterns as `PDBList`:

| Feature | PDBList | AlphaFoldList |
|---------|---------|---------------|
| Single structure download | ✓ | ✓ |
| Bulk download | ✓ | ✓ |
| Multiple formats | ✓ | ✓ |
| Directory tree | ✓ | ✓ |
| Threading support | ✓ | ✓ |
| API queries | Limited | ✓ |
| Confidence data | N/A | ✓ |

## References

1. **AlphaFold Database:**
   - Website: https://alphafold.ebi.ac.uk/
   - API Docs: https://alphafold.ebi.ac.uk/api-docs

2. **Publications:**
   - Jumper, J. et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583–589.
   - Varadi, M. et al. (2022). AlphaFold Protein Structure Database. *Nucleic Acids Research* 50, D439–D444.

3. **License:**
   - AlphaFold Database: CC-BY-4.0
   - BioPython: Biopython License Agreement / BSD 3-Clause

## Troubleshooting

### Connection errors
```python
# Check internet connection
# Verify API server is accessible:
import urllib.request

urllib.request.urlopen("https://alphafold.ebi.ac.uk/").read()
```

### Large downloads timeout
```python
# Increase timeout or download in smaller batches
afl.download_structures_by_accessions(
    large_list[:100], max_num_threads=2  # First 100  # Reduce load
)
```

### Disk space
Proteome downloads can be very large (up to 7+ GB). Check available space:
```bash
df -h  # Unix/Linux
```

## Contributing

This module addresses Issue #3990. For bugs or feature requests, please open an issue on the BioPython GitHub repository.

## See Also

- `Bio.PDB.PDBList` - Similar interface for experimental PDB structures
- `Bio.PDB.alphafold_db` - Minimal AlphaFold database interface
- `Bio.PDB.PDBParser` - Parse PDB format structures
- `Bio.PDB.MMCIFParser` - Parse mmCIF format structures
