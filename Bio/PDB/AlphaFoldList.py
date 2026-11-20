#!/usr/bin/env python
# Copyright 2024, BioPython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# AlphaFoldList.py
#
# A tool for accessing and downloading AlphaFold protein structure predictions.
#
# Developed to integrate with BioPython's PDB module
# Based on the PDBList.py architecture

"""Access the AlphaFold Protein Structure Database over the internet."""

import contextlib
import functools
import gzip
import json
import os
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor
from typing import Optional
from urllib.request import Request
from urllib.request import urlcleanup
from urllib.request import urlopen
from urllib.request import urlretrieve


class AlphaFoldList:
    """Quick access to AlphaFold predicted protein structures.

    This class provides access to the AlphaFold Protein Structure Database,
    allowing users to download predicted protein structures by UniProt accession,
    organism/proteome, or in bulk.

    The AlphaFold database contains AI-predicted protein structures for millions
    of proteins across various organisms. Each structure includes:
    - 3D coordinates (PDB or mmCIF format)
    - Per-residue confidence scores (pLDDT)
    - Predicted aligned error (PAE) matrices
    - Multiple sequence alignments (MSAs)

    All structures are available under CC-BY-4.0 license.

    References:
        Jumper et al. (2021) Nature 596, 583–589
        Varadi et al. (2022) Nucleic Acids Research 50, D439–D444
        https://alphafold.ebi.ac.uk/
    """

    ALPHAFOLD_REF = """
    Highly accurate protein structure prediction with AlphaFold.
    Jumper, J et al. (2021)
    Nature 596, 583–589
    https://doi.org/10.1038/s41586-021-03819-2

    AlphaFold Protein Structure Database: massively expanding the structural coverage
    of protein-sequence space with high-accuracy models.
    Varadi, M et al. (2022)
    Nucleic Acids Research 50, D439–D444
    https://doi.org/10.1093/nar/gkab1061
    """

    def __init__(
        self,
        api_server="https://alphafold.ebi.ac.uk",
        ftp_server="https://ftp.ebi.ac.uk/pub/databases/alphafold",
        local_dir=None,
        verbose=True,
    ):
        """Initialize the AlphaFoldList with server URLs and local directory.

        :param api_server: AlphaFold API server URL
        :type api_server: string

        :param ftp_server: AlphaFold FTP server URL for bulk downloads
        :type ftp_server: string

        :param local_dir: Local directory to store downloaded structures
            (defaults to current directory)
        :type local_dir: string

        :param verbose: Enable verbose output
        :type verbose: bool
        """
        self.api_server = api_server
        self.ftp_server = ftp_server

        if local_dir:
            self.local_dir = local_dir
        else:
            self.local_dir = os.getcwd()

        self._verbose = verbose
        self.flat_tree = False  # Use directory tree structure by default

        # Cache for organism proteome listings
        self._proteome_cache = None

    def _log(self, message):
        """Print message if verbose mode is enabled (PRIVATE)."""
        if self._verbose:
            print(message)

    def get_structure_metadata(self, uniprot_id):
        """Retrieve metadata for an AlphaFold structure from the API.

        :param uniprot_id: UniProt accession (e.g., 'P00520' or 'P00520-2' for isoforms)
        :type uniprot_id: string

        :return: List of structure metadata dictionaries (may include multiple
            isoforms/fragments)
        :rtype: list of dict

        Example:
            >>> afl = AlphaFoldList()
            >>> metadata = afl.get_structure_metadata('P00520')
            >>> print(metadata[0]['entryId'])
            AF-P00520-F1
        """
        url = f"{self.api_server}/api/prediction/{uniprot_id}"
        self._log(f"Querying AlphaFold API for {uniprot_id}...")

        try:
            with contextlib.closing(urlopen(url)) as response:
                data = json.loads(response.read().decode("utf-8"))
                return data
        except Exception as e:
            if self._verbose:
                print(f"Error retrieving metadata for {uniprot_id}: {e}")
            return None

    def retrieve_alphafold_structure(
        self,
        uniprot_id,
        pdir=None,
        file_format="mmCif",
        version=None,
        overwrite=False,
        include_pae=False,
        include_confidence=False,
    ):
        """Fetch AlphaFold structure file and store it locally.

        Downloads the predicted structure for a given UniProt accession.
        By default, downloads the latest version in mmCIF format.

        :param uniprot_id: UniProt accession (e.g., 'P00520')
        :type uniprot_id: string

        :param pdir: Directory to save the file (default: local_dir with tree structure)
        :type pdir: string

        :param file_format: File format - "mmCif" (default) or "pdb"
        :type file_format: string

        :param version: Specific model version to download (default: latest)
            AlphaFold regularly updates predictions, numbered v1, v2, v3, etc.
        :type version: int

        :param overwrite: Overwrite existing files
        :type overwrite: bool

        :param include_pae: Also download predicted aligned error (PAE) JSON file
        :type include_pae: bool

        :param include_confidence: Also download per-residue confidence (pLDDT) JSON
        :type include_confidence: bool

        :return: Path to downloaded structure file, or None if failed
        :rtype: string or None

        Example:
            >>> afl = AlphaFoldList()
            >>> filename = afl.retrieve_alphafold_structure('P00520', file_format='pdb')
            >>> print(filename)
            /path/to/AF-P00520-F1-model_v4.pdb
        """
        # Get metadata to find the correct entry ID and latest version
        metadata = self.get_structure_metadata(uniprot_id)
        if not metadata:
            self._log(f"Could not retrieve metadata for {uniprot_id}")
            return None

        # Use the first entry (canonical isoform)
        entry_data = metadata[0]
        entry_id = entry_data["entryId"]  # e.g., AF-P00520-F1

        # Determine version to download
        if version is None:
            version = entry_data.get("latestVersion", 4)

        # Construct file URLs based on format
        file_format_lower = file_format.lower()
        if file_format_lower in ("mmcif", "cif"):
            archive_name = f"{entry_id}-model_v{version}.cif.gz"
            final_name = f"{entry_id}-model_v{version}.cif"
            url = f"{self.api_server}/files/{entry_id}-model_v{version}.cif"
        elif file_format_lower == "pdb":
            archive_name = f"{entry_id}-model_v{version}.pdb.gz"
            final_name = f"{entry_id}-model_v{version}.pdb"
            url = f"{self.api_server}/files/{entry_id}-model_v{version}.pdb"
        else:
            raise ValueError(
                f"Unsupported file format: {file_format}. Use 'mmCif' or 'pdb'."
            )

        # Determine save path
        if pdir is None:
            path = self.local_dir
            if not self.flat_tree:
                # Create subdirectory based on first two characters of entry
                path = os.path.join(path, entry_id[3:5].lower())
        else:
            path = pdir

        os.makedirs(path, exist_ok=True)

        final_file = os.path.join(path, final_name)

        # Skip if file exists and not overwriting
        if not overwrite and os.path.exists(final_file):
            self._log(f"Structure exists: '{final_file}'")
            return final_file

        # Download the structure file
        self._log(f"Downloading {entry_id} version {version}...")
        try:
            urlcleanup()
            # AlphaFold files may or may not be gzipped depending on source
            # Try direct download first
            urlretrieve(url, final_file)
            self._log(f"Downloaded: {final_file}")

            # Download additional files if requested
            if include_pae:
                self._download_pae(entry_id, version, path, overwrite)
            if include_confidence:
                self._download_confidence(entry_id, version, path, overwrite)

            return final_file

        except Exception as e:
            self._log(f"Download failed for {uniprot_id}: {e}")
            if os.path.exists(final_file):
                os.remove(final_file)
            return None

    def _download_pae(self, entry_id, version, path, overwrite=False):
        """Download predicted aligned error (PAE) JSON file (PRIVATE)."""
        pae_filename = f"{entry_id}-predicted_aligned_error_v{version}.json"
        pae_file = os.path.join(path, pae_filename)

        if not overwrite and os.path.exists(pae_file):
            self._log(f"PAE file exists: '{pae_file}'")
            return pae_file

        url = f"{self.api_server}/files/{pae_filename}"
        try:
            urlretrieve(url, pae_file)
            self._log(f"Downloaded PAE: {pae_file}")
            return pae_file
        except Exception as e:
            self._log(f"Failed to download PAE: {e}")
            return None

    def _download_confidence(self, entry_id, version, path, overwrite=False):
        """Download per-residue confidence (pLDDT) JSON file (PRIVATE)."""
        conf_filename = f"{entry_id}-confidence_v{version}.json"
        conf_file = os.path.join(path, conf_filename)

        if not overwrite and os.path.exists(conf_file):
            self._log(f"Confidence file exists: '{conf_file}'")
            return conf_file

        url = f"{self.api_server}/files/{conf_filename}"
        try:
            urlretrieve(url, conf_file)
            self._log(f"Downloaded confidence: {conf_file}")
            return conf_file
        except Exception as e:
            self._log(f"Failed to download confidence: {e}")
            return None

    def download_structures_by_accessions(
        self,
        uniprot_ids: list[str],
        pdir: str | None = None,
        file_format: str = "mmCif",
        overwrite: bool = False,
        max_num_threads: int | None = None,
    ):
        """Download multiple AlphaFold structures in parallel.

        :param uniprot_ids: List of UniProt accessions
        :type uniprot_ids: list of strings

        :param pdir: Directory to save files
        :type pdir: string

        :param file_format: "mmCif" (default) or "pdb"
        :type file_format: string

        :param overwrite: Overwrite existing files
        :type overwrite: bool

        :param max_num_threads: Maximum number of parallel downloads
        :type max_num_threads: int

        Example:
            >>> afl = AlphaFoldList()
            >>> ids = ['P00520', 'P12345', 'Q9Y6K9']
            >>> afl.download_structures_by_accessions(ids)
        """
        self._log(f"Downloading {len(uniprot_ids)} structures...")

        with ThreadPoolExecutor(max_num_threads) as executor:
            executor.map(
                functools.partial(
                    self.retrieve_alphafold_structure,
                    pdir=pdir,
                    file_format=file_format,
                    overwrite=overwrite,
                ),
                uniprot_ids,
            )

    def get_available_proteomes(self):
        """Retrieve list of available organism proteomes from FTP server.

        Returns information about bulk proteome downloads available
        (e.g., human, mouse, yeast, E. coli, etc.)

        :return: List of proteome dictionaries with 'filename', 'organism', 'size'
        :rtype: list of dict

        Example:
            >>> afl = AlphaFoldList()
            >>> proteomes = afl.get_available_proteomes()
            >>> human = [p for p in proteomes if 'HUMAN' in p['filename']]
            >>> print(human[0]['filename'])
            UP000005640_9606_HUMAN_v6.tar
        """
        if self._proteome_cache:
            return self._proteome_cache

        url = f"{self.ftp_server}/latest/"
        self._log("Retrieving list of available proteomes...")

        try:
            with contextlib.closing(urlopen(url)) as response:
                html = response.read().decode("utf-8")

                # Parse HTML to extract .tar files
                import re

                pattern = r'<a href="(UP\d+_\d+_\w+_v\d+\.tar)">'
                matches = re.findall(pattern, html)

                proteomes = []
                for filename in matches:
                    # Extract organism code from filename
                    parts = filename.replace(".tar", "").split("_")
                    proteomes.append(
                        {
                            "filename": filename,
                            "uniprot_id": parts[0],
                            "taxid": parts[1],
                            "organism": parts[2] if len(parts) > 2 else "",
                            "version": parts[3] if len(parts) > 3 else "",
                        }
                    )

                self._proteome_cache = proteomes
                self._log(f"Found {len(proteomes)} available proteomes")
                return proteomes

        except Exception as e:
            self._log(f"Error retrieving proteome list: {e}")
            return []

    def download_proteome(
        self, organism_code, pdir=None, extract=True, overwrite=False
    ):
        """Download complete proteome for a model organism.

        Downloads bulk proteome archives from the AlphaFold FTP server.
        These contain all predicted structures for organisms like human, mouse, etc.

        :param organism_code: Organism code (e.g., 'HUMAN', 'MOUSE', 'ECOLI')
            or full filename (e.g., 'UP000005640_9606_HUMAN_v6.tar')
        :type organism_code: string

        :param pdir: Directory to save and extract files
        :type pdir: string

        :param extract: Extract tar archive after download
        :type extract: bool

        :param overwrite: Overwrite existing files
        :type overwrite: bool

        :return: Path to downloaded/extracted files
        :rtype: string

        Example:
            >>> afl = AlphaFoldList()
            >>> afl.download_proteome('HUMAN', extract=True)
            Downloading human proteome (4.8 GB)...
        """
        # Get available proteomes
        proteomes = self.get_available_proteomes()

        # Find matching proteome
        proteome = None
        for p in proteomes:
            if organism_code.upper() in p["filename"].upper():
                proteome = p
                break

        if not proteome:
            self._log(f"Organism '{organism_code}' not found in available proteomes")
            self._log("Use get_available_proteomes() to see available options")
            return None

        filename = proteome["filename"]
        url = f"{self.ftp_server}/latest/{filename}"

        # Determine save path
        if pdir is None:
            path = self.local_dir
        else:
            path = pdir

        os.makedirs(path, exist_ok=True)

        tar_file = os.path.join(path, filename)

        # Check if already exists
        if not overwrite and os.path.exists(tar_file):
            self._log(f"Archive exists: '{tar_file}'")
            if extract:
                self._extract_proteome(tar_file, path)
            return path

        # Download the archive
        self._log(f"Downloading {organism_code} proteome: {filename}")
        self._log("This may take a while for large proteomes (several GB)...")

        try:
            urlcleanup()
            urlretrieve(url, tar_file)
            self._log(f"Downloaded: {tar_file}")

            if extract:
                self._extract_proteome(tar_file, path)

            return path

        except Exception as e:
            self._log(f"Download failed: {e}")
            if os.path.exists(tar_file):
                os.remove(tar_file)
            return None

    def _extract_proteome(self, tar_file, extract_path):
        """Extract proteome tar archive (PRIVATE)."""
        self._log(f"Extracting {tar_file}...")
        try:
            with tarfile.open(tar_file, "r") as tar:
                # Extract .gz files
                tar.extractall(path=extract_path)

                # Decompress .gz files
                for member in tar.getmembers():
                    if member.name.endswith(".gz"):
                        gz_path = os.path.join(extract_path, member.name)
                        out_path = gz_path[:-3]  # Remove .gz extension

                        if os.path.exists(gz_path):
                            with gzip.open(gz_path, "rb") as f_in:
                                with open(out_path, "wb") as f_out:
                                    f_out.write(f_in.read())
                            os.remove(gz_path)

            self._log(f"Extraction complete: {extract_path}")
        except Exception as e:
            self._log(f"Extraction failed: {e}")


if __name__ == "__main__":
    doc = """AlphaFoldList.py
    Access and download AlphaFold predicted protein structures

    Usage:
        AlphaFoldList.py <uniprot_id> [options]    - Download single structure
        AlphaFoldList.py proteomes                 - List available proteomes
        AlphaFoldList.py proteome <organism> [dir] - Download organism proteome

    Options:
        -d <dir>    Save to specified directory
        -pdb        Download in PDB format (default: mmCIF)
        -pae        Include predicted aligned error
        -conf       Include confidence scores
        -v <n>      Download specific version number

    Examples:
        AlphaFoldList.py P00520
        AlphaFoldList.py P00520 -pdb -pae -conf
        AlphaFoldList.py proteome HUMAN
    """

    if len(sys.argv) < 2:
        print(doc)
        sys.exit(0)

    afl = AlphaFoldList()

    if sys.argv[1] == "proteomes":
        # List available proteomes
        proteomes = afl.get_available_proteomes()
        print(f"\n{'Organism':<15} {'UniProt ID':<15} {'Tax ID':<10} {'Filename'}")
        print("=" * 80)
        for p in proteomes:
            print(
                f"{p['organism']:<15} {p['uniprot_id']:<15} {p['taxid']:<10} {p['filename']}"
            )

    elif sys.argv[1] == "proteome" and len(sys.argv) >= 3:
        # Download proteome
        organism = sys.argv[2]
        directory = sys.argv[3] if len(sys.argv) > 3 else None
        afl.download_proteome(organism, pdir=directory)

    else:
        # Download single structure
        uniprot_id = sys.argv[1]
        pdir = None
        file_format = "mmCif"
        include_pae = False
        include_conf = False
        version = None

        # Parse options
        i = 2
        while i < len(sys.argv):
            if sys.argv[i] == "-d" and i + 1 < len(sys.argv):
                pdir = sys.argv[i + 1]
                i += 2
            elif sys.argv[i] == "-pdb":
                file_format = "pdb"
                i += 1
            elif sys.argv[i] == "-pae":
                include_pae = True
                i += 1
            elif sys.argv[i] == "-conf":
                include_conf = True
                i += 1
            elif sys.argv[i] == "-v" and i + 1 < len(sys.argv):
                version = int(sys.argv[i + 1])
                i += 2
            else:
                i += 1

        result = afl.retrieve_alphafold_structure(
            uniprot_id,
            pdir=pdir,
            file_format=file_format,
            version=version,
            include_pae=include_pae,
            include_confidence=include_conf,
        )

        if result:
            print(f"\nSuccess! Downloaded to: {result}")
        else:
            print(f"\nFailed to download structure for {uniprot_id}")
            sys.exit(1)
