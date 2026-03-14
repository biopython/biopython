"""A module for interacting with the AlphaFold Protein Structure Database.

See the `database website <https://alphafold.ebi.ac.uk/>`_ and the
`API docs <https://alphafold.ebi.ac.uk/api-docs/>`_.
"""

import json
import os
from os import PathLike
from collections.abc import Iterator
from typing import Optional
from typing import Union
from urllib.request import urlopen
from urllib.request import urlretrieve
from urllib.request import urlcleanup

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure as StructuralModel


def get_predictions(qualifier: str) -> Iterator[dict]:
    """Get all AlphaFold predictions for a UniProt accession.

    :param qualifier: A UniProt accession, e.g. P00520
    :type qualifier: str
    :return: The AlphaFold predictions
    :rtype: Iterator[dict]
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{qualifier}"
    # Retrieve the AlphaFold predictions with urllib
    with urlopen(url) as response:
        yield from json.loads(response.read().decode())


def _get_mmcif_file_path_for(
    prediction: dict, directory: str | bytes | PathLike | None = None
) -> str:
    """Get the path to the mmCIF file for an AlphaFold prediction.

    :param prediction: An AlphaFold prediction
    :type prediction: dict
    :param directory: The directory that stores the mmCIF file, defaults to the current working directory
    :type directory: Union[int, str, bytes, PathLike], optional
    :return: The path to the mmCIF file
    :rtype: str
    """
    if directory is None:
        directory = os.getcwd()

    cif_url = prediction["cifUrl"]
    # Get the file name from the URL
    file_name = cif_url.split("/")[-1]
    return str(os.path.join(directory, file_name))


def download_cif_for(
    prediction: dict, directory: str | bytes | PathLike | None = None
) -> str:
    """Download the mmCIF file for an AlphaFold prediction.

    Downloads the file to the current working directory if no destination is specified.

    :param prediction: An AlphaFold prediction
    :type prediction: dict
    :param directory: The directory to write the mmCIF data to, defaults to the current working directory
    :type directory: Union[int, str, bytes, PathLike], optional
    :return: The path to the mmCIF file
    :rtype: str
    """
    if directory is None:
        directory = os.getcwd()

    cif_url = prediction["cifUrl"]
    # Create the directory in case it does not exist
    os.makedirs(directory, exist_ok=True)
    file_path = _get_mmcif_file_path_for(prediction, directory)

    if os.path.exists(file_path):
        print(f"File {file_path} already exists, skipping download.")
    else:
        with urlopen(cif_url) as response:
            data = response.read()
        # Write the data to destination
        with open(file_path, "wb") as file:
            file.write(data)

    return file_path


def get_structural_models_for(
    qualifier: str,
    mmcif_parser: MMCIFParser | None = None,
    directory: str | bytes | PathLike | None = None,
) -> Iterator[StructuralModel]:
    """Get the PDB structures for a UniProt accession.

    Downloads the mmCIF files to the directory if they are not present.

    :param qualifier: A UniProt accession, e.g. P00520
    :type qualifier: str
    :param mmcif_parser: The mmCIF parser to use, defaults to ``MMCIFParser()``
    :type mmcif_parser: MMCIFParser, optional
    :param directory: The directory to store the mmCIF data, defaults to the current working directory
    :type directory: Union[int, str, bytes, PathLike], optional
    :return: An iterator over the PDB structures
    :rtype: Iterator[PDBStructure]
    """
    if mmcif_parser is None:
        mmcif_parser = MMCIFParser()
    if directory is None:
        directory = os.getcwd()

    for prediction in get_predictions(qualifier):
        mmcif_path = _get_mmcif_file_path_for(prediction, directory)

        if not os.path.exists(mmcif_path):
            mmcif_path = download_cif_for(prediction, directory)

        yield mmcif_parser.get_structure(qualifier, mmcif_path)


class AFDBList:
    """Download AlphaFold structures by UniProt ID, similar to PDBList for PDB.

    This class provides a PDBList-like interface for the AlphaFold Protein
    Structure Database. You can download predicted structures for any UniProt
    accession that has an AlphaFold model.

    Example::

        from Bio.PDB import AFDBList

        afdb = AFDBList()
        # Download the canonical model for P00520 (ABL1) as mmCIF
        path = afdb.retrieve_pdb_file("P00520")
        # Or as PDB format
        path = afdb.retrieve_pdb_file("P00520", file_format="pdb")
    """

    def __init__(
        self,
        server: str = "https://alphafold.ebi.ac.uk",
        pdb: str | None = None,
        verbose: bool = True,
    ):
        """Set up the download interface.

        :param server: Base URL for the AlphaFold API (default: EBI)
        :param pdb: Local directory to store files (default: current directory)
        :param verbose: Whether to print progress messages
        """
        self.server = server.rstrip("/")
        self.local_pdb = pdb if pdb else os.getcwd()
        self._verbose = verbose

    def retrieve_pdb_file(
        self,
        uniprot_id: str,
        pdir: str | None = None,
        file_format: str | None = None,
        overwrite: bool = False,
        model_index: int = 0,
    ) -> str | None:
        """Download an AlphaFold structure for a UniProt accession.

        :param uniprot_id: UniProt accession (e.g. P00520, Q9Y2X3)
        :param pdir: Directory to save the file (default: local_pdb)
        :param file_format: "mmCif" (default) or "pdb"
        :param overwrite: If True, re-download even if file exists
        :param model_index: Which model to use when multiple exist (0 = canonical)
        :return: Path to the downloaded file, or None if not found
        """
        if file_format is None:
            file_format = "mmCif"

        if file_format not in ("mmCif", "pdb"):
            raise ValueError(
                f"file_format must be 'mmCif' or 'pdb', got {file_format!r}"
            )

        url = f"{self.server}/api/prediction/{uniprot_id}"
        try:
            with urlopen(url) as response:
                predictions = json.loads(response.read().decode())
        except OSError as e:
            if self._verbose:
                print(f"Could not fetch predictions for {uniprot_id}: {e}")
            return None

        if not predictions:
            if self._verbose:
                print(f"No AlphaFold model found for {uniprot_id}")
            return None

        # Prefer the canonical model (exact uniprot_id match) when available
        canonical = [
            p for p in predictions if p.get("uniprotAccession") == uniprot_id
        ]
        prediction_list = canonical if canonical else predictions
        try:
            prediction = prediction_list[model_index]
        except IndexError:
            raise ValueError(
                f"model_index {model_index} out of range "
                f"(max {len(prediction_list) - 1})"
            )

        if file_format == "mmCif":
            download_url = prediction["cifUrl"]
            filename = download_url.split("/")[-1]
        else:
            download_url = prediction["pdbUrl"]
            filename = download_url.split("/")[-1]

        path = pdir if pdir else self.local_pdb
        os.makedirs(path, exist_ok=True)
        filepath = os.path.join(path, filename)

        if os.path.exists(filepath) and not overwrite:
            if self._verbose:
                print(f"File already exists: {filepath}")
            return filepath

        if self._verbose:
            print(f"Downloading AlphaFold model for {uniprot_id}...")

        try:
            urlcleanup()
            urlretrieve(download_url, filepath)
        except OSError as e:
            if self._verbose:
                print(f"Download failed for {uniprot_id}: {e}")
            return None

        return filepath
