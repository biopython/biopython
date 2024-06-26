"""A module for interacting with the AlphaFold Protein Structure Database.

See the `database website <https://alphafold.com/>`_ and the `API docs <https://alphafold.com/api-docs/>`_.
"""

import json
import os
from os import PathLike
from collections.abc import Iterator
from typing import Optional
from typing import Union
from urllib.request import urlopen

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Structure import Structure as StructuralModel


def get_predictions(qualifier: str) -> Iterator[dict]:
    """Get all AlphaFold predictions for a UniProt accession.

    :param qualifier: A UniProt accession, e.g. P00520
    :type qualifier: str
    :return: The AlphaFold predictions
    :rtype: Iterator[dict]
    """
    url = f"https://alphafold.com/api/prediction/{qualifier}"
    # Retrieve the AlphaFold predictions with urllib
    with urlopen(url) as response:
        yield from json.loads(response.read().decode())


def _get_mmcif_file_path_for(
    prediction: dict, directory: Optional[Union[str, bytes, PathLike]] = None
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
    prediction: dict, directory: Optional[Union[str, bytes, PathLike]] = None
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
    mmcif_parser: Optional[MMCIFParser] = None,
    directory: Optional[Union[str, bytes, PathLike]] = None,
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
