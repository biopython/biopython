"""A module for interacting with the AlphaFold Protein Structure Database.

See the `database website <https://alphafold.com/>`_ and the `API docs <https://alphafold.com/api-docs/>`_.
"""

from functools import lru_cache
import json
import os
import gzip
import tarfile
from os import PathLike
from collections.abc import Iterator
#from typing import Optional
#from typing import Union
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

@lru_cache(maxsize=1)
def list_proteomes() -> list[dict]:
    """List all available proteomes in the AlphaFold Database.

    Downloads the metadata index from the EBI FTP server and returns a list
    of proteome entries. Each entry is a dict with keys such as ``species``,
    ``common_name``, ``reference_proteome``, ``archive_name``,
    ``num_predicted_structures``, and ``size_bytes``.

    :return: Available proteome entries
    :rtype: list[dict]
    """
    url = "https://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json"
    with urlopen(url) as response:
        return json.loads(response.read().decode())
    
    

def download_proteome(species: str, directory: str | bytes| PathLike | None = None,
    ) -> str: 
    """Download all AlphaFold predictions for a reference proteome.

    Downloads the proteome archive (.tar) from the EMBL-EBI FTP server.
    The species can be matched by scientific name, common name, or
    UniProt reference proteome ID.

    Example usage::

        path = download_proteome("Homo sapiens", directory="~/alphafold")
        path = download_proteome("Human", directory="~/alphafold")
        path = download_proteome("UP000005640", directory="~/alphafold")

    :param species: Species identifier - scientific name (e.g. 'Homo sapiens'),
        common name (e.g. 'Human'), or UniProt proteome ID (e.g. 'UP000005640')
    :type species: str
    :param directory: The directory to write the archive to, defaults to
        the current working directory
    :type directory: Union[str, bytes, PathLike], optional
    :return: The path to the downloaded tar archive
    :rtype: str
    :raises ValueError: If the species cannot be found in the AlphaFold database
    """

    if directory is None:
        directory = os.getcwd()

    #find the matching proteome entry
    proteomes = list_proteomes()
    match = None
    species_lower = species.lower()
    for proteome in proteomes:
        if proteome.get("type") not in ("proteome", "global_health"):
            continue # this skips Swiss-Prot and other non species entries
        if species_lower in (
            proteome["species"].lower(),
            proteome["common_name"].lower(),
            proteome["reference_proteome"].lower(),
        ):
            match = proteome
            break
    if match is None:
        raise ValueError(f"Could not find proteome for '{species}'." f"Use list_proteomes() to see available species.")

    #warn user about download size 
    size_gb = match["size_bytes"] / 1e9
    print(f"Downloading {match['num_predicted_structures']} structures "
        f"for {match['species']} ({size_gb:.1f} GB)...")

    os.makedirs(directory, exist_ok=True)
    archive_name = match["archive_name"]
    file_path = os.path.join(directory, archive_name)

    if os.path.exists(file_path):
        print(f"Archive {file_path} already exists, skipping download.")
        return str(file_path)

    url = f"https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/{archive_name}"
    try:
        with urlopen(url) as response:
            with open(file_path, "wb") as out_file:
                while chunk := response.read(8192):
                    out_file.write(chunk)
    except Exception:
        if os.path.exists(file_path):
            os.remove(file_path)
        raise

    print(f"Downloaded to {file_path}")
    return str(file_path)


def extract_proteome(
    archive_path: str | bytes | PathLike,
    directory: str | bytes | PathLike | None = None,
    file_format: str = "mmCif",
) -> list[str]:
    """Extract an AlphaFold proteome archive downloaded with download_proteome().

    Extracts and decompresses structures of the requested format only.

    :param archive_path: Path to the downloaded .tar archive
    :type archive_path: Union[str, bytes, PathLike]
    :param directory: Directory to extract files to, defaults to the same
        directory as the archive
    :type directory: Union[str, bytes, PathLike], optional
    :param file_format: Structure file format to extract, either 'mmCif'
        or 'pdb', defaults to 'mmCif'
    :type file_format: str
    :return: List of paths to the extracted files
    :rtype: list[str]
    :raises ValueError: If an invalid file_format is specified
    """
    fmt_map = {"mmCif": ".cif.gz", "pdb": ".pdb.gz"}
    if file_format not in fmt_map:
        raise ValueError(
            f"Invalid file_format '{file_format}'. Choose 'mmCif' or 'pdb'."
        )
    extension = fmt_map[file_format]

    if directory is None:
        directory = os.path.dirname(archive_path)
    os.makedirs(directory, exist_ok=True)

    extracted_paths = []
    with tarfile.open(archive_path) as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(extension)]
        print(f"Extracting {len(members)} {file_format} files...")
        for member in members:
            # Decompress and write without the .gz extension
            out_name = member.name.replace(".gz", "")
            out_path = os.path.join(directory, out_name)
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with tar.extractfile(member) as gz_file:
                with gzip.open(gz_file) as decompressed:
                    with open(out_path, "wb") as out_file:
                        out_file.write(decompressed.read())
            extracted_paths.append(str(out_path))

    print(f"Extracted {len(extracted_paths)} files to {directory}")
    return extracted_paths


class AlphaFoldDB:
    """A PDBList-like interface for the AlphaFold Protein Structure Database.

    Provides a high-level, database-style API for bulk downloading and local
    management of AlphaFold predictions.  Inspired by
    ``Bio.PDB.PDBList.PDBList``.

    Basic usage::

        from Bio.PDB.alphafold_db import AlphaFoldDB

        afdb = AlphaFoldDB(directory="~/alphafold")

        # Download an entire human proteome
        afdb.download_proteome("Homo sapiens")

        # Extract only mmCIF files
        afdb.extract_proteome("~/alphafold/UP000005640_9606_HUMAN.tar")

        # Fetch and parse structures for a specific protein (uses local cache)
        for structure in afdb.get_structural_models_for("P00520"):
            print(structure)
    """

    def __init__(
        self,
        directory: str | bytes | PathLike | None = None,
        mmcif_parser: MMCIFParser | None = None,
    ) -> None:
        """Initialize AlphaFoldDB.

        :param directory: Base directory for downloaded files, defaults to
            the current working directory
        :type directory: Union[str, bytes, PathLike], optional
        :param mmcif_parser: Pre-configured MMCIFParser instance
        :type mmcif_parser: MMCIFParser, optional
        """
        self.directory = directory
        self.mmcif_parser = mmcif_parser

    def list_proteomes(self) -> list[dict]:
        """Wraps :func:`list_proteomes`."""
        return list_proteomes()

    def download_proteome(
        self,
        species: str,
        directory: str | bytes | PathLike | None = None,
    ) -> str:
        """Wraps :func:`download_proteome`.

        If *directory* is not given, uses the instance default.
        """
        if directory is None:
            directory = self.directory
        return download_proteome(species, directory=directory)

    def extract_proteome(
        self,
        archive_path: str | bytes | PathLike,
        directory: str | bytes | PathLike | None = None,
        file_format: str = "mmCif",
    ) -> list[str]:
        """Wraps :func:`extract_proteome`."""
        return extract_proteome(
            archive_path, directory=directory, file_format=file_format
        )

    def get_structural_models_for(
        self,
        qualifier: str,
        directory: str | bytes | PathLike | None = None,
    ) -> Iterator[StructuralModel]:
        """Wraps :func:`get_structural_models_for`.

        Uses the instance default directory and mmcif_parser when available.
        """
        if directory is None:
            directory = self.directory
        parser = self.mmcif_parser if self.mmcif_parser is not None else MMCIFParser()
        dir_str = str(directory) if directory is not None else None
        return get_structural_models_for(
            qualifier,
            mmcif_parser=parser,
            directory=dir_str,
        )
    