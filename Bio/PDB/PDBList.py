#!/usr/bin/env python
# Copyright 2003, by Kristian Rother. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# PDBList.py
#
# A tool for tracking changes in the PDB Protein Structure Database.
#
# (c) 2003 Kristian Rother
# This work was supported by the German Ministry of Education
# and Research (BMBF). Project http://www.bcbio.de
#
# Contact the author
#    homepage : http://www.rubor.de/bioinf
#    email    : krother@genesilico.pl
#
#
# (c) 2016 Wiktoria Karwicka & Jacek Smietanski
#   - updated and Python 3.x compatible code
#   - new options to enable download PDBx/mmCif, PDBML and mmtf formatted
#       files as well as large PDB bundles
#   - unit tests for the module
#
# Contact the corresponding author
#   homepage : http://jaceksmietanski.net
#   email    : jacek.smietanski@ii.uj.edu.pl
#
#
# (c) 2022 Julian Maurin
#   - typing and refactoring
#   - introduced PDBFileFormat
#   - handle multiple format moving obsolete file (PDBList.update_pdb)
#
# It may be distributed freely with respect to the original authors.
# Any maintainer of the Biopython code may change this notice
# when appropriate.

"""Access the PDB over the internet (e.g. to download structures)."""

from __future__ import annotations

import contextlib
import enum
import gzip
import os
import pathlib
import re
import shutil
import sys
from urllib.request import urlcleanup
from urllib.request import urlopen
from urllib.request import urlretrieve

MMTF_SERVER_URL = "mmtf.rcsb.org"


class PDBListError(Exception):
    """Generic exception for PDBList module."""


@enum.unique
class PDBFileFormat(enum.Enum):
    """Enum regrouping file format information."""

    PDB = ("pdb", "pdb", "pdb", ".ent")  # PDBx/mmCif (default)
    MMCIF = ("mmCif", "mmCIF", "", ".cif")  # PDB
    XML = ("xml", "XML", "", ".xml")  # PDBML/XML
    MMTF = ("mmtf", "", "", ".mmtf")  # highly compressed
    BUNDLE = (
        "bundle",
        "pdb_bundle",
        "",
        "-pdb-bundle.tar",
    )  # PDB formatted archive for large structure

    @property
    def label(self) -> str:
        """Allow backward compatibility with file_format as string.

        NOTE: To be remove once the last usage of file_format as string is updated.
        """
        return self.value[0]

    @property
    def directory(self) -> str:
        """Name of the directory where the file is stored on the server."""
        return self.value[1]

    @property
    def filename_prefix(self) -> str:
        """Return the string to add at the beginning of the code to build the filename."""
        return self.value[2]

    @property
    def filename_suffix(self) -> str:
        """Return the string to add at the end of the code to build the filename."""
        return self.value[3]

    @classmethod
    def from_label(cls, label) -> PDBFileFormat:
        """Get PDBFileFormat from label."""
        for file_format in cls:
            if file_format.label == label:
                return file_format
        raise PDBListError(
            "Specified file_format does not exist or is not supported, maybe a typo "
            f" (file format: {label}, handled file formats: {','.join([file_format.label for file_format in PDBFileFormat])}."
        )


class PDBList:
    """Quick access to the structure lists on the PDB or its mirrors.

    This class provides quick access to the structure lists on the
    PDB server or its mirrors. The structure lists contain
    four-letter PDB codes, indicating that structures are
    new, have been modified or are obsolete. The lists are released
    on a weekly basis.

    It also provides a function to retrieve PDB files from the server.
    To use it properly, prepare a directory /pdb or the like,
    where PDB files are stored.

    All available file formats (PDB, PDBx/mmCif, PDBML, mmtf) are supported.
    Please note that large structures (containing >62 chains
    and/or 99999 ATOM lines) are no longer stored as a single PDB file
    and by default (when PDB format selected) are not downloaded.

    Large structures can be downloaded in other formats, including PDBx/mmCif
    or as a .tar file (a collection of PDB-like formatted files for a given
    structure).

    If you want to use this module from inside a proxy, add
    the proxy variable to your environment, e.g. in Unix:
    export HTTP_PROXY='http://realproxy.charite.de:888'
    (This can also be added to ~/.bashrc)
    """

    PDB_REF = """
    The Protein Data Bank: a computer-based archival file for macromolecular structures.
    F.C.Bernstein, T.F.Koetzle, G.J.B.Williams, E.F.Meyer Jr, M.D.Brice, J.R.Rodgers, O.Kennard, T.Shimanouchi, M.Tasumi
    J. Mol. Biol. 112 pp. 535-542 (1977)
    http://www.pdb.org/.
    """

    def __init__(
        self,
        server: str = "ftp://ftp.wwpdb.org",
        pdb: str | None = None,
        obsolete_pdb: str | None = None,
        verbose: bool = True,
    ):
        """Initialize the class with the default server or a custom one.

        Argument pdb is the local path to use, defaulting to the current
        directory at the moment of initialisation.
        """
        self.pdb_server = server  # remote pdb server

        self.local_pdb = pdb or os.getcwd()  # local pdb file tree
        pathlib.Path(self.local_pdb).mkdir(parents=True, exist_ok=True)

        self.obsolete_pdb = (
            obsolete_pdb if obsolete_pdb else os.path.join(self.local_pdb, "obsolete")
        )  # local file tree for obsolete pdb files
        pathlib.Path(self.obsolete_pdb).mkdir(parents=True, exist_ok=True)

        # enable or disable verbose
        self._verbose = verbose

        # variable for command-line option
        self.flat_tree = False

    @staticmethod
    def _print_default_format_warning(file_format):
        """Print a warning to stdout (PRIVATE).

        Temporary warning (similar to a deprecation warning) that files
        are being downloaded in mmCIF.
        """
        if file_format is None:
            sys.stderr.write(
                "WARNING: The default download format has changed from PDB to PDBx/mmCif\n"
            )
            return PDBFileFormat.MMCIF.label
        return file_format

    @staticmethod
    def get_status_list(url: str) -> list[str]:
        """Retrieve a list of pdb codes in the weekly pdb status file from given URL.

        Used by get_recent_changes. Typical contents of the list files parsed
        by this method is now very simply - one PDB name per line.
        """
        pdb_codes = []
        with contextlib.closing(urlopen(url)) as handle:
            for line in handle:
                pdb_code = line.strip()
                if len(pdb_code) != 4:
                    raise PDBListError(
                        f"Status list contains unexpected PDB code value (url: {url}, value: {pdb_code})."
                    )
                pdb_codes.append(pdb_code.decode())
        return pdb_codes

    def get_recent_changes(self) -> tuple[list[str], list[str], list[str]]:
        """Return three lists of the newest weekly files (added,mod,obsolete).

        Reads the directories with changed entries from the PDB server and
        returns a tuple of three URL's to the files of new, modified and
        obsolete entries from the most recent list. The directory with the
        largest numerical name is used.
        Returns None if something goes wrong.

        Contents of the data/status dir (20031013 would be used);:

            drwxrwxr-x   2 1002     sysadmin     512 Oct  6 18:28 20031006
            drwxrwxr-x   2 1002     sysadmin     512 Oct 14 02:14 20031013
            -rw-r--r--   1 1002     sysadmin    1327 Mar 12  2001 README

        """
        path = self.pdb_server + "/pub/pdb/data/status/latest/"

        # Retrieve the lists
        added = self.get_status_list(path + "added.pdb")
        modified = self.get_status_list(path + "modified.pdb")
        obsolete = self.get_status_list(path + "obsolete.pdb")
        return (added, modified, obsolete)

    def get_all_entries(self) -> list[str]:
        """Retrieve the big file containing all the PDB entries and some annotation.

        Returns a list of PDB codes in the index file.
        """
        url = self.pdb_server + "/pub/pdb/derived_data/index/entries.idx"
        if self._verbose:
            print("Retrieving index file. Takes about 27 MB.")
        with contextlib.closing(urlopen(url)) as handle:
            all_entries = [
                line[:4].decode() for line in handle.readlines()[2:] if len(line) > 4
            ]
        return all_entries

    def get_all_obsolete(self) -> list[str]:
        """Return a list of all obsolete entries ever in the PDB.

        Returns a list of all obsolete pdb codes that have ever been
        in the PDB.

        Gets and parses the file from the PDB server in the format
        (the first pdb_code column is the one used). The file looks
        like this::

             LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS
            OBSLTE    31-JUL-94 116L     216L
            ...
            OBSLTE    29-JAN-96 1HFT     2HFT
            OBSLTE    21-SEP-06 1HFV     2J5X
            OBSLTE    21-NOV-03 1HG6
            OBSLTE    18-JUL-84 1HHB     2HHB 3HHB
            OBSLTE    08-NOV-96 1HID     2HID
            OBSLTE    01-APR-97 1HIU     2HIU
            OBSLTE    14-JAN-04 1HKE     1UUZ
            ...

        """
        url = self.pdb_server + "/pub/pdb/data/status/obsolete.dat"
        obsolete_pdb_codes = []
        with contextlib.closing(urlopen(url)) as handle:
            # Extract pdb codes. Could use a list comprehension, but I want
            # to raise an Exception when mis-reading the data.
            for line in handle:
                if not line.startswith(b"OBSLTE "):
                    continue
                pdb_code = line.split()[2]
                if len(pdb_code) != 4:
                    raise PDBListError(
                        f"Obsolete list contains unexpected PDB code value (url: {url}, value: {pdb_code})."
                    )
                obsolete_pdb_codes.append(pdb_code.decode())
        return obsolete_pdb_codes

    def retrieve_pdb_file(
        self,
        pdb_code: str,
        obsolete: bool = False,
        pdir: str | None = None,
        file_format: str | PDBFileFormat | None = None,
        overwrite: bool = False,
    ) -> str:
        """Fetch PDB structure file from PDB server, and store it locally.

        The PDB structure's file name is returned as a single string.
        If obsolete ``==`` True, the file will be saved in a special file tree.

        NOTE. The default download format has changed from PDB to PDBx/mmCif

        :param pdb_code: 4-symbols structure Id from PDB (e.g. 3J92).
        :type pdb_code: string

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PDBML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure}

        :type file_format: string | PDBFileFormat

        :param overwrite: if set to True, existing structure files will be overwritten. Default: False
        :type overwrite: bool

        :param obsolete:
            Has a meaning only for obsolete structures. If True, download the obsolete structure
            to 'obsolete' folder, otherwise download won't be performed.
            This option doesn't work for mmtf format as obsoleted structures aren't stored in mmtf.
            Also doesn't have meaning when parameter pdir is specified.
            Note: make sure that you are about to download the really obsolete structure.
            Trying to download non-obsolete structure into obsolete folder will not work
            and you face the "structure doesn't exists" error.
            Default: False

        :type obsolete: bool

        :param pdir: put the file in this directory (default: create a PDB-style directory tree)
        :type pdir: string

        :return: filename
        :rtype: string
        """
        if isinstance(file_format, PDBFileFormat):
            file_format_enum = file_format
        else:
            file_format = self._print_default_format_warning(file_format)
            file_format_enum = PDBFileFormat.from_label(file_format)

        code = pdb_code.lower()
        short_code = code[1:3]

        if file_format_enum != PDBFileFormat.MMTF:
            filename = f"{file_format_enum.filename_prefix}{code}{file_format_enum.filename_suffix}"
            archive_filename = f"{filename}.gz"
        else:
            filename = archive_filename = code

        output_directory = self.get_output_directory(pdir, obsolete, short_code)
        output_archive_filepath = pathlib.Path(output_directory, archive_filename)
        output_extracted_filename = (
            f"{filename}{file_format_enum.filename_suffix}"
            if file_format_enum == PDBFileFormat.MMTF
            else filename
        )
        output_extracted_filepath = pathlib.Path(
            output_directory, output_extracted_filename
        )

        # Skip download if the file already exists
        if not overwrite and output_extracted_filepath.exists():
            if self._verbose:
                print(
                    f"Structure exists, skip download and extract (path: {output_extracted_filepath}')."
                )
            return str(output_extracted_filepath)

        # Retrieve the file
        if self._verbose:
            print(f"Downloading and extracting PDB structure (code: {pdb_code}).")
        archive_url = self.build_archive_url(
            file_format_enum, obsolete, code, short_code, archive_filename
        )
        self.download_and_extract_archive(
            archive_url, output_archive_filepath, output_extracted_filepath
        )
        return str(output_extracted_filepath)

    def get_output_directory(
        self, output_directory: str | None, obsolete: bool, short_code: str
    ) -> pathlib.Path:
        """Get the output directory path, creating it if does not exist."""
        if output_directory is None:
            path = pathlib.Path(self.obsolete_pdb if obsolete else self.local_pdb)
            if not self.flat_tree:  # Put in PDB-style directory tree
                path = pathlib.Path(path, short_code)
        else:
            path = pathlib.Path(output_directory)

        path.mkdir(parents=True, exist_ok=True)
        return path

    def build_archive_url(
        self,
        file_format: PDBFileFormat,
        obsolete: bool,
        code: str,
        short_code: str,
        archive_filename: str,
    ) -> str:
        """Build archive URL according to the file format."""
        if file_format in (PDBFileFormat.PDB, PDBFileFormat.MMCIF, PDBFileFormat.XML):
            return self.pdb_server + "/pub/pdb/data/structures/%s/%s/%s/%s" % (
                "divided" if not obsolete else "obsolete",
                file_format.directory,
                short_code,
                archive_filename,
            )
        elif file_format == PDBFileFormat.BUNDLE:
            return self.pdb_server + "/pub/pdb/compatible/%s/%s/%s/%s" % (
                file_format.directory,
                short_code,
                code,
                archive_filename,
            )
        elif file_format == PDBFileFormat.MMTF:
            return f"http://{MMTF_SERVER_URL}/v1.0/full/{code}"
        raise PDBListError(f"Unhandled file format (format: {file_format.label}).")

    def download_and_extract_archive(
        self,
        archive_url: str,
        output_archive_filepath: pathlib.Path,
        output_extracted_filepath: pathlib.Path,
    ) -> None:
        """
        Download archive from server and extract it to the output directory.

        The archive is removed once extracted.
        """
        try:
            urlcleanup()
            urlretrieve(archive_url, str(output_archive_filepath))
        except OSError:
            print(f"PDB file not found on remote server (url: {archive_url}).")
        else:
            with gzip.open(output_archive_filepath, "rb") as gzip_stream:
                with open(output_extracted_filepath, "wb") as output_stream:
                    output_stream.writelines(gzip_stream)
            output_archive_filepath.unlink()

    def update_pdb(self, file_format=None):
        """Update your local copy of the PDB files.

        I guess this is the 'most wanted' function from this module.
        It gets the weekly lists of new and modified pdb entries and
        automatically downloads the according PDB files.
        You can call this module as a weekly cron job.
        """
        if isinstance(file_format, PDBFileFormat):
            file_format_enum = file_format
        else:
            file_format = self._print_default_format_warning(file_format)
            file_format_enum = PDBFileFormat.from_label(file_format)

        new, modified, obsolete = self.get_recent_changes()
        for pdb_code in new + modified:
            try:
                self.retrieve_pdb_file(
                    pdb_code,
                    file_format=file_format_enum,
                    overwrite=pdb_code in modified,
                )
            except Exception as err:
                print(
                    f"Error retrieving pdb file (code: {pdb_code}, exception: {err.__class__.__name__}).\n"
                )

        # Move the obsolete files to a special folder
        if file_format in (PDBFileFormat.PDB, PDBFileFormat.MMCIF, PDBFileFormat.XML):
            for pdb_code in obsolete:
                self.move_obsolete_file(pdb_code, file_format)

    def move_obsolete_file(self, pdb_code: str, file_format: PDBFileFormat) -> None:
        """Move obsolete file to obsolete directory."""
        short_code = pdb_code[1:3]
        filename = (
            f"{file_format.filename_prefix}{pdb_code}{file_format.filename_suffix}"
        )

        if self.flat_tree:
            old_file = pathlib.Path(self.local_pdb, filename)
            new_dir = pathlib.Path(self.obsolete_pdb)
        else:
            old_file = pathlib.Path(self.local_pdb, short_code, filename)
            new_dir = pathlib.Path(self.obsolete_pdb, short_code)

        new_file = pathlib.Path(new_dir, filename)

        if old_file.is_file():
            new_dir.mkdir(parents=True, exist_ok=True)
            try:
                shutil.move(old_file, new_file)
            except Exception:
                print(
                    f"Could not move file to obsolete directory (source path: {old_file}, destination path: {new_file})."
                )
        elif new_file.is_file():
            if self._verbose:
                print(f"Obsolete file already moved (path: {old_file}).")
        else:
            if self._verbose:
                print(f"Obsolete file is missing (path: {old_file}).")

    def download_pdb_files(
        self, pdb_codes, obsolete=False, pdir=None, file_format=None, overwrite=False
    ):
        """Fetch set of PDB structure files from the PDB server and stores them locally.

        The PDB structure's file name is returned as a single string.
        If obsolete ``==`` True, the files will be saved in a special file tree.

        :param pdb_codes: a list of 4-symbols structure Ids from PDB
        :type pdb_codes: list of strings

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PMDML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure}

        :param overwrite: if set to True, existing structure files will be overwritten. Default: False
        :type overwrite: bool

        :param obsolete:
            Has a meaning only for obsolete structures.
            If True, download the obsolete structure
            to 'obsolete' folder, otherwise download won't be performed.
            This option doesn't work for mmtf format as obsoleted structures are not available as mmtf.
            (default: False)

        :type obsolete: bool

        :param pdir: put the file in this directory (default: create a PDB-style directory tree)
        :type pdir: string

        :return: filenames
        :rtype: string
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        for pdb_code in pdb_codes:
            self.retrieve_pdb_file(
                pdb_code,
                obsolete=obsolete,
                pdir=pdir,
                file_format=file_format,
                overwrite=overwrite,
            )

    def download_entire_pdb(self, listfile=None, file_format=None):
        """Retrieve all PDB entries not present in the local PDB copy.

        :param listfile: filename to which all PDB codes will be written (optional)

        :param file_format:
            File format. Available options:

            * "mmCif" (default, PDBx/mmCif file),
            * "pdb" (format PDB),
            * "xml" (PMDML/XML format),
            * "mmtf" (highly compressed),
            * "bundle" (PDB formatted archive for large structure}

        NOTE. The default download format has changed from PDB to PDBx/mmCif
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        entries = self.get_all_entries()
        for pdb_code in entries:
            self.retrieve_pdb_file(pdb_code, file_format=file_format)
        # Write the list
        if listfile:
            with open(listfile, "w") as outfile:
                outfile.writelines(x + "\n" for x in entries)

    def download_obsolete_entries(self, listfile=None, file_format=None):
        """Retrieve all obsolete PDB entries not present in local obsolete PDB copy.

        :param listfile: filename to which all PDB codes will be written (optional)

        :param file_format: file format. Available options:
            "mmCif" (default, PDBx/mmCif file),
            "pdb" (format PDB),
            "xml" (PMDML/XML format),

        NOTE. The default download format has changed from PDB to PDBx/mmCif
        """
        # Deprecation warning
        file_format = self._print_default_format_warning(file_format)
        entries = self.get_all_obsolete()
        for pdb_code in entries:
            self.retrieve_pdb_file(pdb_code, obsolete=True, file_format=file_format)

        # Write the list
        if listfile:
            with open(listfile, "w") as outfile:
                outfile.writelines(x + "\n" for x in entries)

    def get_seqres_file(self, savefile="pdb_seqres.txt"):
        """Retrieve and save a (big) file containing all the sequences of PDB entries."""
        if self._verbose:
            print("Retrieving sequence file (takes over 110 MB).")
        url = self.pdb_server + "/pub/pdb/derived_data/pdb_seqres.txt"
        urlretrieve(url, savefile)


if __name__ == "__main__":

    doc = """PDBList.py
    (c) Kristian Rother 2003, Wiktoria Karwicka & Jacek Smietanski 2016
    Contributed to Biopython

    Usage::

        PDBList.py update <pdb_path> [options]   - write weekly PDB updates to
                                                   local pdb tree.
        PDBList.py all    <pdb_path> [options]   - write all PDB entries to
                                                   local pdb tree.
        PDBList.py obsol  <pdb_path> [options]   - write all obsolete PDB
                                                   entries to local pdb tree.
        PDBList.py <PDB-ID> <pdb_path> [options] - retrieve single structure
        PDBList.py (<PDB-ID1>,<PDB-ID2>,...) <pdb_path> [options] - retrieve a set
                                                   of structures

    Options:
     -d       A single directory will be used as <pdb_path>, not a tree.
     -o       Overwrite existing structure files.
     -pdb     Downloads structures in PDB format
     -xml     Downloads structures in PDBML (XML) format
     -mmtf    Downloads structures in mmtf format

    Maximum one format can be specified simultaneously (if more selected, only
    the last will be considered). By default (no format specified) structures are
    downloaded as PDBx/mmCif files.
    """
    print(doc)

    file_format = "mmCif"
    overwrite = False

    if len(sys.argv) > 2:
        pdb_path = sys.argv[2]
        pl = PDBList(pdb=pdb_path)
        if len(sys.argv) > 3:
            for option in sys.argv[3:]:
                if option == "-d":
                    pl.flat_tree = True
                elif option == "-o":
                    overwrite = True
                elif option in ("-pdb", "-xml", "-mmtf"):
                    file_format = option[1:]
    else:
        pdb_path = os.getcwd()
        pl = PDBList()
        pl.flat_tree = True

    if len(sys.argv) > 1:
        if sys.argv[1] == "update":
            # update PDB
            print("updating local PDB at " + pdb_path)
            pl.update_pdb(file_format=file_format)

        elif sys.argv[1] == "all":
            # get the entire PDB
            pl.download_entire_pdb(file_format=file_format)

        elif sys.argv[1] == "obsol":
            # get all obsolete entries
            pl.download_obsolete_entries(pdb_path, file_format=file_format)

        elif len(sys.argv[1]) == 4 and sys.argv[1][0].isdigit():
            # get single PDB entry
            pl.retrieve_pdb_file(
                sys.argv[1], pdir=pdb_path, file_format=file_format, overwrite=overwrite
            )

        elif sys.argv[1][0] == "(":
            # get a set of PDB entries
            pdb_ids = re.findall("[0-9A-Za-z]{4}", sys.argv[1])
            for pdb_id in pdb_ids:
                pl.retrieve_pdb_file(
                    pdb_id, pdir=pdb_path, file_format=file_format, overwrite=overwrite
                )
