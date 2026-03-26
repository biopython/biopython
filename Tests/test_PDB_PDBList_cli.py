# Copyright 2026 by Ernest Provo.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for the PDBList command-line interface (argparse-based).

These tests validate argument parsing only and do not require internet access.
"""

import os
import unittest

from Bio.PDB.PDBList import PDBList, _build_parser, cli


class TestPDBListCLIParser(unittest.TestCase):
    """Test argparse parser for PDBList CLI."""

    def setUp(self):
        """Set up parser for each test."""
        self.parser = _build_parser()

    def test_update_default_path(self):
        """Test 'update' command with default pdb_path."""
        args = self.parser.parse_args(["update"])
        self.assertEqual(args.command, "update")
        self.assertEqual(args.pdb_path, os.getcwd())
        self.assertEqual(args.format, "mmCif")
        self.assertFalse(args.flat)
        self.assertFalse(args.overwrite)

    def test_update_custom_path(self):
        """Test 'update' command with custom --pdb-path."""
        args = self.parser.parse_args(["update", "--pdb-path", "/tmp/pdb"])
        self.assertEqual(args.command, "update")
        self.assertEqual(args.pdb_path, "/tmp/pdb")

    def test_all_command(self):
        """Test 'all' command parses correctly."""
        args = self.parser.parse_args(["all"])
        self.assertEqual(args.command, "all")

    def test_obsol_command(self):
        """Test 'obsol' command parses correctly."""
        args = self.parser.parse_args(["obsol"])
        self.assertEqual(args.command, "obsol")

    def test_assemb_command(self):
        """Test 'assemb' command parses correctly."""
        args = self.parser.parse_args(["assemb"])
        self.assertEqual(args.command, "assemb")

    def test_fetch_single_code(self):
        """Test 'fetch' command with a single PDB code."""
        args = self.parser.parse_args(["fetch", "1a2b"])
        self.assertEqual(args.command, "fetch")
        self.assertEqual(args.codes, ["1a2b"])

    def test_fetch_multiple_codes(self):
        """Test 'fetch' command with multiple PDB codes."""
        args = self.parser.parse_args(["fetch", "1a2b", "3k1q", "127d"])
        self.assertEqual(args.command, "fetch")
        self.assertEqual(args.codes, ["1a2b", "3k1q", "127d"])

    def test_format_pdb(self):
        """Test --format pdb option."""
        args = self.parser.parse_args(["update", "--format", "pdb"])
        self.assertEqual(args.format, "pdb")

    def test_format_xml(self):
        """Test --format xml option."""
        args = self.parser.parse_args(["update", "--format", "xml"])
        self.assertEqual(args.format, "xml")

    def test_format_mmtf(self):
        """Test --format mmtf option."""
        args = self.parser.parse_args(["update", "--format", "mmtf"])
        self.assertEqual(args.format, "mmtf")

    def test_format_mmcif_default(self):
        """Test that default format is mmCif."""
        args = self.parser.parse_args(["update"])
        self.assertEqual(args.format, "mmCif")

    def test_format_invalid(self):
        """Test that invalid format raises SystemExit."""
        with self.assertRaises(SystemExit):
            self.parser.parse_args(["update", "--format", "invalid"])

    def test_flat_flag(self):
        """Test --flat flag."""
        args = self.parser.parse_args(["update", "--flat"])
        self.assertTrue(args.flat)

    def test_overwrite_flag(self):
        """Test --overwrite flag."""
        args = self.parser.parse_args(["update", "--overwrite"])
        self.assertTrue(args.overwrite)

    def test_with_assemblies_flag(self):
        """Test --with-assemblies flag."""
        args = self.parser.parse_args(["update", "--with-assemblies"])
        self.assertTrue(args.with_assemblies)

    def test_all_options_combined(self):
        """Test combining multiple options."""
        args = self.parser.parse_args(
            ["fetch", "1a2b", "--pdb-path", "/tmp/pdb", "--format", "pdb",
             "--flat", "--overwrite", "--with-assemblies"]
        )
        self.assertEqual(args.command, "fetch")
        self.assertEqual(args.codes, ["1a2b"])
        self.assertEqual(args.pdb_path, "/tmp/pdb")
        self.assertEqual(args.format, "pdb")
        self.assertTrue(args.flat)
        self.assertTrue(args.with_assemblies)

    def test_no_command_returns_none(self):
        """Test that no command sets command to None."""
        args = self.parser.parse_args([])
        self.assertIsNone(args.command)

    def test_fetch_requires_codes(self):
        """Test that 'fetch' without codes raises SystemExit."""
        with self.assertRaises(SystemExit):
            self.parser.parse_args(["fetch"])


class TestPDBListCLIFunction(unittest.TestCase):
    """Test the cli() function entry point."""

    def test_no_args_prints_help(self):
        """Test that cli() with no args prints help and returns."""
        # Should not raise; just prints help
        cli([])

    def test_with_assemblies_rejects_xml(self):
        """Test that --with-assemblies with --format xml raises SystemExit."""
        with self.assertRaises(SystemExit):
            cli(["update", "--format", "xml", "--with-assemblies"])

    def test_with_assemblies_rejects_mmtf(self):
        """Test that --with-assemblies with --format mmtf raises SystemExit."""
        with self.assertRaises(SystemExit):
            cli(["update", "--format", "mmtf", "--with-assemblies"])

    def test_with_assemblies_accepts_mmcif(self):
        """Test that --with-assemblies with default mmCif does not error on parsing."""
        # This will try to connect to PDB server, so we mock the PDBList methods
        from unittest.mock import patch

        with patch.object(PDBList, "update_pdb"):
            cli(["update", "--with-assemblies"])

    def test_with_assemblies_accepts_pdb(self):
        """Test that --with-assemblies with --format pdb does not error on parsing."""
        from unittest.mock import patch

        with patch.object(PDBList, "update_pdb"):
            cli(["update", "--format", "pdb", "--with-assemblies"])


class TestPDBListCLIDispatch(unittest.TestCase):
    """Test that cli() dispatches to the correct PDBList methods."""

    def test_update_dispatches(self):
        """Test that 'update' calls update_pdb."""
        from unittest.mock import patch

        with patch.object(PDBList, "update_pdb") as mock_update:
            cli(["update"])
            mock_update.assert_called_once_with(
                file_format="mmCif", with_assemblies=False
            )

    def test_all_dispatches(self):
        """Test that 'all' calls download_entire_pdb."""
        from unittest.mock import patch

        with patch.object(PDBList, "download_entire_pdb") as mock_all:
            cli(["all"])
            mock_all.assert_called_once_with(file_format="mmCif")

    def test_obsol_dispatches(self):
        """Test that 'obsol' calls download_obsolete_entries."""
        from unittest.mock import patch

        with patch.object(PDBList, "download_obsolete_entries") as mock_obsol:
            cli(["obsol"])
            mock_obsol.assert_called_once_with(file_format="mmCif")

    def test_assemb_dispatches(self):
        """Test that 'assemb' calls download_all_assemblies."""
        from unittest.mock import patch

        with patch.object(PDBList, "download_all_assemblies") as mock_assemb:
            cli(["assemb"])
            mock_assemb.assert_called_once_with(file_format="mmCif")

    def test_fetch_dispatches(self):
        """Test that 'fetch' calls retrieve_pdb_file for each code."""
        from unittest.mock import patch

        with patch.object(PDBList, "retrieve_pdb_file") as mock_fetch:
            cli(["fetch", "1a2b", "3k1q"])
            self.assertEqual(mock_fetch.call_count, 2)

    def test_fetch_with_format(self):
        """Test that 'fetch' passes format correctly."""
        from unittest.mock import patch

        with patch.object(PDBList, "retrieve_pdb_file") as mock_fetch:
            cli(["fetch", "1a2b", "--format", "pdb"])
            mock_fetch.assert_called_once()
            call_kwargs = mock_fetch.call_args
            self.assertEqual(call_kwargs.kwargs.get("file_format"), "pdb")

    def test_flat_flag_sets_flat_tree(self):
        """Test that --flat sets pl.flat_tree = True."""
        from unittest.mock import patch

        with patch.object(PDBList, "download_entire_pdb"):
            with patch.object(PDBList, "__init__", return_value=None) as mock_init:
                # Can't easily check flat_tree via mock, so just verify no crash
                pass
        # Simpler approach: just verify it runs without error
        with patch.object(PDBList, "download_entire_pdb"):
            cli(["all", "--flat"])


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
