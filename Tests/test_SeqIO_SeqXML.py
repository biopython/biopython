# Copyright 2010 by Thomas Schmitt.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for SeqIO SeqXML module."""
import sys
import unittest

from io import BytesIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestSimpleRead(unittest.TestCase):
    def test_check_SeqIO(self):
        """Files readable using parser via SeqIO."""
        records = list(SeqIO.parse("SeqXML/dna_example.xml", "seqxml"))
        self.assertEqual(len(records), 4)
        records = list(SeqIO.parse("SeqXML/rna_example.xml", "seqxml"))
        self.assertEqual(len(records), 5)
        records = list(SeqIO.parse("SeqXML/protein_example.xml", "seqxml"))
        self.assertEqual(len(records), 5)
        records = list(SeqIO.parse("SeqXML/global_species_example.xml", "seqxml"))
        self.assertEqual(len(records), 2)


class TestDetailedRead(unittest.TestCase):

    records = {}

    def setUp(self):
        self.records["dna"] = list(SeqIO.parse("SeqXML/dna_example.xml", "seqxml"))
        self.records["rna"] = list(SeqIO.parse("SeqXML/rna_example.xml", "seqxml"))
        self.records["protein"] = list(
            SeqIO.parse("SeqXML/protein_example.xml", "seqxml")
        )
        self.records["globalSpecies"] = list(
            SeqIO.parse("SeqXML/global_species_example.xml", "seqxml")
        )

    def test_special_characters_desc(self):
        """Read special XML characters in description."""
        self.assertEqual(
            self.records["dna"][2].description,
            'some special characters in the description\n<tag> "quoted string"',
        )

    def test_unicode_characters_desc(self):
        """Test special unicode characters in the description."""
        self.assertEqual(self.records["rna"][2].description, "åÅüöÖßøä¢£$€香肠")

    def test_full_characters_set_read(self):
        """Read full characters set for each type."""
        self.assertEqual(self.records["dna"][1].seq, "ACGTMRWSYKVHDBXN.-")
        self.assertEqual(self.records["rna"][1].seq, "ACGUMRWSYKVHDBXN.-")
        self.assertEqual(
            self.records["protein"][1].seq, "ABCDEFGHIJKLMNOPQRSTUVWXYZ.-*"
        )

    def test_duplicated_property(self):
        """Read property with multiple values."""
        self.assertEqual(
            self.records["protein"][2].annotations["test"], ["1", "2", "3"]
        )

    def test_duplicated_dbxref(self):
        """Read multiple cross references to a single source."""
        self.assertEqual(
            self.records["protein"][2].dbxrefs, ["someDB:G001", "someDB:G002"]
        )

    def test_read_minimal_required(self):
        """Check minimal record."""
        minimalRecord = SeqRecord(id="test", seq=Seq("abc"))
        minimalRecord.annotations["source"] = "Ensembl"
        minimalRecord.annotations["molecule_type"] = "DNA"

        self.assertEqual(self.records["rna"][3].name, minimalRecord.name)
        self.assertEqual(self.records["dna"][3].annotations, minimalRecord.annotations)
        self.assertEqual(self.records["rna"][3].dbxrefs, minimalRecord.dbxrefs)
        self.assertEqual(
            self.records["protein"][3].description, minimalRecord.description
        )

    def test_local_species(self):
        """Check local species."""
        self.assertEqual(self.records["rna"][1].annotations["organism"], "Mus musculus")
        self.assertEqual(self.records["rna"][1].annotations["ncbi_taxid"], "10090")

        self.assertEqual(
            self.records["rna"][0].annotations["organism"], "Gallus gallus"
        )
        self.assertEqual(self.records["rna"][0].annotations["ncbi_taxid"], "9031")

    def test_global_species(self):
        """Check global species."""
        self.assertEqual(
            self.records["globalSpecies"][0].annotations["organism"], "Mus musculus"
        )
        self.assertEqual(
            self.records["globalSpecies"][0].annotations["ncbi_taxid"], "10090"
        )

        self.assertEqual(
            self.records["globalSpecies"][1].annotations["organism"], "Homo sapiens"
        )
        self.assertEqual(
            self.records["globalSpecies"][1].annotations["ncbi_taxid"], "9606"
        )

    def test_local_source_definition(self):
        """Check local source."""
        self.assertEqual(self.records["protein"][4].annotations["source"], "Uniprot")

    def test_empty_description(self):
        """Check empty description."""
        self.assertEqual(
            self.records["rna"][4].description,
            SeqRecord(id="", seq=Seq("")).description,
        )


class TestReadHeader(unittest.TestCase):
    def test_check_dna_header(self):
        """Check if the header information is parsed."""
        records = SeqIO.parse("SeqXML/dna_example.xml", "seqxml")
        self.assertEqual(records.source, "Ensembl")
        self.assertEqual(records.sourceVersion, "56")
        self.assertEqual(records.seqXMLversion, "0.4")

    def test_check_rna_header(self):
        """Check if the header information is parsed."""
        records = SeqIO.parse("SeqXML/rna_example.xml", "seqxml")
        self.assertEqual(records.source, "Ensembl")
        self.assertEqual(records.sourceVersion, "56")
        self.assertEqual(records.seqXMLversion, "0.3")

    def test_check_protein_header(self):
        """Check if the header information is parsed."""
        records = SeqIO.parse("SeqXML/protein_example.xml", "seqxml")
        self.assertEqual(records.source, "Ensembl")
        self.assertEqual(records.sourceVersion, "56")
        self.assertEqual(records.seqXMLversion, "0.4")

    def test_check_global_species_example_header(self):
        """Check if the header information is parsed."""
        records = SeqIO.parse("SeqXML/global_species_example.xml", "seqxml")
        self.assertEqual(records.speciesName, "Mus musculus")
        self.assertEqual(records.ncbiTaxID, "10090")
        self.assertEqual(records.source, "Ensembl")
        self.assertEqual(records.sourceVersion, "56")
        self.assertEqual(records.seqXMLversion, "0.4")


class TestReadAndWrite(unittest.TestCase):
    def test_read_write_rna(self):
        """Read and write RNA."""
        read1_records = list(SeqIO.parse("SeqXML/rna_example.xml", "seqxml"))
        self._write_parse_and_compare(read1_records)

    def test_read_write_dna(self):
        """Read and write DNA."""
        read1_records = list(SeqIO.parse("SeqXML/dna_example.xml", "seqxml"))
        self._write_parse_and_compare(read1_records)

    def test_read_write_protein(self):
        """Read and write protein."""
        read1_records = list(SeqIO.parse("SeqXML/protein_example.xml", "seqxml"))
        self._write_parse_and_compare(read1_records)

    def test_read_write_globalSpecies(self):
        """Read and write global species."""
        read1_records = list(SeqIO.parse("SeqXML/global_species_example.xml", "seqxml"))
        self._write_parse_and_compare(read1_records)

    def _write_parse_and_compare(self, read1_records):

        handle = BytesIO()

        SeqIO.write(read1_records, handle, "seqxml")

        handle.seek(0)
        read2_records = list(SeqIO.parse(handle, "seqxml"))

        self.assertEqual(len(read1_records), len(read2_records))

        for record1, record2 in zip(read1_records, read2_records):
            self.assertEqual(record1.id, record2.id)
            self.assertEqual(record1.name, record2.name)
            self.assertEqual(record1.description, record2.description)
            self.assertEqual(record1.seq, record2.seq)
            self.assertEqual(record1.dbxrefs, record2.dbxrefs)
            self.assertEqual(record1.annotations, record2.annotations)

    def test_write_species(self):
        """Test writing species from annotation tags."""
        record = SeqIO.read("SwissProt/sp016", "swiss")
        self.assertEqual(record.annotations["organism"], "Homo sapiens (Human)")
        self.assertEqual(record.annotations["ncbi_taxid"], ["9606"])
        handle = BytesIO()
        SeqIO.write(record, handle, "seqxml")
        handle.seek(0)
        output = handle.getvalue()
        text = output.decode("UTF-8")
        self.assertIn("Homo sapiens (Human)", text)
        self.assertIn("9606", text)
        self.assertIn(
            '<species name="Homo sapiens (Human)" ncbiTaxID="9606"></species>',
            text,
            msg=f"Missing expected <species> tag: {text!r}",
        )


class TestReadCorruptFiles(unittest.TestCase):
    def test_for_errors(self):
        """Handling of corrupt files."""
        # SeqIO.parse reads the file in blocks until it finds the seqXML
        # element with global information such as the source and sourceVersion.
        # Since one block is likely large enough to cover the first few
        # entries in the file, the ValueError may be raised after we call
        # SeqIO.parse, before we start iterating over the file.
        def f(path):
            records = SeqIO.parse(path, "seqxml")
            for record in records:
                pass

        self.assertRaises(ValueError, f, "SeqXML/corrupt_example1.xml")
        self.assertRaises(ValueError, f, "SeqXML/corrupt_example2.xml")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
