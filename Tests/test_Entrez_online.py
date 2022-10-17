# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for online Entrez access.

This file include tests for accessing the online Entrez service and parsing the
returned results. Note that we are merely testing the access and whether the
results are parseable. Detailed tests on each Entrez service are not within the
scope of this file as they are already covered in test_Entrez.py.

"""
from Bio import Entrez
from Bio import Medline
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import doctest
import os
import sys
import unittest

import requires_internet

requires_internet.check()


# This lets us set the email address to be sent to NCBI Entrez:
Entrez.email = "biopython@biopython.org"
Entrez.api_key = "5cfd4026f9df285d6cfc723c662d74bcbe09"

URL_HEAD = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


class EntrezOnlineCase(unittest.TestCase):
    def test_no_api_key(self):
        """Test Entrez.read without API key."""
        cached = Entrez.api_key
        Entrez.api_key = None  # default
        try:
            handle = Entrez.einfo()
        finally:
            # Do not want any failure here to break other tests
            Entrez.api_key = cached
        self.assertNotIn("api_key=", handle.url)
        rec = Entrez.read(handle)
        handle.close()
        self.assertIsInstance(rec, dict)
        self.assertIn("DbList", rec)
        # arbitrary number, just to make sure that DbList has contents
        self.assertGreater(len(rec["DbList"]), 5)

    def test_read_from_url(self):
        """Test Entrez.read from URL."""
        handle = Entrez.einfo()
        rec = Entrez.read(handle)
        handle.close()
        self.assertIsInstance(rec, dict)
        self.assertIn("DbList", rec)
        # arbitrary number, just to make sure that DbList has contents
        self.assertGreater(len(rec["DbList"]), 5)

    def test_parse_from_url(self):
        """Test Entrez.parse from URL."""
        handle = Entrez.efetch(
            db="protein", id="15718680, 157427902, 119703751", retmode="xml"
        )
        recs = list(Entrez.parse(handle))
        handle.close()
        self.assertEqual(3, len(recs))
        # arbitrary number, just to make sure the parser works
        self.assertTrue(all(len(rec).keys > 5) for rec in recs)

    def test_webenv_search(self):
        """Test Entrez.search from link webenv history."""
        handle = Entrez.elink(
            db="nucleotide",
            dbfrom="protein",
            id="22347800,48526535",
            webenv=None,
            query_key=None,
            cmd="neighbor_history",
        )
        recs = Entrez.read(handle)
        handle.close()
        record = recs.pop()

        webenv = record["WebEnv"]
        query_key = record["LinkSetDbHistory"][0]["QueryKey"]
        handle = Entrez.esearch(
            db="nucleotide",
            term=None,
            retstart=0,
            retmax=10,
            webenv=webenv,
            query_key=query_key,
            usehistory="y",
        )
        search_record = Entrez.read(handle)
        handle.close()
        self.assertEqual(2, len(search_record["IdList"]))

    def test_seqio_from_url(self):
        """Test Entrez into SeqIO.read from URL."""
        handle = Entrez.efetch(
            db="nucleotide", id="186972394", rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
        self.assertIsInstance(record, SeqRecord)
        self.assertEqual("EU490707.1", record.id)
        self.assertEqual(1302, len(record))

    def test_medline_from_url(self):
        """Test Entrez into Medline.read from URL."""
        handle = Entrez.efetch(
            db="pubmed", id="19304878", rettype="medline", retmode="text"
        )
        record = Medline.read(handle)
        handle.close()
        self.assertIsInstance(record, dict)
        self.assertEqual("19304878", record["PMID"])
        self.assertEqual("10.1093/bioinformatics/btp163 [doi]", record["LID"])

    def test_efetch_taxonomy_xml(self):
        """Test Entrez using a integer id - like a taxon id."""
        handle = Entrez.efetch(db="taxonomy", id=3702, retmode="XML")
        taxon_record = Entrez.read(handle)
        self.assertTrue(1, len(taxon_record))
        self.assertIn("TaxId", taxon_record[0])
        self.assertTrue("3702", taxon_record[0]["TaxId"])
        handle.close()

    def test_elink(self):
        """Test Entrez.elink with multiple ids, both comma separated and as list.

        This is tricky because the ELink tool treats the "id" parameter
        differently than the others (see docstring of the elink() function).
        """
        params = {"db": "gene", "dbfrom": "protein"}
        ids = ["15718680", "157427902", "119703751"]

        # Pass list argument - one-to-one
        with Entrez.elink(id=ids, **params) as handle:
            result1 = Entrez.read(handle)

        self.assertEqual(len(result1), len(ids))

        id_map = {}  # Dictionary mapping each gene ID to some number of protein IDs

        for linkset in result1:
            (from_id,) = linkset["IdList"]
            (linksetdb,) = linkset["LinkSetDb"]
            to_ids = [link["Id"] for link in linksetdb["Link"]]
            # Failure here could indicate we used an invalid ID
            self.assertGreater(len(to_ids), 0)
            id_map[from_id] = to_ids

        self.assertCountEqual(id_map.keys(), ids)

        # Pass string argument - single LinkSet
        with Entrez.elink(id=",".join(ids), **params) as handle:
            result2 = Entrez.read(handle)

        (linkset,) = result2
        self.assertCountEqual(linkset["IdList"], ids)

        # Check we got the same set of IDs as in the last request
        (linksetdb,) = linkset["LinkSetDb"]
        to_ids = [link["Id"] for link in linksetdb["Link"]]
        prev_to_ids = set().union(*id_map.values())
        self.assertCountEqual(to_ids, prev_to_ids)

    def test_epost(self):
        """Test Entrez.epost with multiple ids, both comma separated and as list."""
        # TODO - check results
        handle = Entrez.epost("nuccore", id="186972394,160418")
        handle.close()
        handle = Entrez.epost("nuccore", id=["160418", "160351"])
        handle.close()

    def test_egquery(self):
        """Test Entrez.egquery.

        which searches in all Entrez databases for a single text query.
        """
        handle = Entrez.egquery(term="biopython")
        record = Entrez.read(handle)
        handle.close()

        done = False
        for row in record["eGQueryResult"]:
            if "pmc" in row["DbName"]:
                self.assertGreater(int(row["Count"]), 60)
                done = True
        self.assertTrue(done)

    def test_espell(self):
        """Test misspellings with Entrez.espell."""
        handle = Entrez.espell(term="biopythooon")
        record = Entrez.read(handle)
        handle.close()

        self.assertEqual(record["Query"], "biopythooon")
        self.assertEqual(record["CorrectedQuery"], "biopython")

    def test_ecitmatch(self):
        """Test Entrez.ecitmatch to search for a citation."""
        citation = {
            "journal_title": "proc natl acad sci u s a",
            "year": "1991",
            "volume": "88",
            "first_page": "3248",
            "author_name": "mann bj",
            "key": "citation_1",
        }
        handle = Entrez.ecitmatch(db="pubmed", bdata=[citation])
        result = handle.read()
        expected_result = (
            "proc natl acad sci u s a|1991|88|3248|mann bj|citation_1|2014248\n"
        )
        self.assertEqual(result, expected_result)
        handle.close()

    def test_efetch_ids(self):
        """Test different options to supply ids."""
        id_sets = [[15718680, 157427902], [15718680]]
        id_vals = [
            [
                [15718680, 157427902],
                (15718680, 157427902),
                {15718680, 157427902},
                ["15718680", "157427902"],
                ("15718680", "157427902"),
                {15718680, "157427902"},
                "15718680, 157427902",
            ],
            [[15718680], (15718680), {15718680}, 15718680, "15718680", "15718680,"],
        ]

        for ids, vals in zip(id_sets, id_vals):
            for _id in vals:
                with Entrez.efetch(db="protein", id=_id, retmode="xml") as handle:
                    recs = list(Entrez.parse(handle))

                # Extract the numerical IDs of returned records
                rec_ids = [
                    int(seqid[3:])
                    for rec in recs
                    for seqid in rec["GBSeq_other-seqids"]
                    if seqid.startswith("gi|")
                ]
                self.assertCountEqual(rec_ids, ids)

    def test_efetch_gds_utf8(self):
        """Test correct handling of encodings in Entrez.efetch."""
        # See issue #1402 in case any encoding issues occur
        handle = Entrez.efetch(db="gds", id="200079209")
        text = handle.read()
        # Use of Unicode double quotation marks U+201C and U+201D
        expected_phrase = "“field of injury”"
        self.assertEqual(text[342:359], expected_phrase)
        handle.close()

    def test_fetch_xml_schemas(self):
        handle = Entrez.efetch("protein", id="783730874", rettype="ipg", retmode="xml")
        record = Entrez.read(handle, validate=False)
        handle.close()
        self.assertEqual(len(record), 1)
        self.assertIn("IPGReport", record)
        self.assertIn("Product", record["IPGReport"])
        self.assertIn("Statistics", record["IPGReport"])
        self.assertIn("ProteinList", record["IPGReport"])


if __name__ == "__main__":
    # When running test_Entrez.py directly, will also include the
    # Bio.Entrez doctests.
    # TODO: Include the doctests via run_tests.py when online.
    unittest_suite = unittest.TestLoader().loadTestsFromName("test_Entrez_online")
    doctest_suite = doctest.DocTestSuite(Entrez)
    suite = unittest.TestSuite((unittest_suite, doctest_suite))
    runner = unittest.TextTestRunner(sys.stdout, verbosity=2)
    runner.run(suite)
