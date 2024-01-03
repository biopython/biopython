# Copyright 2008-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Testing online code for fetching NCBI qblast.

Uses Bio.Blast.NCBIWWW.qblast() to run some online blast queries, get XML
blast results back, and then checks Bio.Blast.NCBIXML.parse() can read them.

Goals:
    - Make sure that all retrieval is working as expected.
    - Make sure we can parse the latest XML format being used by the NCBI.

If an internet connection is available and run_tests.py has not given the argument
'--offline', these tests will run online (and can take a long time). Otherwise, the
module runs offline using a mocked version of urllib.request.urlopen, which returns
file contents instead of internet responses.

IMPORTANT:
If you add new tests (or change existing tests) you must provide new 'mock' xml
files which fulfill the requirements of the respective tests. These need to be
added to the 'response_list' within the 'mock_response' function. Note: The
tests are run in alphabetical order, so you must place your mock file at the
correct position.

"""
import unittest
from unittest import mock

from urllib.error import HTTPError
from io import BytesIO

from Bio import MissingExternalDependencyError

# We want to test these:
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

import requires_internet


NCBIWWW.email = "biopython@biopython.org"


try:
    requires_internet.check()
    # this will set 'require_internet.check.available' (if not already done before)
    # it is False if we don't have internet or if we run the tests with --offline
except MissingExternalDependencyError:
    # need to catch this error, because we want to run the tests even offline
    pass  # 'requires_internet_check.available' is now False

if not requires_internet.check.available:
    # if the tests are running with --offline or if no internet is available,
    # we mock urlopen, so that it returns files from the blast folder

    def mock_response():
        """Mimic an NCBI qblast response."""
        # Each use of NCBIWWW.qblast makes two urlopen calls with different responses:
        # a. the 'wait' page, and b. the result.
        wait = ["Blast/mock_wait.html"]  # This mimics the 'wait' page

        # These mimic the results. Add new mock files here, if you add more tests.
        # Note: The test are run in alphabetical order, so place new files at the
        # correct position.
        response_list = [
            "Blast/mock_actin.xml",  # result for test_blastp_nr_actin
            "Blast/mock_disco.xml",  # result for test_discomegalast
            "Blast/mock_orchid.xml",  # result for test_orchid_est
            "Blast/mock_pcr.xml",  # result for test_pcr_primers
            "Blast/mock_short_empty.xml",  # result for test_short_query # 1
            "Blast/mock_short_result.xml",  # result for test_short_query # 2
            "Blast/mock_short_result.xml",  # result for test_short_query # 3
        ]

        # Generate a list of responses with the structure wait|result|wait|result...
        responses = (
            BytesIO(open(a, "rb").read())
            for b in zip(len(response_list) * wait, response_list)
            for a in b
        )
        return responses

    NCBIWWW.time.sleep = mock.Mock()  # we don't want to wait 20 sec for each test...
    NCBIWWW.urlopen = mock.Mock(side_effect=mock_response())


#####################################################################

# List of qblast requests stored as a tuple of parameters:
# - program
# - database
# - query identifier or sequence
# - expectation value threshold
# - Entrez filter string (or None)
# - list of hit identifiers expected to be found (or None if expect 0)


class TestQblast(unittest.TestCase):
    def test_blastp_nr_actin(self):
        # Simple protein blast filtered for rat only, using protein
        # GI:160837788 aka NP_075631.2
        # the actin related protein 2/3 complex, subunit 1B [Mus musculus]
        self.run_qblast(
            "blastp",
            "nr",
            "NP_075631.2",
            0.001,
            "rat [ORGN]",
            {"megablast": "FALSE"},
            [
                "NP_112408.1",
                "AAH59131.1",
                "EDM14357.1",
                "NP_001008766.1",
                "NP_001102411.1",
                "EDL80109.1",
                "EDL80106.1",
                "NP_001100434.1",
                "AAI67084.1",
            ],
        )

    def test_pcr_primers(self):
        # This next example finds PCR primer matches in Chimpanzees, e.g. BRCA1:
        self.run_qblast(
            "blastn",
            "nr",
            "GTACCTTGATTTCGTATTC" + ("N" * 30) + "GACTCTACTACCTTTACCC",
            10,
            "pan [ORGN]",
            {"megablast": "FALSE"},
            [
                "XM_034941187.1",
                "XM_034941186.1",
                "XM_034941185.1",
                "XM_034941184.1",
                "XM_034941183.1",
                "XM_034941182.1",
                "XM_034941180.1",
                "XM_034941179.1",
                "XM_034941178.1",
                "XM_034941177.1",
            ],
        )

    def test_orchid_est(self):
        # Try an orchid EST (nucleotide) sequence against NR using BLASTX
        self.run_qblast(
            "blastx",
            "nr",
            """>gi|116660609|gb|EG558220.1|EG558220 CR02019H04 Leaf CR02 cDNA library Catharanthus roseus cDNA clone CR02019H04 5', mRNA sequence
               CTCCATTCCCTCTCTATTTTCAGTCTAATCAAATTAGAGCTTAAAAGAATGAGATTTTTAACAAATAAAA
               AAACATAGGGGAGATTTCATAAAAGTTATATTAGTGATTTGAAGAATATTTTAGTCTATTTTTTTTTTTT
               TCTTTTTTTGATGAAGAAAGGGTATATAAAATCAAGAATCTGGGGTGTTTGTGTTGACTTGGGTCGGGTG
               TGTATAATTCTTGATTTTTTCAGGTAGTTGAAAAGGTAGGGAGAAAAGTGGAGAAGCCTAAGCTGATATT
               GAAATTCATATGGATGGAAAAGAACATTGGTTTAGGATTGGATCAAAAAATAGGTGGACATGGAACTGTA
               CCACTACGTCCTTACTATTTTTGGCCGAGGAAAGATGCTTGGGAAGAACTTAAAACAGTTTTAGAAAGCA
               AGCCATGGATTTCTCAGAAGAAAATGATTATACTTCTTAATCAGGCAACTGATATTATCAATTTATGGCA
               GCAGAGTGGTGGCTCCTTGTCCCAGCAGCAGTAATTACTTTTTTTTCTCTTTTTGTTTCCAAATTAAGAA
               ACATTAGTATCATATGGCTATTTGCTCAATTGCAGATTTCTTTCTTTTGTGAATG""",
            0.0000001,
            None,
            {"megablast": "FALSE"},
            [
                "XP_021665344.1",
                "XP_021615158.1",
                "XP_017223689.1",
                "OMP06800.1",
                "XP_021634873.1",
                "XP_021299673.1",
                "XP_002311451.2",
                "XP_021976565.1",
                "OMO90244.1",
            ],
        )

    def test_discomegablast(self):
        self.run_qblast(
            "blastn",
            "nr",
            """>some sequence
               ATGAAGATCTTCCAGATCCAGTGCAGCAGCTTCAAGGAGAGCAGGTGGCAGAAGAGCAAGTGCGACAACT
               GCCTGAAGTTCCACATCGACATCAACAACAACAGCAAGACCAGCAACACCGACACCGACTTCGACGCCAA
               CACCAACATCAACAGCAACATCAACAGCAACATCAACAGCAACATCAACATCAACAACAGCGGCAACAAC
               AACAAGAACAGCAACAACATCGAGATCACCGAGAACATCGACAACAAGGCCAAGATCATCAACAAGCACA
               TCAAGACCATCACCAACAGCAAGCCCATCCCCATCCCCATCCCCACCCCCACCCCCATCAGCATCAAGGA
               GAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGATG
               AAGAGCACCATCAACCTGAGGAGCGAGGACACCACCAGCAACAAGAGCACCATCGTGTTCACCGAGTGCC
               TGGAGTACAAGGGCCACCAGTGGAGGCCCAACATCTGCGTGACCTGCTTCAGCCCCAAGAACAAGCACAA
               GAACGTGCTGCCCGAGACCAGCACCCCCCTGATCAGCCAGAGCAGCCAGACCAGCACCATCACCCCCAGC
               AGCAGCAGCACCAGCACCAGCACCAGCAGCATCAGCACCCACAAGACCGCCAACAACAAGACCGTGATCA
               CCTACATCAGCAGCACCACCACCACCACCACCACCAGCAGCAGCAGCAGCAGCCCCCCCAGCAGCAGCAT
               CGCCGGCATCACCAACCCCACCAGCAGGAGCAGCAGCCCCATCCTGAAGAGCGTGCCCCCCAGCGCCTAC
               AGCAACGTGGTGATCCCCATCAACAACATCAACAACAGCAACAGCAACAGCAGCAGCGGCGGCGGCAACA
               ACAACAACAAGAGCATCAGCACCCCCAGCAGCCCCATCATCAGCAGGCCCATCACCAACAAGATCAACAA
               CAACAACAACAACAACCAGCCCCAGCTGCACTACAACCAGCCCCAGAGCAGCAGCGTGAGCACCACCAGC
               AGCCCCATCATCAGGCCCGTGCTGAGGAGGCAGTTCCAGAGCTTCCCCAGCAACCCCAAGATCAGCAAGG
               CCATCCTGGAGCAGTGCAACATCATCAACAACAACAGCAACAGCAACAACAGCAACAACAAGGACCCCGT
               GATCCTGTGCAAGTACACCATCGAGAGCCAGCCCAAGAGCAACATCAGCGTGCTGAAGCCCACCCTGGTG
               GAGTTCATCAACCAGCCCGACAGCAAGGACGACGAGAGCAGCGTGAAGAGCCCCCCCCTGCCCGTGGAGA
               GCCAGCCCATCTTCAACAGCAAGCAGAGCGCCACCATGGACGGCATCACCACCCACAAGAGCGTGAGCAT
               CACCATCAGCACCAGCACCAGCCCCAGCAGCACCACCACCACCACCAGCACCACCACCAGCATCATCGCC
               GAGGAGCCCAGCAGCCCCATCCTGCCCACCGCCAGCCCCAGCAGCAGCAGCAGCAGCATCATCACCACCG
               CCACCGCCAGCACCATCCCCATGAGCCCCAGCCTGCCCAGCATCCCCTTCCACGAGTTCGAGACCATGGA
               GAGCAGCACCACCACCACCCTGCTGAGCGAGAACAACGGCGGCGGCGGCGGCAGCAGCTGCAACGACAAC
               AGCAGGAGGAACAGCCTGAACATCCTGCCCCTGAGGCTGAAGAGCTTCAGCTTCAGCGCCCCCCAGAGCG
               ACAGCATGATCGAGCAGCCCGAGGACGACCCCTTCTTCGACTTCGAGGACCTGAGCGACGACGACGACAG
               CAACGACAACGACGACGAGGAGCTGAAGGAGATCAACGGCGAGAAGATCATCCAGCAGAACGACCTGACC
               CCCACCACCACCATCACCAGCACCACCACCATCCTGCAGAGCCCCACCCTGGAGAAGACCCTGAGCACCA
               CCACCACCACCACCATCCCCAGCCCCAGCACCAACAGCAGGAGCATCTGCAACACCCTGATGGACAGCAC
               CGACAGCATCAACAACACCAACACCAACACCAACACCAACACCAACACCAACACCAACACCAACACCAAC
               ACCAACACCAACACCAACACCAACGCCAACATCAACAACAAGGTGAGCACCACCACCACCACCACCACCA
               CCAAGAGGAGGAGCCTGAAGATGGACCAGTTCAAGGAGAAGGAGGACGAGTGGGACCAGGGCGTGGACCT
               GACCAGCTTCCTGAAGAGGAAGCCCACCCTGCAGAGGGACTTCAGCTACTGCAACAACAAGGTGATGGAG
               ATCAGCAGCGTGAAGGAGGAGGCCAAGAGGCTGCACGGCGGCACCGGCTACATCCACCAGTTCGCCTTCG
               AGGCCTTCAAGGACATCCTGGAGGCCAAGCAGACCCAGATCAACAGGGCCTTCTGCAGCCAGAAGATCGA
               CGCCCCCGACTGCGAGATGCTGATCAACGAGATCAACACCGCCAAGAAGCTGCTGGAGGACCTGCTGGAG
               CTGAACAGCAACAGCAGCGGCAGCGGCAACAACAGCAACGACAACAGCGGCAGCAGCAGCCCCAGCAGCA
               GCAAGACCAACACCCTGAACCAGCAGAGCATCTGCATCAAGAGCGAGATCCAACGATACGTTGAAATTCG
               CTTGTGTGCCACTGGTAAATCCACCCCCCCTAAGCCTCTAATAGGGAGACCTTAG""",
            0.0000001,
            None,
            {"template_type": 0, "template_length": 18, "megablast": "on"},  # noqa 231
            ["XM_635681.1", "XM_008496783.1"],
        )

    def run_qblast(
        self,
        program,
        database,
        query,
        e_value,
        entrez_filter,
        additional_args,
        expected_hits,
    ):
        """Do qblast searches with given parameters and analyze results."""
        try:
            if program == "blastn":
                # Check the megablast parameter is accepted
                handle = NCBIWWW.qblast(
                    program,
                    database,
                    query,
                    alignments=10,
                    descriptions=10,
                    hitlist_size=10,
                    entrez_query=entrez_filter,
                    expect=e_value,
                    **additional_args,
                )
            else:
                handle = NCBIWWW.qblast(
                    program,
                    database,
                    query,
                    alignments=10,
                    descriptions=10,
                    hitlist_size=10,
                    entrez_query=entrez_filter,
                    expect=e_value,
                    **additional_args,
                )
        except HTTPError:
            # e.g. a proxy error
            raise MissingExternalDependencyError("internet connection failed") from None

        record = NCBIXML.read(handle)

        if record.query == "No definition line":
            # We used a sequence as the query
            self.assertEqual(len(query), record.query_letters)
        elif query.startswith(">"):
            # We used a FASTA record as the query
            expected = query[1:].split("\n", 1)[0]
            self.assertEqual(expected, record.query)
        elif (
            record.query_id.startswith("Query_") and len(query) == record.query_letters
        ):
            # We used a sequence as the entry and it was given a placeholder name
            pass
        else:
            # We used an identifier as the query
            self.assertIn(
                query,
                record.query_id.split("|"),
                f"Expected {query!r} within query_id {record.query_id!r}",
            )

        # Check the recorded input parameters agree with those requested
        self.assertEqual(float(record.expect), e_value)
        self.assertEqual(record.application.lower(), program)
        self.assertLessEqual(len(record.alignments), 10)
        self.assertLessEqual(len(record.descriptions), 10)

        # Check the expected result(s) are found in the alignments
        if expected_hits is None:
            self.assertEqual(len(record.alignments), 0)  # Expected no alignments!
        else:
            self.assertGreater(len(record.alignments), 0)  # Expected some alignments!
            found_result = False
            for expected_hit in expected_hits:
                for alignment in record.alignments:
                    if expected_hit in alignment.hit_id.split("|"):
                        found_result = True
                        break
            self.assertTrue(
                found_result,
                "Missing all expected hits (%s), instead have: %s"
                % (
                    ", ".join(expected_hits),
                    ", ".join(a.hit_id for a in record.alignments),
                ),
            )

        # Check the expected result(s) are found in the descriptions
        if expected_hits is None:
            # Expected no descriptions!
            self.assertEqual(len(record.descriptions), 0)
        else:
            # Expected some descriptions!
            self.assertGreater(len(record.descriptions), 0)
            found_result = False
            for expected_hit in expected_hits:
                for descr in record.descriptions:
                    if (
                        expected_hit == descr.accession
                        or expected_hit in descr.title.split(None, 1)[0].split("|")
                    ):
                        found_result = True
                        break
            msg = f"Missing all of {expected_hit} in descriptions"
            self.assertTrue(found_result, msg=msg)

    def test_parse_qblast_ref_page(self):
        with open("Blast/html_msgid_29_blastx_001.html", "rb") as f:
            handle = BytesIO(f.read())
        self.assertRaises(ValueError, NCBIWWW._parse_qblast_ref_page, handle)

    def test_short_query(self):
        """Test SHORT_QUERY_ADJUST parameter."""
        # Should give no hits:
        my_search = NCBIWWW.qblast("blastp", "nr", "ICWENRMP", hitlist_size=5)
        my_hits = NCBIXML.read(my_search)
        my_search.close()
        self.assertEqual(len(my_hits.alignments), 0)

        # Should give hits:
        my_search = NCBIWWW.qblast(
            "blastp", "nr", "ICWENRMP", hitlist_size=5, short_query=True
        )
        my_hits = NCBIXML.read(my_search)
        my_search.close()
        self.assertNotEqual(len(my_hits.alignments), 0)

        # Disabled for 1.81 release
        # Query does not seem to be returning anything so
        # maybe the test needs to be updated?
        #
        # with warnings.catch_warnings(record=True) as w:
        #     # Cause all warnings to always be triggered.
        #     warnings.simplefilter("always")
        #     # Trigger a warning.
        #     my_search = NCBIWWW.qblast(
        #         "blastn", "nt", "ATGTCAACTTCAGAA", hitlist_size=5, short_query=True
        #     )
        #     # Verify some things
        #     self.assertEqual(w[0].category, BiopythonWarning)
        #     self.assertIn("blastn", str(w[0].message))
        #     my_hits = NCBIXML.read(my_search)
        #     my_search.close()
        #     self.assertNotEqual(len(my_hits.alignments), 0)

    def test_error_conditions(self):
        """Test if exceptions were properly handled."""
        self.assertRaises(
            ValueError,
            NCBIWWW.qblast,
            "megablast",
            "nt",
            "ATGCGTACGCAGCTAAAGTAAACCTATCGCGTCTCCT",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
