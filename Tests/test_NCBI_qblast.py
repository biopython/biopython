# Copyright 2008-2016 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Testing online code for fetching NCBI qblast.

Uses Bio.Blast.NCBIWWW.qblast() to run some online blast queries, get XML
blast results back, and then checks Bio.Blast.NCBIXML.parse() can read them.

Goals:
    Make sure that all retrieval is working as expected.
    Make sure we can parse the latest XML format being used by the NCBI.
"""
from __future__ import print_function
import unittest

from Bio._py3k import HTTPError
from Bio._py3k import StringIO

from Bio import MissingExternalDependencyError

# We want to test these:
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

import requires_internet
requires_internet.check()

#####################################################################

# List of qblast requests stored as a tuple of parameters:
# - program
# - database
# - query identifier or sequence
# - expectation value threshold
# - Entrez filter string (or None)
# - list of hit identifiers expected to be found (or None if expect 0)

print("Checking Bio.Blast.NCBIWWW.qblast() with various queries")


class TestQblast(unittest.TestCase):

    def test_blastp_nr_actin(self):
        # Simple protein blast filtered for rat only, using protein
        # GI:160837788 aka NP_075631.2
        # the actin related protein 2/3 complex, subunit 1B [Mus musculus]
        self.run_qblast("blastp", "nr", "NP_075631.2", 0.001,
                        "rat [ORGN]", ['9506405', '13592137', '37589612', '149064087', '56912225'])

    def test_pcr_primers(self):
        # This next example finds PCR primer matches in Chimpanzees, e.g. BRCA1:
        self.run_qblast("blastn", "nr", "GTACCTTGATTTCGTATTC" + ("N" * 30) + "GACTCTACTACCTTTACCC",
                        10, "pan [ORGN]", ["XM_009432096.2", "XM_009432102.2", "XM_009432101.2",
                                           "XM_016930487.1", "XM_009432104.2", "XM_009432099.2",
                                           "XR_001710553.1", "XM_016930485.1", "XM_009432089.2",
                                           "XM_016930484.1"])

    def test_orchid_est(self):
        # Try an orchid EST (nucleotide) sequence against NR using BLASTX
        self.run_qblast("blastx", "nr",
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
                        0.0000001, None, ["21554275", "18409071", "296087288", "566183510"])

    def run_qblast(self, program, database, query, e_value, entrez_filter, expected_hits):
        try:
            if program == "blastn":
                # Check the megablast parameter is accepted
                handle = NCBIWWW.qblast(program, database, query,
                                        alignments=10, descriptions=10,
                                        hitlist_size=10,
                                        entrez_query=entrez_filter,
                                        expect=e_value, megablast="FALSE")
            else:
                handle = NCBIWWW.qblast(program, database, query,
                                        alignments=10, descriptions=10,
                                        hitlist_size=10,
                                        entrez_query=entrez_filter,
                                        expect=e_value)
        except HTTPError:
            # e.g. a proxy error
            raise MissingExternalDependencyError("internet connection failed")
        record = NCBIXML.read(handle)

        if record.query == "No definition line":
            # We used a sequence as the query
            self.assertEqual(len(query), record.query_letters)
        elif query.startswith(">"):
            # We used a FASTA record as the query
            expected = query[1:].split("\n", 1)[0]
            self.assertEqual(expected, record.query)
        elif record.query_id.startswith("Query_") and len(query) == record.query_letters:
            # We used a sequence as the entry and it was given a placeholder name
            pass
        else:
            # We used an identifier as the query
            self.assertIn(query, record.query_id.split("|"),
                          "Expected %r within query_id %r" % (query, record.query_id))

        # Check the recorded input parameters agree with those requested
        self.assertEqual(float(record.expect), e_value)
        self.assertEqual(record.application.lower(), program)
        self.assertTrue(len(record.alignments) <= 10)
        self.assertTrue(len(record.descriptions) <= 10)

        # Check the expected result(s) are found in the alignments
        if expected_hits is None:
            self.assertEqual(len(record.alignments), 0)  # Expected no alignments!
        else:
            self.assertTrue(len(record.alignments) > 0)  # Expected some alignments!
            found_result = False
            for expected_hit in expected_hits:
                for alignment in record.alignments:
                    if expected_hit in alignment.hit_id.split("|"):
                        found_result = True
                        break
            if len(expected_hits) == 1:
                print("Update this test to have some redundancy...")
                for alignment in record.alignments:
                    print(alignment.hit_id)
            self.assertTrue(found_result,
                            "Missing all expected hits (%s), instead have: %s"
                            % (", ".join(expected_hits),
                               ", ".join(a.hit_id for a in record.alignments)))

        # Check the expected result(s) are found in the descriptions
        if expected_hits is None:
            self.assertEqual(len(record.descriptions), 0)  # Expected no descriptions!
        else:
            self.assertTrue(len(record.descriptions) > 0)  # Expected some descriptions!
            found_result = False
            for expected_hit in expected_hits:
                for descr in record.descriptions:
                    if expected_hit == descr.accession \
                            or expected_hit in descr.title.split(None, 1)[0].split("|"):
                        found_result = True
                        break
            assert found_result, "Missing all of %s in descriptions" % expected_hit
            self.assertTrue(found_result)

    def test_parse_qblast_ref_page(self):
        with open("Blast/html_msgid_29_blastx_001.html", "r") as f:
            handle = StringIO(f.read())
        self.assertRaises(ValueError, NCBIWWW._parse_qblast_ref_page, handle)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
