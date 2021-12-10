# Copyright 2019-2020 by Sergio Valqui. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.UniProt.GOA Module.

GOA files can be found here ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/
"""

import unittest
import tempfile
import os
from Bio.UniProt import GOA


class GoaTests(unittest.TestCase):
    """Test for UniProt GOA Files."""

    def test_gaf_iterator(self):
        """Test GOA GAF file iterator."""
        # Test GAF 2.0
        recs = []
        with open("UniProt/goa_yeast.gaf") as handle:
            for rec in GOA.gafiterator(handle):
                recs.append(rec)

        # Check number of records
        self.assertEqual(len(recs), 587)
        # Check keys are same as predefined fields
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GAF20FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "A0A023PXA5")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "YAL019W-A")
        self.assertEqual(recs[0]["Qualifier"], [""])
        self.assertEqual(recs[0]["GO_ID"], "GO:0003674")
        self.assertEqual(recs[0]["DB:Reference"], ["GO_REF:0000015"])
        self.assertEqual(recs[0]["Evidence"], "ND")
        self.assertEqual(recs[0]["With"], [""])

        # Test GAF 2.1, it has the same fields as GAF 2.0
        recs = []
        with open("UniProt/gene_association.goa_yeast.1.gaf") as handle:
            for rec in GOA.gafiterator(handle):
                recs.append(rec)

        # Check number of records
        self.assertEqual(len(recs), 300)
        # Check keys are same as predefined fields
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GAF20FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "P17536")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "TPM1")
        self.assertEqual(recs[0]["Qualifier"], [""])
        self.assertEqual(recs[0]["GO_ID"], "GO:0000001")
        self.assertEqual(recs[0]["DB:Reference"], ["PMID:10652251"])
        self.assertEqual(recs[0]["Evidence"], "TAS")
        self.assertEqual(recs[0]["With"], [""])

    def test_gpa_iterator(self):
        """Test GOA GPA file iterator."""
        recs = []
        with open("UniProt/goa_yeast.gpa.59.gpa") as handle:
            for rec in GOA.gpa_iterator(handle):
                recs.append(rec)
        self.assertEqual(len(recs), 300)
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GPA11FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "A0A023PXA5")
        self.assertEqual(recs[0]["Qualifier"], ["enables"])
        self.assertEqual(recs[0]["GO_ID"], "GO:0003674")
        self.assertEqual(recs[0]["DB:Reference"], ["GO_REF:0000015"])
        self.assertEqual(recs[0]["ECO_Evidence_code"], "ECO:0000307")
        self.assertEqual(recs[0]["With"], [""])
        self.assertEqual(recs[0]["Interacting_taxon_ID"], "")
        self.assertEqual(recs[0]["Date"], "20030730")
        self.assertEqual(recs[0]["Assigned_by"], "SGD")
        self.assertEqual(recs[0]["Annotation Extension"], [""])
        self.assertEqual(recs[0]["Annotation_Properties"], "go_evidence=ND")

    def test_gpi_iterator(self):
        """Test GOA GPI file iterator, gpi-version: 1.1."""
        recs = []
        with open("UniProt/gp_information.goa_yeast.28.gpi") as handle:
            for rec in GOA.gpi_iterator(handle):
                recs.append(rec)
        self.assertEqual(len(recs), 300)
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GPI11FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB_Object_ID"], "A2P2R3")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "YMR084W")
        self.assertEqual(
            recs[0]["DB_Object_Name"],
            [
                "Putative glutamine--fructose"
                "-6-phosphate aminotransferase"
                " [isomerizing]"
            ],
        )
        self.assertEqual(recs[0]["DB_Object_Synonym"], ["YM084_YEAST", "YMR084W"])
        self.assertEqual(recs[0]["DB_Object_Type"], "protein")
        self.assertEqual(recs[0]["Taxon"], "taxon:559292")
        self.assertEqual(recs[0]["Parent_Object_ID"], "")
        self.assertEqual(recs[0]["DB_Xref"], [""])
        self.assertEqual(recs[0]["Gene_Product_Properties"], ["db_subset=Swiss-Prot"])

    def test_gpi_iterator_one_two(self):
        """Test GOA GPI file iterator, gpi-version: 1.2."""
        recs = []
        with open("UniProt/goa_human_sample.gpi") as handle:
            for rec in GOA.gpi_iterator(handle):
                recs.append(rec)
        self.assertEqual(len(recs), 9)
        self.assertEqual(sorted(recs[0].keys()), sorted(GOA.GPI12FIELDS))
        # Check values of first record
        self.assertEqual(recs[0]["DB"], "UniProtKB")
        self.assertEqual(recs[0]["DB_Object_ID"], "A0A024R1R8")
        self.assertEqual(recs[0]["DB_Object_Symbol"], "hCG_2014768")
        self.assertEqual(recs[0]["DB_Object_Name"], ["HCG2014768, isoform CRA_a"])
        self.assertEqual(recs[0]["DB_Object_Synonym"], ["hCG_2014768"])
        self.assertEqual(recs[0]["DB_Object_Type"], "protein")
        self.assertEqual(recs[0]["Taxon"], "taxon:9606")
        self.assertEqual(recs[0]["Parent_Object_ID"], "")
        self.assertEqual(recs[0]["DB_Xref"], [""])
        self.assertEqual(recs[0]["Gene_Product_Properties"], ["db_subset=TrEMBL"])

    def test_selection_writing(self):
        """Test record_has, and writerec.

        Adapted from Bio.UniProt.GOA.py by Iddo Friedberg idoerg@gmail.com.
        """
        recs = []
        filtered = []

        # Fields to filter
        evidence = {"Evidence": {"ND"}}
        synonym = {"Synonym": {"YA19A_YEAST", "YAL019W-A"}}
        taxon_id = {"Taxon_ID": {"taxon:559292"}}

        # Temporal file to test writerec
        f_number, f_filtered = tempfile.mkstemp()
        os.close(f_number)

        # Open a file and select records as per filter
        with open("UniProt/goa_yeast.gaf") as handle:
            for rec in GOA.gafiterator(handle):
                recs.append(rec)
                # Filtering
                if (
                    GOA.record_has(rec, taxon_id)
                    and GOA.record_has(rec, evidence)
                    and GOA.record_has(rec, synonym)
                ):
                    filtered.append(rec)

        # Check number of filtered records
        self.assertEqual(len(filtered), 3)

        # Write the filtered records to a file using writerec
        with open(f_filtered, "w") as handle:
            # '!gaf-version: 2.1'
            handle.write("!gaf-version: 2.1 \n")  # Adding file header
            for rec in filtered:
                GOA.writerec(rec, handle)

        # Open and read the file containing the filtered records
        recs_ff = []  # Records from filtered file
        with open(f_filtered) as handle:
            for rec in GOA.gafiterator(handle):
                recs_ff.append(rec)

        # Delete test file
        os.remove(f_filtered)

        # Compare, recs saved by writerec and filtered recs
        self.assertEqual(filtered, recs_ff)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
