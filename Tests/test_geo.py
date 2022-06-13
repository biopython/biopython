# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Tests the basic functionality of the GEO parsers."""

import unittest

from Bio import Geo


class TestGeo(unittest.TestCase):
    def test_soft_ex_dual(self):
        path = "Geo/soft_ex_dual.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(
                record.entity_id, "Control Embyronic Stem Cell Replicate 1"
            )
            self.assertEqual(len(record.entity_attributes), 24)
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch2"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "TriZol procedure",
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL3759")
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "Control Embyronic Stem Cell Replicate 1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"], "file1.gpr"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch2"], "Mus musculus"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"], "Mus musculus"
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "Cy5")
            self.assertEqual(len(record.entity_attributes["Sample_scan_protocol"]), 2)
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][0],
                "Scanned on an Agilent G2565AA scanner.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][1],
                "Images were quantified using Agilent Feature Extraction Software (version A.7.5).",
            )
            self.assertEqual(record.entity_attributes["sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch2"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal.",
            )
            self.assertEqual(record.entity_attributes["sample_table_end"], "")
            self.assertEqual(record.entity_attributes["Sample_label_ch2"], "Cy3")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Biological replicate 1 of 4. Control embryonic stem cells, untreated, harvested after several passages.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch2"],
                "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch2"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 4
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "ES-D3 cell line (CRL-1934)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1], "Age: day 4"
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][2],
                "Tissue: blastocytes",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][3],
                "Strain: 129/Sv mice",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch2"]), 3
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][0],
                "Strain: C57BL/6",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][1],
                "Age: e17.5 d",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][2],
                "Tissue: whole embryo",
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"],
                "log(REDsignal/GREENsignal) per feature (processed signals used).",
            )
            self.assertEqual(
                record.col_defs["gProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.',
            )
            self.assertEqual(
                record.col_defs["LogRatioError"],
                "error of the log ratio calculated according to the error model chosen.",
            )
            self.assertEqual(
                record.col_defs["PValueLogRatio"],
                "Significance level of the Log Ratio computed for a feature.",
            )
            self.assertEqual(
                record.col_defs["rProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LogRatioError")
            self.assertEqual(record.table_rows[0][3], "PValueLogRatio")
            self.assertEqual(record.table_rows[0][4], "gProcessedSignal")
            self.assertEqual(record.table_rows[0][5], "rProcessedSignal")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "-1.6274758")
            self.assertEqual(record.table_rows[1][2], "1.36E-01")
            self.assertEqual(record.table_rows[1][3], "6.41E-33")
            self.assertEqual(record.table_rows[1][4], "9.13E+03")
            self.assertEqual(record.table_rows[1][5], "2.15E+02")
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "0.1412248")
            self.assertEqual(record.table_rows[2][2], "1.34E+00")
            self.assertEqual(record.table_rows[2][3], "1.00E+00")
            self.assertEqual(record.table_rows[2][4], "4.14E+01")
            self.assertEqual(record.table_rows[2][5], "5.72E+01")
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.1827684")
            self.assertEqual(record.table_rows[3][2], "5.19E-02")
            self.assertEqual(record.table_rows[3][3], "4.33E-04")
            self.assertEqual(record.table_rows[3][4], "5.13E+03")
            self.assertEqual(record.table_rows[3][5], "7.81E+03")
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.3932267")
            self.assertEqual(record.table_rows[4][2], "6.08E-02")
            self.assertEqual(record.table_rows[4][3], "1.02E-10")
            self.assertEqual(record.table_rows[4][4], "4.65E+03")
            self.assertEqual(record.table_rows[4][5], "1.88E+03")
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "-0.9865994")
            self.assertEqual(record.table_rows[5][2], "1.05E-01")
            self.assertEqual(record.table_rows[5][3], "6.32E-21")
            self.assertEqual(record.table_rows[5][4], "2.91E+03")
            self.assertEqual(record.table_rows[5][5], "3.01E+02")
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "0.0238812")
            self.assertEqual(record.table_rows[6][2], "1.02E-01")
            self.assertEqual(record.table_rows[6][3], "8.15E-01")
            self.assertEqual(record.table_rows[6][4], "7.08E+02")
            self.assertEqual(record.table_rows[6][5], "7.48E+02")
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-1.4841822")
            self.assertEqual(record.table_rows[7][2], "1.25E-01")
            self.assertEqual(record.table_rows[7][3], "1.42E-32")
            self.assertEqual(record.table_rows[7][4], "1.02E+04")
            self.assertEqual(record.table_rows[7][5], "3.36E+02")
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "-1.8261356")
            self.assertEqual(record.table_rows[8][2], "4.15E-01")
            self.assertEqual(record.table_rows[8][3], "1.10E-05")
            self.assertEqual(record.table_rows[8][4], "7.19E+02")
            self.assertEqual(record.table_rows[8][5], "1.07E+01")
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "-1.0344779")
            self.assertEqual(record.table_rows[9][2], "1.78E+00")
            self.assertEqual(record.table_rows[9][3], "1.00E+00")
            self.assertEqual(record.table_rows[9][4], "9.62E+01")
            self.assertEqual(record.table_rows[9][5], "8.89E+00")
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "0.2405891")
            self.assertEqual(record.table_rows[10][2], "3.09E-01")
            self.assertEqual(record.table_rows[10][3], "4.36E-01")
            self.assertEqual(record.table_rows[10][4], "1.61E+02")
            self.assertEqual(record.table_rows[10][5], "2.80E+02")
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "0.3209366")
            self.assertEqual(record.table_rows[11][2], "3.59E-01")
            self.assertEqual(record.table_rows[11][3], "3.71E-01")
            self.assertEqual(record.table_rows[11][4], "1.25E+02")
            self.assertEqual(record.table_rows[11][5], "2.61E+02")
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "0.358304")
            self.assertEqual(record.table_rows[12][2], "2.06E+00")
            self.assertEqual(record.table_rows[12][3], "1.00E+00")
            self.assertEqual(record.table_rows[12][4], "2.04E+01")
            self.assertEqual(record.table_rows[12][5], "4.66E+01")
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "-0.0122072")
            self.assertEqual(record.table_rows[13][2], "3.64E-01")
            self.assertEqual(record.table_rows[13][3], "9.73E-01")
            self.assertEqual(record.table_rows[13][4], "1.84E+02")
            self.assertEqual(record.table_rows[13][5], "1.79E+02")
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "-1.5480396")
            self.assertEqual(record.table_rows[14][2], "1.30E-01")
            self.assertEqual(record.table_rows[14][3], "7.21E-33")
            self.assertEqual(record.table_rows[14][4], "1.02E+04")
            self.assertEqual(record.table_rows[14][5], "2.90E+02")
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "0.0073419")
            self.assertEqual(record.table_rows[15][2], "2.98E-01")
            self.assertEqual(record.table_rows[15][3], "9.80E-01")
            self.assertEqual(record.table_rows[15][4], "2.21E+02")
            self.assertEqual(record.table_rows[15][5], "2.25E+02")
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "-0.2267015")
            self.assertEqual(record.table_rows[16][2], "9.44E-01")
            self.assertEqual(record.table_rows[16][3], "8.10E-01")
            self.assertEqual(record.table_rows[16][4], "8.90E+01")
            self.assertEqual(record.table_rows[16][5], "5.28E+01")
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "-0.1484023")
            self.assertEqual(record.table_rows[17][2], "8.01E-01")
            self.assertEqual(record.table_rows[17][3], "8.53E-01")
            self.assertEqual(record.table_rows[17][4], "9.65E+01")
            self.assertEqual(record.table_rows[17][5], "6.86E+01")
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.6122195")
            self.assertEqual(record.table_rows[18][2], "1.28E-01")
            self.assertEqual(record.table_rows[18][3], "1.69E-06")
            self.assertEqual(record.table_rows[18][4], "1.12E+03")
            self.assertEqual(record.table_rows[18][5], "2.73E+02")
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.0796905")
            self.assertEqual(record.table_rows[19][2], "8.78E-02")
            self.assertEqual(record.table_rows[19][3], "3.64E-01")
            self.assertEqual(record.table_rows[19][4], "8.21E+02")
            self.assertEqual(record.table_rows[19][5], "9.87E+02")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "-0.084895")
            self.assertEqual(record.table_rows[20][2], "9.38E-01")
            self.assertEqual(record.table_rows[20][3], "9.28E-01")
            self.assertEqual(record.table_rows[20][4], "7.68E+01")
            self.assertEqual(record.table_rows[20][5], "6.32E+01")
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(
                record.entity_id, "Control Embyronic Stem Cell Replicate 2"
            )
            self.assertEqual(len(record.entity_attributes), 24)
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch2"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "TriZol procedure",
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL3759")
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "Control Embyronic Stem Cell Replicate 2",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"], "file2.gpr"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch2"], "Mus musculus"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"], "Mus musculus"
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "Cy5")
            self.assertEqual(len(record.entity_attributes["Sample_scan_protocol"]), 2)
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][0],
                "Scanned on an Agilent G2565AA scanner.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][1],
                "Images were quantified using Agilent Feature Extraction Software (version A.7.5).",
            )
            self.assertEqual(record.entity_attributes["sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch2"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal.",
            )
            self.assertEqual(record.entity_attributes["sample_table_end"], "")
            self.assertEqual(record.entity_attributes["Sample_label_ch2"], "Cy3")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Biological replicate 2 of 4. Control embryonic stem cells, untreated, harvested after several passages.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch2"],
                "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch2"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 4
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "ES-D3 cell line (CRL-1934)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1], "Age: day 4"
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][2],
                "Tissue: blastocytes",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][3],
                "Strain: 129/Sv mice",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch2"]), 3
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][0],
                "Strain: C57BL/6",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][1],
                "Age: e17.5 d",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][2],
                "Tissue: whole embryo",
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"],
                "log(REDsignal/GREENsignal) per feature (processed signals used).",
            )
            self.assertEqual(
                record.col_defs["gProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.',
            )
            self.assertEqual(
                record.col_defs["LogRatioError"],
                "error of the log ratio calculated according to the error model chosen.",
            )
            self.assertEqual(
                record.col_defs["PValueLogRatio"],
                "Significance level of the Log Ratio computed for a feature.",
            )
            self.assertEqual(
                record.col_defs["rProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LogRatioError")
            self.assertEqual(record.table_rows[0][3], "PValueLogRatio")
            self.assertEqual(record.table_rows[0][4], "gProcessedSignal")
            self.assertEqual(record.table_rows[0][5], "rProcessedSignal")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "-1.1697263")
            self.assertEqual(record.table_rows[1][2], "1.23E-01")
            self.assertEqual(record.table_rows[1][3], "2.14E-21")
            self.assertEqual(record.table_rows[1][4], "3.17E+03")
            self.assertEqual(record.table_rows[1][5], "2.14E+02")
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "-0.1111353")
            self.assertEqual(record.table_rows[2][2], "1.63E+00")
            self.assertEqual(record.table_rows[2][3], "9.46E-01")
            self.assertEqual(record.table_rows[2][4], "5.43E+01")
            self.assertEqual(record.table_rows[2][5], "4.20E+01")
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.1400597")
            self.assertEqual(record.table_rows[3][2], "5.11E-02")
            self.assertEqual(record.table_rows[3][3], "6.17E-03")
            self.assertEqual(record.table_rows[3][4], "6.72E+03")
            self.assertEqual(record.table_rows[3][5], "9.28E+03")
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.4820633")
            self.assertEqual(record.table_rows[4][2], "6.38E-02")
            self.assertEqual(record.table_rows[4][3], "4.06E-14")
            self.assertEqual(record.table_rows[4][4], "6.46E+03")
            self.assertEqual(record.table_rows[4][5], "2.13E+03")
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "-1.2116196")
            self.assertEqual(record.table_rows[5][2], "1.22E-01")
            self.assertEqual(record.table_rows[5][3], "2.31E-23")
            self.assertEqual(record.table_rows[5][4], "3.62E+03")
            self.assertEqual(record.table_rows[5][5], "2.22E+02")
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "-0.0230528")
            self.assertEqual(record.table_rows[6][2], "1.04E-01")
            self.assertEqual(record.table_rows[6][3], "8.24E-01")
            self.assertEqual(record.table_rows[6][4], "8.76E+02")
            self.assertEqual(record.table_rows[6][5], "8.31E+02")
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-1.1380152")
            self.assertEqual(record.table_rows[7][2], "1.13E-01")
            self.assertEqual(record.table_rows[7][3], "9.23E-24")
            self.assertEqual(record.table_rows[7][4], "3.94E+03")
            self.assertEqual(record.table_rows[7][5], "2.86E+02")
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "-1.834596")
            self.assertEqual(record.table_rows[8][2], "5.40E-01")
            self.assertEqual(record.table_rows[8][3], "6.74E-04")
            self.assertEqual(record.table_rows[8][4], "6.44E+02")
            self.assertEqual(record.table_rows[8][5], "9.43E+00")
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "-0.9747637")
            self.assertEqual(record.table_rows[9][2], "2.14E+00")
            self.assertEqual(record.table_rows[9][3], "1.00E+00")
            self.assertEqual(record.table_rows[9][4], "9.17E+01")
            self.assertEqual(record.table_rows[9][5], "9.72E+00")
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "0.3874005")
            self.assertEqual(record.table_rows[10][2], "2.92E-01")
            self.assertEqual(record.table_rows[10][3], "1.85E-01")
            self.assertEqual(record.table_rows[10][4], "1.69E+02")
            self.assertEqual(record.table_rows[10][5], "4.11E+02")
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "0.5340442")
            self.assertEqual(record.table_rows[11][2], "3.29E-01")
            self.assertEqual(record.table_rows[11][3], "1.04E-01")
            self.assertEqual(record.table_rows[11][4], "1.23E+02")
            self.assertEqual(record.table_rows[11][5], "4.20E+02")
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "0.3260696")
            self.assertEqual(record.table_rows[12][2], "1.92E+00")
            self.assertEqual(record.table_rows[12][3], "8.65E-01")
            self.assertEqual(record.table_rows[12][4], "2.73E+01")
            self.assertEqual(record.table_rows[12][5], "5.77E+01")
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "0.3010618")
            self.assertEqual(record.table_rows[13][2], "2.84E-01")
            self.assertEqual(record.table_rows[13][3], "2.90E-01")
            self.assertEqual(record.table_rows[13][4], "1.93E+02")
            self.assertEqual(record.table_rows[13][5], "3.87E+02")
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "-1.0760413")
            self.assertEqual(record.table_rows[14][2], "1.08E-01")
            self.assertEqual(record.table_rows[14][3], "1.63E-23")
            self.assertEqual(record.table_rows[14][4], "4.06E+03")
            self.assertEqual(record.table_rows[14][5], "3.41E+02")
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "-0.1167371")
            self.assertEqual(record.table_rows[15][2], "3.87E-01")
            self.assertEqual(record.table_rows[15][3], "7.63E-01")
            self.assertEqual(record.table_rows[15][4], "2.32E+02")
            self.assertEqual(record.table_rows[15][5], "1.77E+02")
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "-0.1936322")
            self.assertEqual(record.table_rows[16][2], "9.44E-01")
            self.assertEqual(record.table_rows[16][3], "8.38E-01")
            self.assertEqual(record.table_rows[16][4], "1.02E+02")
            self.assertEqual(record.table_rows[16][5], "6.56E+01")
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "-0.3275898")
            self.assertEqual(record.table_rows[17][2], "7.87E-01")
            self.assertEqual(record.table_rows[17][3], "6.77E-01")
            self.assertEqual(record.table_rows[17][4], "1.41E+02")
            self.assertEqual(record.table_rows[17][5], "6.65E+01")
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.4805853")
            self.assertEqual(record.table_rows[18][2], "1.14E-01")
            self.assertEqual(record.table_rows[18][3], "2.41E-05")
            self.assertEqual(record.table_rows[18][4], "1.34E+03")
            self.assertEqual(record.table_rows[18][5], "4.42E+02")
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.1109524")
            self.assertEqual(record.table_rows[19][2], "9.56E-02")
            self.assertEqual(record.table_rows[19][3], "2.46E-01")
            self.assertEqual(record.table_rows[19][4], "8.38E+02")
            self.assertEqual(record.table_rows[19][5], "1.08E+03")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "0.1677912")
            self.assertEqual(record.table_rows[20][2], "6.51E-01")
            self.assertEqual(record.table_rows[20][3], "7.97E-01")
            self.assertEqual(record.table_rows[20][4], "9.84E+01")
            self.assertEqual(record.table_rows[20][5], "1.45E+02")
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(
                record.entity_id,
                "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1",
            )
            self.assertEqual(len(record.entity_attributes), 25)
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch2"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "TriZol procedure",
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL3759")
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"], "file3.gpr"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch2"], "Mus musculus"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"], "Mus musculus"
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "Cy5")
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "Scanned on an Agilent G2565AA scanner.",
            )
            self.assertEqual(record.entity_attributes["sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch2"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal.",
            )
            self.assertEqual(record.entity_attributes["sample_table_end"], "")
            self.assertEqual(record.entity_attributes["Sample_label_ch2"], "Cy3")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Biological replicate 1 of 3. Stable triple-fusion-reporter-gene transfected embryonic stem cells, harvested after several passages.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Total RNA from murine ES-D3 triple-transfected embryonic stem cells labeled with Cyanine-5 (red).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch2"],
                "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch2"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 5
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "ES-D3 cell line (CRL-1934)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Transfected with pUb-fluc-mrfp-ttk triple fusion reporter gene.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][2], "Age: day 4"
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][3],
                "Tissue: blastocytes",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][4],
                "Strain: 129/Sv mice",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch2"]), 3
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][0],
                "Strain: C57BL/6",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][1],
                "Age: e17.5 d",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][2],
                "Tissue: whole embryo",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "PCR amplification and standard cloning techniques were used to insert fluc and mrfp genes from plasmids pCDNA 3.1-CMV-fluc (Promega, Madison, WI) and pCDNA3.1-CMV-mrfp in frame with the ttk gene into the pCDNA3.1-truncated	sr39tk. This triple fusion (TF) reporter gene fragment (3.3 kbp) was released from the plasmid with Not1 and BamH1 restriction enzymes before blunt-end ligation into the multiple cloning site of lentiviral transfer vector, FUG, driven by the human ubiquitin-C promoter. Self-inactivating (SIN) lentivirus was prepared by transient transfection of 293T cells. Briefly, pFUG-TF containing the triple fusion reporter gene was co-transfected into 293T cells with HIV-1 packaging vector (?8.9) and vesicular stomatitis virus G glycoprotein-pseudotyped envelop vector (pVSVG). Lentivirus supernatant was concentrated by sediment centrifugation using a SW29 rotor at 50,000 x g for two hours. Concentrated virus was titered on 293T cells. Murine ES cells were transfected with LV-pUb-fluc-mrfp-ttk at a multiplicity of infection (MOI) of 10.",
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"],
                "log(REDsignal/GREENsignal) per feature (processed signals used).",
            )
            self.assertEqual(
                record.col_defs["gProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.',
            )
            self.assertEqual(
                record.col_defs["LogRatioError"],
                "error of the log ratio calculated according to the error model chosen.",
            )
            self.assertEqual(
                record.col_defs["PValueLogRatio"],
                "Significance level of the Log Ratio computed for a feature.",
            )
            self.assertEqual(
                record.col_defs["rProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LogRatioError")
            self.assertEqual(record.table_rows[0][3], "PValueLogRatio")
            self.assertEqual(record.table_rows[0][4], "gProcessedSignal")
            self.assertEqual(record.table_rows[0][5], "rProcessedSignal")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "-0.7837546")
            self.assertEqual(record.table_rows[1][2], "1.30E-01")
            self.assertEqual(record.table_rows[1][3], "1.70E-09")
            self.assertEqual(record.table_rows[1][4], "2.10E+03")
            self.assertEqual(record.table_rows[1][5], "3.46E+02")
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "0.3797837")
            self.assertEqual(record.table_rows[2][2], "1.15E+00")
            self.assertEqual(record.table_rows[2][3], "7.41E-01")
            self.assertEqual(record.table_rows[2][4], "5.59E+01")
            self.assertEqual(record.table_rows[2][5], "1.34E+02")
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.2079269")
            self.assertEqual(record.table_rows[3][2], "5.38E-02")
            self.assertEqual(record.table_rows[3][3], "1.12E-04")
            self.assertEqual(record.table_rows[3][4], "5.04E+03")
            self.assertEqual(record.table_rows[3][5], "8.14E+03")
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.4730291")
            self.assertEqual(record.table_rows[4][2], "6.71E-02")
            self.assertEqual(record.table_rows[4][3], "1.86E-12")
            self.assertEqual(record.table_rows[4][4], "5.66E+03")
            self.assertEqual(record.table_rows[4][5], "1.91E+03")
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "-0.9481128")
            self.assertEqual(record.table_rows[5][2], "1.19E-01")
            self.assertEqual(record.table_rows[5][3], "1.30E-15")
            self.assertEqual(record.table_rows[5][4], "3.10E+03")
            self.assertEqual(record.table_rows[5][5], "3.49E+02")
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "-0.0159867")
            self.assertEqual(record.table_rows[6][2], "1.33E-01")
            self.assertEqual(record.table_rows[6][3], "9.05E-01")
            self.assertEqual(record.table_rows[6][4], "8.45E+02")
            self.assertEqual(record.table_rows[6][5], "8.14E+02")
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-0.819922")
            self.assertEqual(record.table_rows[7][2], "1.14E-01")
            self.assertEqual(record.table_rows[7][3], "7.01E-13")
            self.assertEqual(record.table_rows[7][4], "2.75E+03")
            self.assertEqual(record.table_rows[7][5], "4.16E+02")
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "-0.1559774")
            self.assertEqual(record.table_rows[8][2], "9.16E-01")
            self.assertEqual(record.table_rows[8][3], "8.65E-01")
            self.assertEqual(record.table_rows[8][4], "1.34E+02")
            self.assertEqual(record.table_rows[8][5], "9.34E+01")
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "0.145267")
            self.assertEqual(record.table_rows[9][2], "3.90E+00")
            self.assertEqual(record.table_rows[9][3], "1.00E+00")
            self.assertEqual(record.table_rows[9][4], "2.22E+01")
            self.assertEqual(record.table_rows[9][5], "3.10E+01")
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "0.3611211")
            self.assertEqual(record.table_rows[10][2], "3.40E-01")
            self.assertEqual(record.table_rows[10][3], "2.88E-01")
            self.assertEqual(record.table_rows[10][4], "1.97E+02")
            self.assertEqual(record.table_rows[10][5], "4.52E+02")
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "0.5092089")
            self.assertEqual(record.table_rows[11][2], "4.39E-01")
            self.assertEqual(record.table_rows[11][3], "2.46E-01")
            self.assertEqual(record.table_rows[11][4], "1.24E+02")
            self.assertEqual(record.table_rows[11][5], "4.01E+02")
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "0.3715387")
            self.assertEqual(record.table_rows[12][2], "1.69E+00")
            self.assertEqual(record.table_rows[12][3], "8.26E-01")
            self.assertEqual(record.table_rows[12][4], "3.84E+01")
            self.assertEqual(record.table_rows[12][5], "9.04E+01")
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "0.1734934")
            self.assertEqual(record.table_rows[13][2], "3.57E-01")
            self.assertEqual(record.table_rows[13][3], "6.27E-01")
            self.assertEqual(record.table_rows[13][4], "2.37E+02")
            self.assertEqual(record.table_rows[13][5], "3.53E+02")
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "-0.9340707")
            self.assertEqual(record.table_rows[14][2], "1.20E-01")
            self.assertEqual(record.table_rows[14][3], "6.90E-15")
            self.assertEqual(record.table_rows[14][4], "2.96E+03")
            self.assertEqual(record.table_rows[14][5], "3.45E+02")
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "-0.2956317")
            self.assertEqual(record.table_rows[15][2], "5.78E-01")
            self.assertEqual(record.table_rows[15][3], "6.09E-01")
            self.assertEqual(record.table_rows[15][4], "2.46E+02")
            self.assertEqual(record.table_rows[15][5], "1.25E+02")
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "-0.2321102")
            self.assertEqual(record.table_rows[16][2], "1.22E+00")
            self.assertEqual(record.table_rows[16][3], "8.49E-01")
            self.assertEqual(record.table_rows[16][4], "1.09E+02")
            self.assertEqual(record.table_rows[16][5], "6.37E+01")
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "-0.1603561")
            self.assertEqual(record.table_rows[17][2], "1.16E+00")
            self.assertEqual(record.table_rows[17][3], "8.90E-01")
            self.assertEqual(record.table_rows[17][4], "1.06E+02")
            self.assertEqual(record.table_rows[17][5], "7.34E+01")
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.5063897")
            self.assertEqual(record.table_rows[18][2], "1.63E-01")
            self.assertEqual(record.table_rows[18][3], "1.95E-03")
            self.assertEqual(record.table_rows[18][4], "1.15E+03")
            self.assertEqual(record.table_rows[18][5], "3.58E+02")
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.1990761")
            self.assertEqual(record.table_rows[19][2], "1.32E-01")
            self.assertEqual(record.table_rows[19][3], "1.32E-01")
            self.assertEqual(record.table_rows[19][4], "6.65E+02")
            self.assertEqual(record.table_rows[19][5], "1.05E+03")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "0.2985912")
            self.assertEqual(record.table_rows[20][2], "8.89E-01")
            self.assertEqual(record.table_rows[20][3], "7.37E-01")
            self.assertEqual(record.table_rows[20][4], "8.06E+01")
            self.assertEqual(record.table_rows[20][5], "1.60E+02")

    def test_soft_ex_affy(self):
        path = "Geo/soft_ex_affy.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "Drosophila_T0-1")
            self.assertEqual(len(record.entity_attributes), 18)
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"],
                "Drosophila melanogaster",
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "biotin")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation.",
            )
            self.assertEqual(record.entity_attributes["Sample_table_end"], "")
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 2
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "Genotype: yellow white and Oregon R parents",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "Trizol extraction of total RNA was performed according to the manufacturer's instructions.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Drosophila embryos before nuclear cycle 9 (maternal transcripts)",
            )
            self.assertEqual(record.entity_attributes["Sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "embryo at T0, biological rep1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"],
                "Drosophila_T0-1.CEL",
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL72")
            self.assertEqual(len(record.col_defs), 4)
            self.assertEqual(
                record.col_defs["DETECTION P-VALUE"],
                "'detection p-value', p-value that indicates the significance level of the detection call",
            )
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"], "MAS5-calculated Signal intensity"
            )
            self.assertEqual(
                record.col_defs["ABS_CALL"],
                "the call in an absolute analysis that indicates if the transcript was present (P), absent (A), marginal (M), or no call (NC)",
            )
            self.assertEqual(len(record.table_rows), 22)
            self.assertEqual(len(record.table_rows[0]), 4)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "ABS_CALL")
            self.assertEqual(record.table_rows[0][3], "DETECTION P-VALUE")
            self.assertEqual(len(record.table_rows[1]), 4)
            self.assertEqual(record.table_rows[1][0], "141200_at")
            self.assertEqual(record.table_rows[1][1], "36.6")
            self.assertEqual(record.table_rows[1][2], "A")
            self.assertEqual(record.table_rows[1][3], "0.818657")
            self.assertEqual(len(record.table_rows[2]), 4)
            self.assertEqual(record.table_rows[2][0], "141201_at")
            self.assertEqual(record.table_rows[2][1], "41.5")
            self.assertEqual(record.table_rows[2][2], "A")
            self.assertEqual(record.table_rows[2][3], "0.703191")
            self.assertEqual(len(record.table_rows[3]), 4)
            self.assertEqual(record.table_rows[3][0], "141202_at")
            self.assertEqual(record.table_rows[3][1], "607.3")
            self.assertEqual(record.table_rows[3][2], "P")
            self.assertEqual(record.table_rows[3][3], "0.000944")
            self.assertEqual(len(record.table_rows[4]), 4)
            self.assertEqual(record.table_rows[4][0], "141203_at")
            self.assertEqual(record.table_rows[4][1], "1509.1")
            self.assertEqual(record.table_rows[4][2], "P")
            self.assertEqual(record.table_rows[4][3], "0.000762")
            self.assertEqual(len(record.table_rows[5]), 4)
            self.assertEqual(record.table_rows[5][0], "141204_at")
            self.assertEqual(record.table_rows[5][1], "837.3")
            self.assertEqual(record.table_rows[5][2], "P")
            self.assertEqual(record.table_rows[5][3], "0.000613")
            self.assertEqual(len(record.table_rows[6]), 4)
            self.assertEqual(record.table_rows[6][0], "141205_at")
            self.assertEqual(record.table_rows[6][1], "363.2")
            self.assertEqual(record.table_rows[6][2], "P")
            self.assertEqual(record.table_rows[6][3], "0.003815")
            self.assertEqual(len(record.table_rows[7]), 4)
            self.assertEqual(record.table_rows[7][0], "141206_at")
            self.assertEqual(record.table_rows[7][1], "1193.6")
            self.assertEqual(record.table_rows[7][2], "P")
            self.assertEqual(record.table_rows[7][3], "0.000491")
            self.assertEqual(len(record.table_rows[8]), 4)
            self.assertEqual(record.table_rows[8][0], "141207_at")
            self.assertEqual(record.table_rows[8][1], "346.6")
            self.assertEqual(record.table_rows[8][2], "P")
            self.assertEqual(record.table_rows[8][3], "0.001165")
            self.assertEqual(len(record.table_rows[9]), 4)
            self.assertEqual(record.table_rows[9][0], "141208_at")
            self.assertEqual(record.table_rows[9][1], "257.8")
            self.assertEqual(record.table_rows[9][2], "P")
            self.assertEqual(record.table_rows[9][3], "0.006575")
            self.assertEqual(len(record.table_rows[10]), 4)
            self.assertEqual(record.table_rows[10][0], "141209_at")
            self.assertEqual(record.table_rows[10][1], "337.1")
            self.assertEqual(record.table_rows[10][2], "P")
            self.assertEqual(record.table_rows[10][3], "0.002607")
            self.assertEqual(len(record.table_rows[11]), 4)
            self.assertEqual(record.table_rows[11][0], "141210_at")
            self.assertEqual(record.table_rows[11][1], "48")
            self.assertEqual(record.table_rows[11][2], "A")
            self.assertEqual(record.table_rows[11][3], "0.150145")
            self.assertEqual(len(record.table_rows[12]), 4)
            self.assertEqual(record.table_rows[12][0], "141211_at")
            self.assertEqual(record.table_rows[12][1], "130.7")
            self.assertEqual(record.table_rows[12][2], "P")
            self.assertEqual(record.table_rows[12][3], "0.005504")
            self.assertEqual(len(record.table_rows[13]), 4)
            self.assertEqual(record.table_rows[13][0], "141212_at")
            self.assertEqual(record.table_rows[13][1], "1454.3")
            self.assertEqual(record.table_rows[13][2], "P")
            self.assertEqual(record.table_rows[13][3], "0.000491")
            self.assertEqual(len(record.table_rows[14]), 4)
            self.assertEqual(record.table_rows[14][0], "141213_at")
            self.assertEqual(record.table_rows[14][1], "21.2")
            self.assertEqual(record.table_rows[14][2], "A")
            self.assertEqual(record.table_rows[14][3], "0.635055")
            self.assertEqual(len(record.table_rows[15]), 4)
            self.assertEqual(record.table_rows[15][0], "142121_at")
            self.assertEqual(record.table_rows[15][1], "133.7")
            self.assertEqual(record.table_rows[15][2], "A")
            self.assertEqual(record.table_rows[15][3], "0.889551")
            self.assertEqual(len(record.table_rows[16]), 4)
            self.assertEqual(record.table_rows[16][0], "142122_at")
            self.assertEqual(record.table_rows[16][1], "275.3")
            self.assertEqual(record.table_rows[16][2], "A")
            self.assertEqual(record.table_rows[16][3], "0.611218")
            self.assertEqual(len(record.table_rows[17]), 4)
            self.assertEqual(record.table_rows[17][0], "142123_at")
            self.assertEqual(record.table_rows[17][1], "307.6")
            self.assertEqual(record.table_rows[17][2], "A")
            self.assertEqual(record.table_rows[17][3], "0.611218")
            self.assertEqual(len(record.table_rows[18]), 4)
            self.assertEqual(record.table_rows[18][0], "142124_at")
            self.assertEqual(record.table_rows[18][1], "132.6")
            self.assertEqual(record.table_rows[18][2], "A")
            self.assertEqual(record.table_rows[18][3], "0.437646")
            self.assertEqual(len(record.table_rows[19]), 4)
            self.assertEqual(record.table_rows[19][0], "142125_at")
            self.assertEqual(record.table_rows[19][1], "195.8")
            self.assertEqual(record.table_rows[19][2], "A")
            self.assertEqual(record.table_rows[19][3], "0.110449")
            self.assertEqual(len(record.table_rows[20]), 4)
            self.assertEqual(record.table_rows[20][0], "142126_at")
            self.assertEqual(record.table_rows[20][1], "174.1")
            self.assertEqual(record.table_rows[20][2], "A")
            self.assertEqual(record.table_rows[20][3], "0.681117")
            self.assertEqual(len(record.table_rows[21]), 4)
            self.assertEqual(record.table_rows[21][0], "142127_at")
            self.assertEqual(record.table_rows[21][1], "316.3")
            self.assertEqual(record.table_rows[21][2], "A")
            self.assertEqual(record.table_rows[21][3], "0.65838")
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "Drosophila_T0-2")
            self.assertEqual(len(record.entity_attributes), 18)
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"],
                "Drosophila melanogaster",
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "biotin")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation.",
            )
            self.assertEqual(record.entity_attributes["Sample_table_end"], "")
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 2
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "Genotype: yellow white and Oregon R parents",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "Trizol extraction of total RNA was performed according to the manufacturer's instructions.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Drosophila embryos before nuclear cycle 9 (maternal transcripts)",
            )
            self.assertEqual(record.entity_attributes["Sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "embryo at T0, biological rep2",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"],
                "Drosophila_T0-2.CEL",
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL72")
            self.assertEqual(len(record.col_defs), 4)
            self.assertEqual(
                record.col_defs["DETECTION P-VALUE"],
                "'detection p-value', p-value that indicates the significance level of the detection call",
            )
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"], "MAS5-calculated Signal intensity"
            )
            self.assertEqual(
                record.col_defs["ABS_CALL"],
                "the call in an absolute analysis that indicates if the transcript was present (P), absent (A), marginal (M), or no call (NC)",
            )
            self.assertEqual(len(record.table_rows), 22)
            self.assertEqual(len(record.table_rows[0]), 4)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "ABS_CALL")
            self.assertEqual(record.table_rows[0][3], "DETECTION P-VALUE")
            self.assertEqual(len(record.table_rows[1]), 4)
            self.assertEqual(record.table_rows[1][0], "141200_at")
            self.assertEqual(record.table_rows[1][1], "70.3")
            self.assertEqual(record.table_rows[1][2], "A")
            self.assertEqual(record.table_rows[1][3], "0.216313")
            self.assertEqual(len(record.table_rows[2]), 4)
            self.assertEqual(record.table_rows[2][0], "141201_at")
            self.assertEqual(record.table_rows[2][1], "38")
            self.assertEqual(record.table_rows[2][2], "A")
            self.assertEqual(record.table_rows[2][3], "0.635055")
            self.assertEqual(len(record.table_rows[3]), 4)
            self.assertEqual(record.table_rows[3][0], "141202_at")
            self.assertEqual(record.table_rows[3][1], "831.8")
            self.assertEqual(record.table_rows[3][2], "P")
            self.assertEqual(record.table_rows[3][3], "0.000613")
            self.assertEqual(len(record.table_rows[4]), 4)
            self.assertEqual(record.table_rows[4][0], "141203_at")
            self.assertEqual(record.table_rows[4][1], "2215.5")
            self.assertEqual(record.table_rows[4][2], "P")
            self.assertEqual(record.table_rows[4][3], "0.000944")
            self.assertEqual(len(record.table_rows[5]), 4)
            self.assertEqual(record.table_rows[5][0], "141204_at")
            self.assertEqual(record.table_rows[5][1], "965.6")
            self.assertEqual(record.table_rows[5][2], "P")
            self.assertEqual(record.table_rows[5][3], "0.000491")
            self.assertEqual(len(record.table_rows[6]), 4)
            self.assertEqual(record.table_rows[6][0], "141205_at")
            self.assertEqual(record.table_rows[6][1], "383.2")
            self.assertEqual(record.table_rows[6][2], "P")
            self.assertEqual(record.table_rows[6][3], "0.001432")
            self.assertEqual(len(record.table_rows[7]), 4)
            self.assertEqual(record.table_rows[7][0], "141206_at")
            self.assertEqual(record.table_rows[7][1], "1195")
            self.assertEqual(record.table_rows[7][2], "P")
            self.assertEqual(record.table_rows[7][3], "0.000491")
            self.assertEqual(len(record.table_rows[8]), 4)
            self.assertEqual(record.table_rows[8][0], "141207_at")
            self.assertEqual(record.table_rows[8][1], "413.7")
            self.assertEqual(record.table_rows[8][2], "P")
            self.assertEqual(record.table_rows[8][3], "0.000613")
            self.assertEqual(len(record.table_rows[9]), 4)
            self.assertEqual(record.table_rows[9][0], "141208_at")
            self.assertEqual(record.table_rows[9][1], "447.3")
            self.assertEqual(record.table_rows[9][2], "P")
            self.assertEqual(record.table_rows[9][3], "0.000762")
            self.assertEqual(len(record.table_rows[10]), 4)
            self.assertEqual(record.table_rows[10][0], "141209_at")
            self.assertEqual(record.table_rows[10][1], "294.4")
            self.assertEqual(record.table_rows[10][2], "P")
            self.assertEqual(record.table_rows[10][3], "0.004591")
            self.assertEqual(len(record.table_rows[11]), 4)
            self.assertEqual(record.table_rows[11][0], "141210_at")
            self.assertEqual(record.table_rows[11][1], "81.7")
            self.assertEqual(record.table_rows[11][2], "M")
            self.assertEqual(record.table_rows[11][3], "0.054711")
            self.assertEqual(len(record.table_rows[12]), 4)
            self.assertEqual(record.table_rows[12][0], "141211_at")
            self.assertEqual(record.table_rows[12][1], "84.9")
            self.assertEqual(record.table_rows[12][2], "P")
            self.assertEqual(record.table_rows[12][3], "0.005504")
            self.assertEqual(len(record.table_rows[13]), 4)
            self.assertEqual(record.table_rows[13][0], "141212_at")
            self.assertEqual(record.table_rows[13][1], "1456.4")
            self.assertEqual(record.table_rows[13][2], "P")
            self.assertEqual(record.table_rows[13][3], "0.000491")
            self.assertEqual(len(record.table_rows[14]), 4)
            self.assertEqual(record.table_rows[14][0], "141213_at")
            self.assertEqual(record.table_rows[14][1], "37")
            self.assertEqual(record.table_rows[14][2], "A")
            self.assertEqual(record.table_rows[14][3], "0.122747")
            self.assertEqual(len(record.table_rows[15]), 4)
            self.assertEqual(record.table_rows[15][0], "142121_at")
            self.assertEqual(record.table_rows[15][1], "133.7")
            self.assertEqual(record.table_rows[15][2], "A")
            self.assertEqual(record.table_rows[15][3], "0.889551")
            self.assertEqual(len(record.table_rows[16]), 4)
            self.assertEqual(record.table_rows[16][0], "142122_at")
            self.assertEqual(record.table_rows[16][1], "275.3")
            self.assertEqual(record.table_rows[16][2], "A")
            self.assertEqual(record.table_rows[16][3], "0.611218")
            self.assertEqual(len(record.table_rows[17]), 4)
            self.assertEqual(record.table_rows[17][0], "142123_at")
            self.assertEqual(record.table_rows[17][1], "307.6")
            self.assertEqual(record.table_rows[17][2], "A")
            self.assertEqual(record.table_rows[17][3], "0.611218")
            self.assertEqual(len(record.table_rows[18]), 4)
            self.assertEqual(record.table_rows[18][0], "142124_at")
            self.assertEqual(record.table_rows[18][1], "132.6")
            self.assertEqual(record.table_rows[18][2], "A")
            self.assertEqual(record.table_rows[18][3], "0.437646")
            self.assertEqual(len(record.table_rows[19]), 4)
            self.assertEqual(record.table_rows[19][0], "142125_at")
            self.assertEqual(record.table_rows[19][1], "195.8")
            self.assertEqual(record.table_rows[19][2], "A")
            self.assertEqual(record.table_rows[19][3], "0.110449")
            self.assertEqual(len(record.table_rows[20]), 4)
            self.assertEqual(record.table_rows[20][0], "142126_at")
            self.assertEqual(record.table_rows[20][1], "174.1")
            self.assertEqual(record.table_rows[20][2], "A")
            self.assertEqual(record.table_rows[20][3], "0.681117")
            self.assertEqual(len(record.table_rows[21]), 4)
            self.assertEqual(record.table_rows[21][0], "142127_at")
            self.assertEqual(record.table_rows[21][1], "316.3")
            self.assertEqual(record.table_rows[21][2], "A")
            self.assertEqual(record.table_rows[21][3], "0.65838")
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "Drosophila_T1-1")
            self.assertEqual(len(record.entity_attributes), 18)
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"],
                "Drosophila melanogaster",
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "biotin")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Gene expression data from embryos in slow phase of cellularisation.",
            )
            self.assertEqual(record.entity_attributes["Sample_table_end"], "")
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 2
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "Genotype: yellow white and Oregon R parents",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Age: embryos in slow phase of cellularisation",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "Trizol extraction of total RNA was performed according to the manufacturer's instructions.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Drosophila embryos in slow phase of cellularisation",
            )
            self.assertEqual(record.entity_attributes["Sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "embryo at T1, biological rep1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"],
                "Drosophila_T1-1.CEL",
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL72")
            self.assertEqual(len(record.col_defs), 4)
            self.assertEqual(
                record.col_defs["DETECTION P-VALUE"],
                "'detection p-value', p-value that indicates the significance level of the detection call",
            )
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"], "MAS5-calculated Signal intensity"
            )
            self.assertEqual(
                record.col_defs["ABS_CALL"],
                "the call in an absolute analysis that indicates if the transcript was present (P), absent (A), marginal (M), or no call (NC)",
            )
            self.assertEqual(len(record.table_rows), 22)
            self.assertEqual(len(record.table_rows[0]), 4)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "ABS_CALL")
            self.assertEqual(record.table_rows[0][3], "DETECTION P-VALUE")
            self.assertEqual(len(record.table_rows[1]), 4)
            self.assertEqual(record.table_rows[1][0], "141200_at")
            self.assertEqual(record.table_rows[1][1], "20.8")
            self.assertEqual(record.table_rows[1][2], "A")
            self.assertEqual(record.table_rows[1][3], "0.801637")
            self.assertEqual(len(record.table_rows[2]), 4)
            self.assertEqual(record.table_rows[2][0], "141201_at")
            self.assertEqual(record.table_rows[2][1], "85.8")
            self.assertEqual(record.table_rows[2][2], "A")
            self.assertEqual(record.table_rows[2][3], "0.48748")
            self.assertEqual(len(record.table_rows[3]), 4)
            self.assertEqual(record.table_rows[3][0], "141202_at")
            self.assertEqual(record.table_rows[3][1], "704.8")
            self.assertEqual(record.table_rows[3][2], "P")
            self.assertEqual(record.table_rows[3][3], "0.000613")
            self.assertEqual(len(record.table_rows[4]), 4)
            self.assertEqual(record.table_rows[4][0], "141203_at")
            self.assertEqual(record.table_rows[4][1], "1036.6")
            self.assertEqual(record.table_rows[4][2], "P")
            self.assertEqual(record.table_rows[4][3], "0.000944")
            self.assertEqual(len(record.table_rows[5]), 4)
            self.assertEqual(record.table_rows[5][0], "141204_at")
            self.assertEqual(record.table_rows[5][1], "700.3")
            self.assertEqual(record.table_rows[5][2], "P")
            self.assertEqual(record.table_rows[5][3], "0.000491")
            self.assertEqual(len(record.table_rows[6]), 4)
            self.assertEqual(record.table_rows[6][0], "141205_at")
            self.assertEqual(record.table_rows[6][1], "462.4")
            self.assertEqual(record.table_rows[6][2], "P")
            self.assertEqual(record.table_rows[6][3], "0.003159")
            self.assertEqual(len(record.table_rows[7]), 4)
            self.assertEqual(record.table_rows[7][0], "141206_at")
            self.assertEqual(record.table_rows[7][1], "1301.9")
            self.assertEqual(record.table_rows[7][2], "P")
            self.assertEqual(record.table_rows[7][3], "0.000491")
            self.assertEqual(len(record.table_rows[8]), 4)
            self.assertEqual(record.table_rows[8][0], "141207_at")
            self.assertEqual(record.table_rows[8][1], "454.8")
            self.assertEqual(record.table_rows[8][2], "P")
            self.assertEqual(record.table_rows[8][3], "0.000944")
            self.assertEqual(len(record.table_rows[9]), 4)
            self.assertEqual(record.table_rows[9][0], "141208_at")
            self.assertEqual(record.table_rows[9][1], "438.6")
            self.assertEqual(record.table_rows[9][2], "P")
            self.assertEqual(record.table_rows[9][3], "0.000944")
            self.assertEqual(len(record.table_rows[10]), 4)
            self.assertEqual(record.table_rows[10][0], "141209_at")
            self.assertEqual(record.table_rows[10][1], "264.4")
            self.assertEqual(record.table_rows[10][2], "P")
            self.assertEqual(record.table_rows[10][3], "0.004591")
            self.assertEqual(len(record.table_rows[11]), 4)
            self.assertEqual(record.table_rows[11][0], "141210_at")
            self.assertEqual(record.table_rows[11][1], "65.6")
            self.assertEqual(record.table_rows[11][2], "A")
            self.assertEqual(record.table_rows[11][3], "0.150145")
            self.assertEqual(len(record.table_rows[12]), 4)
            self.assertEqual(record.table_rows[12][0], "141211_at")
            self.assertEqual(record.table_rows[12][1], "72.2")
            self.assertEqual(record.table_rows[12][2], "A")
            self.assertEqual(record.table_rows[12][3], "0.070073")
            self.assertEqual(len(record.table_rows[13]), 4)
            self.assertEqual(record.table_rows[13][0], "141212_at")
            self.assertEqual(record.table_rows[13][1], "1200")
            self.assertEqual(record.table_rows[13][2], "P")
            self.assertEqual(record.table_rows[13][3], "0.000491")
            self.assertEqual(len(record.table_rows[14]), 4)
            self.assertEqual(record.table_rows[14][0], "141213_at")
            self.assertEqual(record.table_rows[14][1], "13.7")
            self.assertEqual(record.table_rows[14][2], "A")
            self.assertEqual(record.table_rows[14][3], "0.635055")
            self.assertEqual(len(record.table_rows[15]), 4)
            self.assertEqual(record.table_rows[15][0], "142121_at")
            self.assertEqual(record.table_rows[15][1], "133.7")
            self.assertEqual(record.table_rows[15][2], "A")
            self.assertEqual(record.table_rows[15][3], "0.889551")
            self.assertEqual(len(record.table_rows[16]), 4)
            self.assertEqual(record.table_rows[16][0], "142122_at")
            self.assertEqual(record.table_rows[16][1], "275.3")
            self.assertEqual(record.table_rows[16][2], "A")
            self.assertEqual(record.table_rows[16][3], "0.611218")
            self.assertEqual(len(record.table_rows[17]), 4)
            self.assertEqual(record.table_rows[17][0], "142123_at")
            self.assertEqual(record.table_rows[17][1], "307.6")
            self.assertEqual(record.table_rows[17][2], "A")
            self.assertEqual(record.table_rows[17][3], "0.611218")
            self.assertEqual(len(record.table_rows[18]), 4)
            self.assertEqual(record.table_rows[18][0], "142124_at")
            self.assertEqual(record.table_rows[18][1], "132.6")
            self.assertEqual(record.table_rows[18][2], "A")
            self.assertEqual(record.table_rows[18][3], "0.437646")
            self.assertEqual(len(record.table_rows[19]), 4)
            self.assertEqual(record.table_rows[19][0], "142125_at")
            self.assertEqual(record.table_rows[19][1], "195.8")
            self.assertEqual(record.table_rows[19][2], "A")
            self.assertEqual(record.table_rows[19][3], "0.110449")
            self.assertEqual(len(record.table_rows[20]), 4)
            self.assertEqual(record.table_rows[20][0], "142126_at")
            self.assertEqual(record.table_rows[20][1], "174.1")
            self.assertEqual(record.table_rows[20][2], "A")
            self.assertEqual(record.table_rows[20][3], "0.681117")
            self.assertEqual(len(record.table_rows[21]), 4)
            self.assertEqual(record.table_rows[21][0], "142127_at")
            self.assertEqual(record.table_rows[21][1], "316.3")
            self.assertEqual(record.table_rows[21][2], "A")
            self.assertEqual(record.table_rows[21][3], "0.65838")
            record = next(records)
            self.assertEqual(record.entity_type, "SERIES")
            self.assertEqual(record.entity_id, "Dros_embryo_timecourse")
            self.assertEqual(len(record.entity_attributes), 6)
            self.assertEqual(len(record.entity_attributes["Series_sample_id"]), 3)
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][0], "Drosophila_T0-1"
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][1], "Drosophila_T0-2"
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][2], "Drosophila_T1-1"
            )
            self.assertEqual(len(record.entity_attributes["Series_contributor"]), 5)
            self.assertEqual(
                record.entity_attributes["Series_contributor"][0], "Jane,Doe"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][1], "John,A,Smith"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][2], "Hans,van Elton"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][3], "John,Smithers Jr"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][4], "Jie,D,Chen"
            )
            self.assertEqual(len(record.entity_attributes["Series_summary"]), 2)
            self.assertEqual(
                record.entity_attributes["Series_summary"][0],
                "Morphogenesis of epithelial tissues relies on the precise developmental control of cell polarity and architecture. In the early Drosophila embryo, the primary epithelium forms during cellularisation, following a tightly controlled genetic programme where specific sets of genes are up-regulated. Some of them, for instance, control membrane invagination between the nuclei anchored at the apical surface of the syncytium.",
            )
            self.assertEqual(
                record.entity_attributes["Series_summary"][1],
                "We used microarrays to detail the global programme of gene expression underlying cellularisation and identified distinct classes of up-regulated genes during this process.",
            )
            self.assertEqual(record.entity_attributes["Series_type"], "time course")
            self.assertEqual(
                record.entity_attributes["Series_title"],
                "Expression data from early Drosophila embryo",
            )
            self.assertEqual(
                record.entity_attributes["Series_overall_design"],
                "Drosophila embryos were selected at successive stages of early development for RNA extraction and hybridization on Affymetrix microarrays. We sought to obtain homogeneous populations of embryos at each developmental stage in order to increase the temporal resolution of expression profiles. To that end, we hand-selected embryos according to morphological criteria at five time-points: before pole cell formation, i.e. before zygotic transcription (T0), during the slow phase (T1) and the fast phase (T2) of cellularisation and at the beginning (T3) and the end (T4) of gastrulation.",
            )
            self.assertEqual(len(record.col_defs), 0)
            self.assertEqual(len(record.table_rows), 0)

    def test_GSE16(self):
        path = "Geo/GSE16.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "GSM804")
            self.assertEqual(len(record.entity_attributes), 18)
            self.assertEqual(record.entity_attributes["Sample_pubmed_id"], "11687795")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_institute"],
                "University of California San Francisco",
            )
            self.assertEqual(len(record.entity_attributes["Sample_author"]), 19)
            self.assertEqual(
                record.entity_attributes["Sample_author"][0], "Antoine,M,Snijders"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][1], "Norma,,Nowak"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][2], "Richard,,Segraves"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][3], "Stephanie,,Blackwood"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][4], "Nils,,Brown"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][5], "Jeffery,,Conroy"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][6], "Greg,,Hamilton"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][7], "Anna,K,Hindle"
            )
            self.assertEqual(record.entity_attributes["Sample_author"][8], "Bing,,Huey")
            self.assertEqual(
                record.entity_attributes["Sample_author"][9], "Karen,,Kimura"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][10], "Sindy,,Law"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][11], "Ken,,Myambo"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][12], "Joel,,Palmer"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][13], "Bauke,,Ylstra"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][14], "Jingzhu,P,Yue"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][15], "Joe,W,Gray"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][16], "Ajay,N,Jain"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][17], "Daniel,,Pinkel"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][18], "Donna,G,Albertson"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_phone"], "415 502-8463"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_department"],
                "Comprehensive Cancer Center",
            )
            self.assertEqual(len(record.entity_attributes["Sample_description"]), 4)
            self.assertEqual(
                record.entity_attributes["Sample_description"][0],
                'Coriell Cell Repositories cell line <a href="http://locus.umdnj.edu/nigms/nigms_cgi/display.cgi?GM05296">GM05296</a>.',
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][1],
                "Fibroblast cell line derived from a 1 month old female with multiple congenital malformations, dysmorphic features, intrauterine growth retardation, heart murmur, cleft palate, equinovarus deformity, microcephaly, coloboma of right iris, clinodactyly, reduced RBC catalase activity, and 1 copy of catalase gene.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][2],
                "Chromosome abnormalities are present.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][3],
                "Karyotype is 46,XX,-11,+der(11)inv ins(11;10)(11pter> 11p13::10q21>10q24::11p13>11qter)mat",
            )
            self.assertEqual(
                record.entity_attributes["Sample_target_source2"],
                "normal male reference genomic DNA",
            )
            self.assertEqual(
                record.entity_attributes["Sample_target_source1"], "Cell line GM05296"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_name"], "Donna,G,Albertson"
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL28")
            self.assertEqual(
                record.entity_attributes["Sample_type"], "dual channel genomic"
            )
            self.assertEqual(
                record.entity_attributes["Sample_status"], "Public on Feb 12 2002"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_email"],
                "albertson@cc.ucsf.edu",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"], "CGH_Albertson_GM05296-001218"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism"], "Homo sapiens"
            )
            self.assertEqual(record.entity_attributes["Sample_series_id"], "GSE16")
            self.assertEqual(
                record.entity_attributes["Sample_submission_date"], "Jan 17 2002"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_city"],
                "San Francisco,CA,94143,USA",
            )
            self.assertEqual(len(record.col_defs), 5)
            self.assertEqual(
                record.col_defs["NO_REPLICATES"],
                "Number of replicate spot measurements",
            )
            self.assertEqual(
                record.col_defs["LOG2STDDEV"], "Standard deviation of VALUE"
            )
            self.assertEqual(
                record.col_defs["ID_REF"],
                "Unique row identifier, genome position order",
            )
            self.assertEqual(
                record.col_defs["VALUE"],
                "aka LOG2RATIO, mean of log base 2 of LINEAR_RATIO",
            )
            self.assertEqual(
                record.col_defs["LINEAR_RATIO"], "Mean of replicate Cy3/Cy5 ratios"
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 5)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LINEAR_RATIO")
            self.assertEqual(record.table_rows[0][3], "LOG2STDDEV")
            self.assertEqual(record.table_rows[0][4], "NO_REPLICATES")
            self.assertEqual(len(record.table_rows[1]), 5)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "")
            self.assertEqual(record.table_rows[1][2], "1.047765")
            self.assertEqual(record.table_rows[1][3], "0.011853")
            self.assertEqual(record.table_rows[1][4], "3")
            self.assertEqual(len(record.table_rows[2]), 5)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "")
            self.assertEqual(record.table_rows[2][2], "")
            self.assertEqual(record.table_rows[2][3], "")
            self.assertEqual(record.table_rows[2][4], "0")
            self.assertEqual(len(record.table_rows[3]), 5)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.008824")
            self.assertEqual(record.table_rows[3][2], "1.006135")
            self.assertEqual(record.table_rows[3][3], "0.00143")
            self.assertEqual(record.table_rows[3][4], "3")
            self.assertEqual(len(record.table_rows[4]), 5)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.000894")
            self.assertEqual(record.table_rows[4][2], "0.99938")
            self.assertEqual(record.table_rows[4][3], "0.001454")
            self.assertEqual(record.table_rows[4][4], "3")
            self.assertEqual(len(record.table_rows[5]), 5)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "0.075875")
            self.assertEqual(record.table_rows[5][2], "1.054")
            self.assertEqual(record.table_rows[5][3], "0.003077")
            self.assertEqual(record.table_rows[5][4], "3")
            self.assertEqual(len(record.table_rows[6]), 5)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "0.017303")
            self.assertEqual(record.table_rows[6][2], "1.012066")
            self.assertEqual(record.table_rows[6][3], "0.005876")
            self.assertEqual(record.table_rows[6][4], "2")
            self.assertEqual(len(record.table_rows[7]), 5)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-0.006766")
            self.assertEqual(record.table_rows[7][2], "0.995321")
            self.assertEqual(record.table_rows[7][3], "0.013881")
            self.assertEqual(record.table_rows[7][4], "3")
            self.assertEqual(len(record.table_rows[8]), 5)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "0.020755")
            self.assertEqual(record.table_rows[8][2], "1.014491")
            self.assertEqual(record.table_rows[8][3], "0.005506")
            self.assertEqual(record.table_rows[8][4], "3")
            self.assertEqual(len(record.table_rows[9]), 5)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "-0.094938")
            self.assertEqual(record.table_rows[9][2], "0.936313")
            self.assertEqual(record.table_rows[9][3], "0.012662")
            self.assertEqual(record.table_rows[9][4], "3")
            self.assertEqual(len(record.table_rows[10]), 5)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "-0.054527")
            self.assertEqual(record.table_rows[10][2], "0.96291")
            self.assertEqual(record.table_rows[10][3], "0.01073")
            self.assertEqual(record.table_rows[10][4], "3")
            self.assertEqual(len(record.table_rows[11]), 5)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "-0.025057")
            self.assertEqual(record.table_rows[11][2], "0.982782")
            self.assertEqual(record.table_rows[11][3], "0.003855")
            self.assertEqual(record.table_rows[11][4], "3")
            self.assertEqual(len(record.table_rows[12]), 5)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "")
            self.assertEqual(record.table_rows[12][2], "")
            self.assertEqual(record.table_rows[12][3], "")
            self.assertEqual(record.table_rows[12][4], "0")
            self.assertEqual(len(record.table_rows[13]), 5)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "0.108454")
            self.assertEqual(record.table_rows[13][2], "1.078072")
            self.assertEqual(record.table_rows[13][3], "0.005196")
            self.assertEqual(record.table_rows[13][4], "3")
            self.assertEqual(len(record.table_rows[14]), 5)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "0.078633")
            self.assertEqual(record.table_rows[14][2], "1.056017")
            self.assertEqual(record.table_rows[14][3], "0.009165")
            self.assertEqual(record.table_rows[14][4], "3")
            self.assertEqual(len(record.table_rows[15]), 5)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "0.098571")
            self.assertEqual(record.table_rows[15][2], "1.070712")
            self.assertEqual(record.table_rows[15][3], "0.007834")
            self.assertEqual(record.table_rows[15][4], "3")
            self.assertEqual(len(record.table_rows[16]), 5)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "0.044048")
            self.assertEqual(record.table_rows[16][2], "1.031003")
            self.assertEqual(record.table_rows[16][3], "0.013651")
            self.assertEqual(record.table_rows[16][4], "3")
            self.assertEqual(len(record.table_rows[17]), 5)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "0.018039")
            self.assertEqual(record.table_rows[17][2], "1.012582")
            self.assertEqual(record.table_rows[17][3], "0.005471")
            self.assertEqual(record.table_rows[17][4], "3")
            self.assertEqual(len(record.table_rows[18]), 5)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.088807")
            self.assertEqual(record.table_rows[18][2], "0.9403")
            self.assertEqual(record.table_rows[18][3], "0.010571")
            self.assertEqual(record.table_rows[18][4], "3")
            self.assertEqual(len(record.table_rows[19]), 5)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.016349")
            self.assertEqual(record.table_rows[19][2], "1.011397")
            self.assertEqual(record.table_rows[19][3], "0.007113")
            self.assertEqual(record.table_rows[19][4], "3")
            self.assertEqual(len(record.table_rows[20]), 5)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "0.030977")
            self.assertEqual(record.table_rows[20][2], "1.021704")
            self.assertEqual(record.table_rows[20][3], "0.016798")
            self.assertEqual(record.table_rows[20][4], "3")

    def test_soft_ex_platform(self):
        path = "Geo/soft_ex_platform.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "PLATFORM")
            self.assertEqual(
                record.entity_id, "Murine 15K long oligo array version 2.0"
            )
            self.assertEqual(len(record.entity_attributes), 12)
            self.assertEqual(
                record.entity_attributes["Platform_title"],
                "Murine 15K long oligo array version 2.0",
            )
            self.assertEqual(
                record.entity_attributes["Platform_web_link"],
                "http://www.microarray.protocols.html",
            )
            self.assertEqual(record.entity_attributes["platform_table_end"], "")
            self.assertEqual(record.entity_attributes["Platform_support"], "glass")
            self.assertEqual(
                record.entity_attributes["Platform_manufacturer"],
                "Un. London microarray facility",
            )
            self.assertEqual(record.entity_attributes["Platform_coating"], "polysine")
            self.assertEqual(
                record.entity_attributes["Platform_technology"],
                "spotted oligonucleotide",
            )
            self.assertEqual(record.entity_attributes["platform_table_begin"], "")
            self.assertEqual(
                len(record.entity_attributes["Platform_manufacture_protocol"]), 12
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][0],
                "1.  Oligos are arrayed in Greiner 384-well flat-bottom plates. Each well contains 600 pmol of 70-mer oligo.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][1],
                "2. Resuspend oligos in water to 20 uM and rearray 5 \xb5L into 384-well, Genetix polystyrene V-bottom plates (cat# X6004).",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][2],
                "3. Allow Genetix plates to dry through passive water evaporation in a protected environment (e.g., chemical hood).",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][3],
                "4. Before printing, add 5 \xb5L of 1X Printing Buffer to each well. This can be done the night before a print run is started.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][4],
                "5. Seal plates with Corning seals.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][5],
                "6. Incubate at 37\xb0C for 30 minutes to aid resuspension of DNA.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][6],
                "7. Shake plates near maximum rotational speed on flat-bed shaker for 1 minute.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][7],
                "8. Centrifuge plates at 2000 rpm for 3 minutes.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][8],
                "9. Remove seals and cover with plate lids. Place in appropriate location of plate cassette. This should be done with first plates just before print run is started to minimize evaporation time before printing. For second and third cassettes, wait until 30 minutes before next cassette is needed to begin centrifugation.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][9],
                "10. Make sure plates rest behind both holding clips in the cassettes. Push plates back into the cassettes as far as they will go, putting them in the proper position for the server arm.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][10],
                "11. After the print run is completed, allow plates to dry through passive evaporation in a protected environment.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][11],
                "12. For each subsequent preparation of these plates for a print run, add water to the wells instead of sodium phosphate buffer. The amount of water should be decreased by 0.25 \xb5L per print run, as this is the amount drawn up by the pin capillary during each dip.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_organism"], "Mus musculus"
            )
            self.assertEqual(len(record.entity_attributes["Platform_contributor"]), 5)
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][0], "Jane,Doe"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][1], "John,A,Smith"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][2], "Hans,van Elton"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][3], "John,Smithers Jr"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][4], "Jie,D,Chen"
            )
            self.assertEqual(
                record.entity_attributes["Platform_distribution"], "non-commercial"
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["Gene_Desc"], "Gene description")
            self.assertEqual(record.col_defs["SEQUENCE"], "Probe sequence information")
            self.assertEqual(record.col_defs["Gene_Sym"], "Gene symbols")
            self.assertEqual(
                record.col_defs["GB_ACC"],
                'GenBank accession number of sequence used to design oligonucleotide probe   LINK_PRE:"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&term="',
            )
            self.assertEqual(record.col_defs["SPOT_ID"], "alternative identifier")
            self.assertEqual(record.col_defs["ID"], "")
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID")
            self.assertEqual(record.table_rows[0][1], "GB_ACC")
            self.assertEqual(record.table_rows[0][2], "Gene_Desc")
            self.assertEqual(record.table_rows[0][3], "Gene_Sym")
            self.assertEqual(record.table_rows[0][4], "SPOT_ID")
            self.assertEqual(record.table_rows[0][5], "SEQUENCE")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "U02079")
            self.assertEqual(
                record.table_rows[1][2],
                "nuclear factor of activated T-cells, cytoplasmic 2",
            )
            self.assertEqual(record.table_rows[1][3], "Nfatc2")
            self.assertEqual(record.table_rows[1][4], "")
            self.assertEqual(
                record.table_rows[1][5],
                "ACCTGGATGACGCAGCCACTTCAGAAAGCTGGGTTGGGACAGAAAGGTATATAGAGAGAAAATTTTGGAA",
            )
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "NM_008154")
            self.assertEqual(record.table_rows[2][2], "G-protein coupled receptor 3")
            self.assertEqual(record.table_rows[2][3], "Gpr3")
            self.assertEqual(record.table_rows[2][4], "")
            self.assertEqual(
                record.table_rows[2][5],
                "CTGTACAATGCTCTCACTTACTACTCAGAGACAACGGTAACTCGGACTTATGTGATGCTGGCCTTGGTGT",
            )
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "AK015719")
            self.assertEqual(record.table_rows[3][2], "tropomodulin 2")
            self.assertEqual(record.table_rows[3][3], "Tmod2")
            self.assertEqual(record.table_rows[3][4], "")
            self.assertEqual(
                record.table_rows[3][5],
                "CACCAGGCTCAGTGCCTAGTATCGGCTTCACCTAGTGTGGTTACTCAGGGCACGCAGAGCTACAGAACAC",
            )
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "AK003367")
            self.assertEqual(
                record.table_rows[4][2], "mitochondrial ribosomal protein L15"
            )
            self.assertEqual(record.table_rows[4][3], "Mrpl15")
            self.assertEqual(record.table_rows[4][4], "")
            self.assertEqual(
                record.table_rows[4][5],
                "CAAGAAGTCTAGAAATTCTGTGCAAGCCTATTCCATTCTTTCTGCGGGGACAACCAATTCCGAAAAGAAT",
            )
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "BC003333")
            self.assertEqual(record.table_rows[5][2], "RIKEN cDNA 0610033I05 gene")
            self.assertEqual(record.table_rows[5][3], "0610033I05Rik")
            self.assertEqual(record.table_rows[5][4], "")
            self.assertEqual(
                record.table_rows[5][5],
                "AGAACTGGGTGGCAGATATCCTAGAGTTTTGACCAACGTTCACAGCACACATATTGATCTTATAGGACCT",
            )
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "NM_008462")
            self.assertEqual(
                record.table_rows[6][2],
                "killer cell lectin-like receptor, subfamily A, member 2",
            )
            self.assertEqual(record.table_rows[6][3], "Klra2")
            self.assertEqual(record.table_rows[6][4], "")
            self.assertEqual(
                record.table_rows[6][5],
                "TGAATTGAAGTTCCTTAAATCCCAACTTCAAAGAAACACATACTGGATTTCACTGACACATCATAAAAGC",
            )
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "NM_008029")
            self.assertEqual(record.table_rows[7][2], "FMS-like tyrosine kinase 4")
            self.assertEqual(record.table_rows[7][3], "Flt4")
            self.assertEqual(record.table_rows[7][4], "")
            self.assertEqual(
                record.table_rows[7][5],
                "GAGGTGCTGTGGGATGACCGCCGGGGCATGCGGGTGCCCACTCAACTGTTGCGCGATGCCCTGTACCTGC",
            )
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "NM_054088")
            self.assertEqual(record.table_rows[8][2], "adiponutrin")
            self.assertEqual(record.table_rows[8][3], "Adpn")
            self.assertEqual(record.table_rows[8][4], "")
            self.assertEqual(
                record.table_rows[8][5],
                "GTCTGAGTTCCATTCCAAAGACGAAGTCGTGGATGCCCTGGTGTGTTCCTGCTTCATTCCCCTCTTCTCT",
            )
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "NM_009750")
            self.assertEqual(
                record.table_rows[9][2],
                "nerve growth factor receptor (TNFRSF16) associated protein 1",
            )
            self.assertEqual(record.table_rows[9][3], "Ngfrap1")
            self.assertEqual(record.table_rows[9][4], "")
            self.assertEqual(
                record.table_rows[9][5],
                "TACAGCTGAGAAATTGTCTACGCATCCTTATGGGGGAGCTGTCTAACCACCACGATCACCATGATGAATT",
            )
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "AB045323")
            self.assertEqual(
                record.table_rows[10][2], "DNA segment, Chr 8, ERATO Doi 594, expressed"
            )
            self.assertEqual(record.table_rows[10][3], "D8Ertd594e")
            self.assertEqual(record.table_rows[10][4], "")
            self.assertEqual(
                record.table_rows[10][5],
                "GATTCAGACTCGGGAGGAGCATCCCAACCTCTCCTTGAGGATAAAGGCCTGAGCGATTGCCCTGGGGAGC",
            )
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "AK005789")
            self.assertEqual(
                record.table_rows[11][2], "dynein, cytoplasmic, light chain 2B"
            )
            self.assertEqual(record.table_rows[11][3], "Dncl2b")
            self.assertEqual(record.table_rows[11][4], "")
            self.assertEqual(
                record.table_rows[11][5],
                "TGCAGAAGGCATTCCAATCCGAACAACCCTGGACAACTCCACAACGGTTCAGTATGCGGGTCTTCTCCAC",
            )
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "NM_010517")
            self.assertEqual(
                record.table_rows[12][2], "insulin-like growth factor binding protein 4"
            )
            self.assertEqual(record.table_rows[12][3], "Igfbp4")
            self.assertEqual(record.table_rows[12][4], "")
            self.assertEqual(
                record.table_rows[12][5],
                "GGAGAAGCTGGCGCGCTGCCGCCCCCCCGTGGGTTGCGAGGAGTTGGTGCGGGAGCCAGGCTGCGGTTGT",
            )
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "AK010722")
            self.assertEqual(record.table_rows[13][2], "RIKEN cDNA 2410075D05 gene")
            self.assertEqual(record.table_rows[13][3], "2410075D05Rik")
            self.assertEqual(record.table_rows[13][4], "")
            self.assertEqual(
                record.table_rows[13][5],
                "GGAGCATCTGGAGTTCCGCTTACCGGAAATAAAGTCTTTACTATCGGTGATTGGAGGGCAGTTCACTAAC",
            )
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "AK003755")
            self.assertEqual(
                record.table_rows[14][2], "DNA segment, Chr 4, ERATO Doi 421, expressed"
            )
            self.assertEqual(record.table_rows[14][3], "D4Ertd421e")
            self.assertEqual(record.table_rows[14][4], "")
            self.assertEqual(
                record.table_rows[14][5],
                "AGCAAAGAGATCTCCCTCAGTGTGCCCATAGGTGGCGGTGCGAGCTTGCGGTTATTGGCCAGTGACTTGC",
            )
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "BC003241")
            self.assertEqual(
                record.table_rows[15][2],
                "cleavage stimulation factor, 3' pre-RNA, subunit 3",
            )
            self.assertEqual(record.table_rows[15][3], "Cstf3")
            self.assertEqual(record.table_rows[15][4], "")
            self.assertEqual(
                record.table_rows[15][5],
                "AAATTAGAAGAAAATCCATATGACCTTGATGCTTGGAGCATTCTCATTCGAGAGGCACAGAATCAACCTA",
            )
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "AK004937")
            self.assertEqual(record.table_rows[16][2], "RIKEN cDNA 1300007O09 gene")
            self.assertEqual(record.table_rows[16][3], "1300007O09Rik")
            self.assertEqual(record.table_rows[16][4], "")
            self.assertEqual(
                record.table_rows[16][5],
                "CAGACACAAACCCTAGGTTGTATTGTAGACCGGAGTTTAAGCAGGCACTACCTGTCTGTCTTTTCTTCAT",
            )
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "AK004524")
            self.assertEqual(
                record.table_rows[17][2],
                "unnamed protein product; hypothetical SOCS domain",
            )
            self.assertEqual(record.table_rows[17][3], "")
            self.assertEqual(record.table_rows[17][4], "")
            self.assertEqual(
                record.table_rows[17][5],
                "CGGAGCCCTGCGCGCCCAGAGCCCCCTCCCACCCGCTTCCACCAAGTGCATGGAGCCAACATCCGCATGG",
            )
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "NM_025999")
            self.assertEqual(record.table_rows[18][2], "RIKEN cDNA 2610110L04 gene")
            self.assertEqual(record.table_rows[18][3], "2610110L04Rik")
            self.assertEqual(record.table_rows[18][4], "")
            self.assertEqual(
                record.table_rows[18][5],
                "TGCATTGATAAATGGAGTGATCGACACAGGAACTGCCCCATTTGTCGCCTACAGATGACTGGAGCAAATG",
            )
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "")
            self.assertEqual(record.table_rows[19][2], "")
            self.assertEqual(record.table_rows[19][3], "")
            self.assertEqual(record.table_rows[19][4], "-- CONTROL")
            self.assertEqual(record.table_rows[19][5], "")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "NM_023120")
            self.assertEqual(
                record.table_rows[20][2],
                "guanine nucleotide binding protein (G protein), beta polypeptide 1-like",
            )
            self.assertEqual(record.table_rows[20][3], "Gnb1l")
            self.assertEqual(record.table_rows[20][4], "")
            self.assertEqual(
                record.table_rows[20][5],
                "ACCGCCTGGTCCCAGATTTGTCCTCCGAGGCACACAGTCGGCTGTGAACACGCTCCATTTCTGCCCACCA",
            )

    def test_GSM700(self):
        path = "Geo/GSM700.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "GSM700")
            self.assertEqual(len(record.entity_attributes), 20)
            self.assertEqual(
                record.entity_attributes["Sample_submitter_institute"],
                "National Cancer Institute",
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_department"],
                "Cancer Genome Anatomy Project",
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_web_link"],
                "http://cgap.nci.nih.gov/",
            )
            self.assertEqual(len(record.entity_attributes["Sample_description"]), 14)
            self.assertEqual(
                record.entity_attributes["Sample_description"][0],
                "This library represents a Cancer Genome Anatomy Project library, which was either produced through CGAP funding, or donated to CGAP.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][1],
                "The Cancer Genome Anatomy Project (CGAP: http://cgap.nci.nih.gov) is an interdisciplinary program established and administered by the National Cancer Institute (NCI: http://www.nci.nih.gov) to generate the information and technological tools needed to decipher the molecular anatomy of the cancer cell.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][2],
                "Cell line grown under 1.5% oxygen conditions for 24 hours prior to harvesting in zinc option media with 10% RBS and harvested at passage 102. Library constructed in the laboratory of G. Riggins, M.D., Ph.D. (Duke University).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][3], "Organ: brain"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][4],
                "Tissue_type: glioblastoma multiforme",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][5], "Cell_line: H247"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][6], "Lab host: DH10B"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][7], "Vector: pZErO-1"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][8],
                "Vector type: plasmid",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][9], "R. Site 1: Sph1"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][10], "R. Site 2: Sph1"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][11],
                "Library treatment: non-normalized",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][12],
                "Tissue description: Brain, Duke glioblastoma multiforme cell line, H247, grown under 1.5% oxygen conditions  for 24 hours prior to harvesting.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][13], "Tissue"
            )
            self.assertEqual(len(record.entity_attributes["Sample_author"]), 2)
            self.assertEqual(
                record.entity_attributes["Sample_author"][0], "Gregory,J,Riggins"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][1], "Robert,L,Strausberg"
            )
            self.assertEqual(
                record.entity_attributes["Sample_web_link"], "http://cgap.nci.nih.gov"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_phone"], "301-496-1550"
            )
            self.assertEqual(record.entity_attributes["Sample_series_id"], "GSE14")
            self.assertEqual(record.entity_attributes["Sample_tag_count"], "72031")
            self.assertEqual(record.entity_attributes["Sample_type"], "sage")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_name"], "Robert,L,Strausberg"
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL4")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_city"],
                "Bethesda,MD,20892,USA",
            )
            self.assertEqual(
                record.entity_attributes["Sample_status"], "Public on Nov 28 2001"
            )
            self.assertEqual(record.entity_attributes["Sample_anchor"], "NlaIII")
            self.assertEqual(
                record.entity_attributes["Sample_title"], "SAGE_Duke_H247_Hypoxia"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism"], "Homo sapiens"
            )
            self.assertEqual(
                record.entity_attributes["Sample_target_source"],
                "Brain, glioblastoma multiforme, cell-line H247",
            )
            self.assertEqual(
                record.entity_attributes["Sample_submission_date"], "Nov 28 2001"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_email"],
                "cgapbs-r@mail.nih.gov",
            )
            self.assertEqual(len(record.col_defs), 3)
            self.assertEqual(record.col_defs["COUNT"], "Absolute tag count")
            self.assertEqual(
                record.col_defs["TPM"],
                "Tags per million, or (1000000*COUNT)/(Total tags)",
            )
            self.assertEqual(
                record.col_defs["TAG"],
                'Ten base SAGE tag, LINK_PRE:"http://www.ncbi.nlm.nih.gov/SAGE/SAGEtag.cgi?tag="',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 3)
            self.assertEqual(record.table_rows[0][0], "TAG")
            self.assertEqual(record.table_rows[0][1], "COUNT")
            self.assertEqual(record.table_rows[0][2], "TPM")
            self.assertEqual(len(record.table_rows[1]), 3)
            self.assertEqual(record.table_rows[1][0], "TCCAAATCGA")
            self.assertEqual(record.table_rows[1][1], "520")
            self.assertEqual(record.table_rows[1][2], "7219.11")
            self.assertEqual(len(record.table_rows[2]), 3)
            self.assertEqual(record.table_rows[2][0], "TACCATCAAT")
            self.assertEqual(record.table_rows[2][1], "434")
            self.assertEqual(record.table_rows[2][2], "6025.18")
            self.assertEqual(len(record.table_rows[3]), 3)
            self.assertEqual(record.table_rows[3][0], "TTGGGGTTTC")
            self.assertEqual(record.table_rows[3][1], "389")
            self.assertEqual(record.table_rows[3][2], "5400.45")
            self.assertEqual(len(record.table_rows[4]), 3)
            self.assertEqual(record.table_rows[4][0], "CCCATCGTCC")
            self.assertEqual(record.table_rows[4][1], "367")
            self.assertEqual(record.table_rows[4][2], "5095.03")
            self.assertEqual(len(record.table_rows[5]), 3)
            self.assertEqual(record.table_rows[5][0], "GTGAAACCCC")
            self.assertEqual(record.table_rows[5][1], "365")
            self.assertEqual(record.table_rows[5][2], "5067.26")
            self.assertEqual(len(record.table_rows[6]), 3)
            self.assertEqual(record.table_rows[6][0], "GGGGAAATCG")
            self.assertEqual(record.table_rows[6][1], "357")
            self.assertEqual(record.table_rows[6][2], "4956.2")
            self.assertEqual(len(record.table_rows[7]), 3)
            self.assertEqual(record.table_rows[7][0], "CCTGTAATCC")
            self.assertEqual(record.table_rows[7][1], "346")
            self.assertEqual(record.table_rows[7][2], "4803.49")
            self.assertEqual(len(record.table_rows[8]), 3)
            self.assertEqual(record.table_rows[8][0], "TGATTTCACT")
            self.assertEqual(record.table_rows[8][1], "334")
            self.assertEqual(record.table_rows[8][2], "4636.89")
            self.assertEqual(len(record.table_rows[9]), 3)
            self.assertEqual(record.table_rows[9][0], "TGTGTTGAGA")
            self.assertEqual(record.table_rows[9][1], "315")
            self.assertEqual(record.table_rows[9][2], "4373.12")
            self.assertEqual(len(record.table_rows[10]), 3)
            self.assertEqual(record.table_rows[10][0], "GCCCCCAATA")
            self.assertEqual(record.table_rows[10][1], "303")
            self.assertEqual(record.table_rows[10][2], "4206.52")
            self.assertEqual(len(record.table_rows[11]), 3)
            self.assertEqual(record.table_rows[11][0], "CTAAGACTTC")
            self.assertEqual(record.table_rows[11][1], "279")
            self.assertEqual(record.table_rows[11][2], "3873.33")
            self.assertEqual(len(record.table_rows[12]), 3)
            self.assertEqual(record.table_rows[12][0], "GCGACCGTCA")
            self.assertEqual(record.table_rows[12][1], "276")
            self.assertEqual(record.table_rows[12][2], "3831.68")
            self.assertEqual(len(record.table_rows[13]), 3)
            self.assertEqual(record.table_rows[13][0], "TTGGTCCTCT")
            self.assertEqual(record.table_rows[13][1], "276")
            self.assertEqual(record.table_rows[13][2], "3831.68")
            self.assertEqual(len(record.table_rows[14]), 3)
            self.assertEqual(record.table_rows[14][0], "CCTAGCTGGA")
            self.assertEqual(record.table_rows[14][1], "268")
            self.assertEqual(record.table_rows[14][2], "3720.62")
            self.assertEqual(len(record.table_rows[15]), 3)
            self.assertEqual(record.table_rows[15][0], "GATGAGGAGA")
            self.assertEqual(record.table_rows[15][1], "251")
            self.assertEqual(record.table_rows[15][2], "3484.61")
            self.assertEqual(len(record.table_rows[16]), 3)
            self.assertEqual(record.table_rows[16][0], "ACTTTTTCAA")
            self.assertEqual(record.table_rows[16][1], "244")
            self.assertEqual(record.table_rows[16][2], "3387.43")
            self.assertEqual(len(record.table_rows[17]), 3)
            self.assertEqual(record.table_rows[17][0], "CCACTGCACT")
            self.assertEqual(record.table_rows[17][1], "223")
            self.assertEqual(record.table_rows[17][2], "3095.89")
            self.assertEqual(len(record.table_rows[18]), 3)
            self.assertEqual(record.table_rows[18][0], "GTGTGTTTGT")
            self.assertEqual(record.table_rows[18][1], "223")
            self.assertEqual(record.table_rows[18][2], "3095.89")
            self.assertEqual(len(record.table_rows[19]), 3)
            self.assertEqual(record.table_rows[19][0], "GAAATACAGT")
            self.assertEqual(record.table_rows[19][1], "218")
            self.assertEqual(record.table_rows[19][2], "3026.47")
            self.assertEqual(len(record.table_rows[20]), 3)
            self.assertEqual(record.table_rows[20][0], "GCTTTATTTG")
            self.assertEqual(record.table_rows[20][1], "218")
            self.assertEqual(record.table_rows[20][2], "3026.47")

    def test_GSM645(self):
        path = "Geo/GSM645.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "GSM645")
            self.assertEqual(len(record.entity_attributes), 17)
            self.assertEqual(
                record.entity_attributes["Sample_submitter_institute"],
                "Max von Pettenkofer Institut",
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_department"], "Bacteriology"
            )
            self.assertEqual(len(record.entity_attributes["Sample_author"]), 4)
            self.assertEqual(
                record.entity_attributes["Sample_author"][0], "Reinhard,,Hoffmann"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][1], "Thomas,,Seidl"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][2], "Ton,,Rolink"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][3], "Fritz,,Melchers"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_phone"], "+49-89-5160-5424"
            )
            self.assertEqual(record.entity_attributes["Sample_series_id"], "GSE13")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "B220+CD25+sIg- Large Pre BII cells sorted out of mouse bone marrow, sort no. 8",
            )
            self.assertEqual(record.entity_attributes["Sample_type"], "single channel")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_name"], "Reinhard,,Hoffmann"
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL22")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_city"],
                "Munich,80336,Germany",
            )
            self.assertEqual(
                record.entity_attributes["Sample_status"], "Public on Dec 17 2001"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_email"],
                "r_hoffmann@m3401.mpk.med.uni-muenchen.de",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"], "Large Pre-BII cells 8b"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism"], "Mus musculus"
            )
            self.assertEqual(
                record.entity_attributes["Sample_target_source"], "Large Pre-BII cells"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submission_date"], "Nov 27 2001"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_address"],
                "Pettenkoferstr. 9a",
            )
            self.assertEqual(len(record.col_defs), 14)
            self.assertEqual(
                record.col_defs["PAIRS"],
                "number of probe set specific probe pairs on the array",
            )
            self.assertEqual(
                record.col_defs["ABS_CALL"],
                "Whether a probe set is present, marginal, or absent; see Affymetrix Literature",
            )
            self.assertEqual(
                record.col_defs["PM Excess"],
                "number of probe pairs  where PM/MM exceeds the ratio limit (10 by default)",
            )
            self.assertEqual(
                record.col_defs["POSITIVE"], "number of poisitive probe pairs"
            )
            self.assertEqual(
                record.col_defs["MM Excess"],
                "Number of probe peirs where MM/PM exceeds 1/ratio limit (10 by default)",
            )
            self.assertEqual(
                record.col_defs["ID_REF"], "Affymetrix Probe Set Identifier"
            )
            self.assertEqual(
                record.col_defs["NEGATIVE"], "number of negative probe pairs"
            )
            self.assertEqual(record.col_defs["VALUE"], "Average Difference Intensity")
            self.assertEqual(record.col_defs["POS_FRACTION"], "Positive/Pairs Used")
            self.assertEqual(record.col_defs["Experiment Name"], "Experiment Name")
            self.assertEqual(record.col_defs["POS/NEG"], "Positive/Negative")
            self.assertEqual(record.col_defs["PAIRS_USED"], "")
            self.assertEqual(record.col_defs["Log Avg"], "")
            self.assertEqual(record.col_defs["PAIRS_IN_AVG"], "Trimmed probe pair set")
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 14)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "Experiment Name")
            self.assertEqual(record.table_rows[0][2], "POSITIVE")
            self.assertEqual(record.table_rows[0][3], "NEGATIVE")
            self.assertEqual(record.table_rows[0][4], "PAIRS")
            self.assertEqual(record.table_rows[0][5], "PAIRS_USED")
            self.assertEqual(record.table_rows[0][6], "PAIRS_IN_AVG")
            self.assertEqual(record.table_rows[0][7], "POS_FRACTION")
            self.assertEqual(record.table_rows[0][8], "Log Avg")
            self.assertEqual(record.table_rows[0][9], "PM Excess")
            self.assertEqual(record.table_rows[0][10], "MM Excess")
            self.assertEqual(record.table_rows[0][11], "POS/NEG")
            self.assertEqual(record.table_rows[0][12], "VALUE")
            self.assertEqual(record.table_rows[0][13], "ABS_CALL")
            self.assertEqual(len(record.table_rows[1]), 14)
            self.assertEqual(record.table_rows[1][0], "IL2_at")
            self.assertEqual(record.table_rows[1][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[1][2], "4")
            self.assertEqual(record.table_rows[1][3], "4")
            self.assertEqual(record.table_rows[1][4], "19")
            self.assertEqual(record.table_rows[1][5], "19")
            self.assertEqual(record.table_rows[1][6], "19")
            self.assertEqual(record.table_rows[1][7], "0.21")
            self.assertEqual(record.table_rows[1][8], "-0.58")
            self.assertEqual(record.table_rows[1][9], "0")
            self.assertEqual(record.table_rows[1][10], "0")
            self.assertEqual(record.table_rows[1][11], "1.0")
            self.assertEqual(record.table_rows[1][12], "-78")
            self.assertEqual(record.table_rows[1][13], "A")
            self.assertEqual(len(record.table_rows[2]), 14)
            self.assertEqual(record.table_rows[2][0], "IL10_at")
            self.assertEqual(record.table_rows[2][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[2][2], "7")
            self.assertEqual(record.table_rows[2][3], "4")
            self.assertEqual(record.table_rows[2][4], "20")
            self.assertEqual(record.table_rows[2][5], "20")
            self.assertEqual(record.table_rows[2][6], "18")
            self.assertEqual(record.table_rows[2][7], "0.35")
            self.assertEqual(record.table_rows[2][8], "1.87")
            self.assertEqual(record.table_rows[2][9], "1")
            self.assertEqual(record.table_rows[2][10], "0")
            self.assertEqual(record.table_rows[2][11], "1.8")
            self.assertEqual(record.table_rows[2][12], "161")
            self.assertEqual(record.table_rows[2][13], "A")
            self.assertEqual(len(record.table_rows[3]), 14)
            self.assertEqual(record.table_rows[3][0], "GMCSF_at")
            self.assertEqual(record.table_rows[3][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[3][2], "4")
            self.assertEqual(record.table_rows[3][3], "4")
            self.assertEqual(record.table_rows[3][4], "20")
            self.assertEqual(record.table_rows[3][5], "20")
            self.assertEqual(record.table_rows[3][6], "19")
            self.assertEqual(record.table_rows[3][7], "0.20")
            self.assertEqual(record.table_rows[3][8], "0.39")
            self.assertEqual(record.table_rows[3][9], "0")
            self.assertEqual(record.table_rows[3][10], "0")
            self.assertEqual(record.table_rows[3][11], "1.0")
            self.assertEqual(record.table_rows[3][12], "-11")
            self.assertEqual(record.table_rows[3][13], "A")
            self.assertEqual(len(record.table_rows[4]), 14)
            self.assertEqual(record.table_rows[4][0], "TNFRII_at")
            self.assertEqual(record.table_rows[4][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[4][2], "2")
            self.assertEqual(record.table_rows[4][3], "2")
            self.assertEqual(record.table_rows[4][4], "20")
            self.assertEqual(record.table_rows[4][5], "20")
            self.assertEqual(record.table_rows[4][6], "18")
            self.assertEqual(record.table_rows[4][7], "0.10")
            self.assertEqual(record.table_rows[4][8], "0.48")
            self.assertEqual(record.table_rows[4][9], "0")
            self.assertEqual(record.table_rows[4][10], "0")
            self.assertEqual(record.table_rows[4][11], "1.0")
            self.assertEqual(record.table_rows[4][12], "52")
            self.assertEqual(record.table_rows[4][13], "A")
            self.assertEqual(len(record.table_rows[5]), 14)
            self.assertEqual(record.table_rows[5][0], "MIP1-B_at")
            self.assertEqual(record.table_rows[5][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[5][2], "6")
            self.assertEqual(record.table_rows[5][3], "4")
            self.assertEqual(record.table_rows[5][4], "20")
            self.assertEqual(record.table_rows[5][5], "20")
            self.assertEqual(record.table_rows[5][6], "19")
            self.assertEqual(record.table_rows[5][7], "0.30")
            self.assertEqual(record.table_rows[5][8], "0.43")
            self.assertEqual(record.table_rows[5][9], "0")
            self.assertEqual(record.table_rows[5][10], "0")
            self.assertEqual(record.table_rows[5][11], "1.5")
            self.assertEqual(record.table_rows[5][12], "373")
            self.assertEqual(record.table_rows[5][13], "A")
            self.assertEqual(len(record.table_rows[6]), 14)
            self.assertEqual(record.table_rows[6][0], "IL4_at")
            self.assertEqual(record.table_rows[6][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[6][2], "3")
            self.assertEqual(record.table_rows[6][3], "3")
            self.assertEqual(record.table_rows[6][4], "20")
            self.assertEqual(record.table_rows[6][5], "20")
            self.assertEqual(record.table_rows[6][6], "19")
            self.assertEqual(record.table_rows[6][7], "0.15")
            self.assertEqual(record.table_rows[6][8], "0.29")
            self.assertEqual(record.table_rows[6][9], "0")
            self.assertEqual(record.table_rows[6][10], "0")
            self.assertEqual(record.table_rows[6][11], "1.0")
            self.assertEqual(record.table_rows[6][12], "27")
            self.assertEqual(record.table_rows[6][13], "A")
            self.assertEqual(len(record.table_rows[7]), 14)
            self.assertEqual(record.table_rows[7][0], "IL12_P40_at")
            self.assertEqual(record.table_rows[7][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[7][2], "3")
            self.assertEqual(record.table_rows[7][3], "5")
            self.assertEqual(record.table_rows[7][4], "20")
            self.assertEqual(record.table_rows[7][5], "20")
            self.assertEqual(record.table_rows[7][6], "19")
            self.assertEqual(record.table_rows[7][7], "0.15")
            self.assertEqual(record.table_rows[7][8], "-0.22")
            self.assertEqual(record.table_rows[7][9], "0")
            self.assertEqual(record.table_rows[7][10], "0")
            self.assertEqual(record.table_rows[7][11], "0.6")
            self.assertEqual(record.table_rows[7][12], "-163")
            self.assertEqual(record.table_rows[7][13], "A")
            self.assertEqual(len(record.table_rows[8]), 14)
            self.assertEqual(record.table_rows[8][0], "TNFa_at")
            self.assertEqual(record.table_rows[8][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[8][2], "3")
            self.assertEqual(record.table_rows[8][3], "4")
            self.assertEqual(record.table_rows[8][4], "20")
            self.assertEqual(record.table_rows[8][5], "20")
            self.assertEqual(record.table_rows[8][6], "20")
            self.assertEqual(record.table_rows[8][7], "0.15")
            self.assertEqual(record.table_rows[8][8], "-0.57")
            self.assertEqual(record.table_rows[8][9], "1")
            self.assertEqual(record.table_rows[8][10], "0")
            self.assertEqual(record.table_rows[8][11], "0.8")
            self.assertEqual(record.table_rows[8][12], "-95")
            self.assertEqual(record.table_rows[8][13], "A")
            self.assertEqual(len(record.table_rows[9]), 14)
            self.assertEqual(record.table_rows[9][0], "TCRa_at")
            self.assertEqual(record.table_rows[9][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[9][2], "1")
            self.assertEqual(record.table_rows[9][3], "4")
            self.assertEqual(record.table_rows[9][4], "20")
            self.assertEqual(record.table_rows[9][5], "20")
            self.assertEqual(record.table_rows[9][6], "19")
            self.assertEqual(record.table_rows[9][7], "0.05")
            self.assertEqual(record.table_rows[9][8], "-0.50")
            self.assertEqual(record.table_rows[9][9], "0")
            self.assertEqual(record.table_rows[9][10], "0")
            self.assertEqual(record.table_rows[9][11], "0.3")
            self.assertEqual(record.table_rows[9][12], "-186")
            self.assertEqual(record.table_rows[9][13], "A")
            self.assertEqual(len(record.table_rows[10]), 14)
            self.assertEqual(record.table_rows[10][0], "AFFX-BioB-5_at")
            self.assertEqual(record.table_rows[10][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[10][2], "0")
            self.assertEqual(record.table_rows[10][3], "1")
            self.assertEqual(record.table_rows[10][4], "20")
            self.assertEqual(record.table_rows[10][5], "20")
            self.assertEqual(record.table_rows[10][6], "19")
            self.assertEqual(record.table_rows[10][7], "0.00")
            self.assertEqual(record.table_rows[10][8], "0.35")
            self.assertEqual(record.table_rows[10][9], "0")
            self.assertEqual(record.table_rows[10][10], "0")
            self.assertEqual(record.table_rows[10][11], "0.0")
            self.assertEqual(record.table_rows[10][12], "120")
            self.assertEqual(record.table_rows[10][13], "A")
            self.assertEqual(len(record.table_rows[11]), 14)
            self.assertEqual(record.table_rows[11][0], "AFFX-BioB-M_at")
            self.assertEqual(record.table_rows[11][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[11][2], "0")
            self.assertEqual(record.table_rows[11][3], "1")
            self.assertEqual(record.table_rows[11][4], "20")
            self.assertEqual(record.table_rows[11][5], "20")
            self.assertEqual(record.table_rows[11][6], "19")
            self.assertEqual(record.table_rows[11][7], "0.00")
            self.assertEqual(record.table_rows[11][8], "0.02")
            self.assertEqual(record.table_rows[11][9], "0")
            self.assertEqual(record.table_rows[11][10], "0")
            self.assertEqual(record.table_rows[11][11], "0.0")
            self.assertEqual(record.table_rows[11][12], "-13")
            self.assertEqual(record.table_rows[11][13], "A")
            self.assertEqual(len(record.table_rows[12]), 14)
            self.assertEqual(record.table_rows[12][0], "AFFX-BioB-3_at")
            self.assertEqual(record.table_rows[12][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[12][2], "2")
            self.assertEqual(record.table_rows[12][3], "0")
            self.assertEqual(record.table_rows[12][4], "20")
            self.assertEqual(record.table_rows[12][5], "20")
            self.assertEqual(record.table_rows[12][6], "19")
            self.assertEqual(record.table_rows[12][7], "0.10")
            self.assertEqual(record.table_rows[12][8], "0.38")
            self.assertEqual(record.table_rows[12][9], "0")
            self.assertEqual(record.table_rows[12][10], "0")
            self.assertEqual(record.table_rows[12][11], "Undef")
            self.assertEqual(record.table_rows[12][12], "136")
            self.assertEqual(record.table_rows[12][13], "A")
            self.assertEqual(len(record.table_rows[13]), 14)
            self.assertEqual(record.table_rows[13][0], "AFFX-BioC-5_at")
            self.assertEqual(record.table_rows[13][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[13][2], "9")
            self.assertEqual(record.table_rows[13][3], "0")
            self.assertEqual(record.table_rows[13][4], "20")
            self.assertEqual(record.table_rows[13][5], "20")
            self.assertEqual(record.table_rows[13][6], "20")
            self.assertEqual(record.table_rows[13][7], "0.45")
            self.assertEqual(record.table_rows[13][8], "1.33")
            self.assertEqual(record.table_rows[13][9], "0")
            self.assertEqual(record.table_rows[13][10], "0")
            self.assertEqual(record.table_rows[13][11], "Undef")
            self.assertEqual(record.table_rows[13][12], "606")
            self.assertEqual(record.table_rows[13][13], "P")
            self.assertEqual(len(record.table_rows[14]), 14)
            self.assertEqual(record.table_rows[14][0], "AFFX-BioC-3_at")
            self.assertEqual(record.table_rows[14][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[14][2], "2")
            self.assertEqual(record.table_rows[14][3], "0")
            self.assertEqual(record.table_rows[14][4], "20")
            self.assertEqual(record.table_rows[14][5], "20")
            self.assertEqual(record.table_rows[14][6], "19")
            self.assertEqual(record.table_rows[14][7], "0.10")
            self.assertEqual(record.table_rows[14][8], "0.64")
            self.assertEqual(record.table_rows[14][9], "0")
            self.assertEqual(record.table_rows[14][10], "0")
            self.assertEqual(record.table_rows[14][11], "Undef")
            self.assertEqual(record.table_rows[14][12], "257")
            self.assertEqual(record.table_rows[14][13], "A")
            self.assertEqual(len(record.table_rows[15]), 14)
            self.assertEqual(record.table_rows[15][0], "AFFX-BioDn-5_at")
            self.assertEqual(record.table_rows[15][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[15][2], "8")
            self.assertEqual(record.table_rows[15][3], "0")
            self.assertEqual(record.table_rows[15][4], "20")
            self.assertEqual(record.table_rows[15][5], "20")
            self.assertEqual(record.table_rows[15][6], "20")
            self.assertEqual(record.table_rows[15][7], "0.40")
            self.assertEqual(record.table_rows[15][8], "1.23")
            self.assertEqual(record.table_rows[15][9], "0")
            self.assertEqual(record.table_rows[15][10], "0")
            self.assertEqual(record.table_rows[15][11], "Undef")
            self.assertEqual(record.table_rows[15][12], "380")
            self.assertEqual(record.table_rows[15][13], "P")
            self.assertEqual(len(record.table_rows[16]), 14)
            self.assertEqual(record.table_rows[16][0], "AFFX-BioDn-3_at")
            self.assertEqual(record.table_rows[16][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[16][2], "16")
            self.assertEqual(record.table_rows[16][3], "0")
            self.assertEqual(record.table_rows[16][4], "20")
            self.assertEqual(record.table_rows[16][5], "20")
            self.assertEqual(record.table_rows[16][6], "19")
            self.assertEqual(record.table_rows[16][7], "0.80")
            self.assertEqual(record.table_rows[16][8], "2.79")
            self.assertEqual(record.table_rows[16][9], "0")
            self.assertEqual(record.table_rows[16][10], "0")
            self.assertEqual(record.table_rows[16][11], "Undef")
            self.assertEqual(record.table_rows[16][12], "2764")
            self.assertEqual(record.table_rows[16][13], "P")
            self.assertEqual(len(record.table_rows[17]), 14)
            self.assertEqual(record.table_rows[17][0], "AFFX-CreX-5_at")
            self.assertEqual(record.table_rows[17][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[17][2], "19")
            self.assertEqual(record.table_rows[17][3], "0")
            self.assertEqual(record.table_rows[17][4], "20")
            self.assertEqual(record.table_rows[17][5], "20")
            self.assertEqual(record.table_rows[17][6], "19")
            self.assertEqual(record.table_rows[17][7], "0.95")
            self.assertEqual(record.table_rows[17][8], "5.65")
            self.assertEqual(record.table_rows[17][9], "0")
            self.assertEqual(record.table_rows[17][10], "0")
            self.assertEqual(record.table_rows[17][11], "Undef")
            self.assertEqual(record.table_rows[17][12], "4391")
            self.assertEqual(record.table_rows[17][13], "P")
            self.assertEqual(len(record.table_rows[18]), 14)
            self.assertEqual(record.table_rows[18][0], "AFFX-CreX-3_at")
            self.assertEqual(record.table_rows[18][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[18][2], "19")
            self.assertEqual(record.table_rows[18][3], "0")
            self.assertEqual(record.table_rows[18][4], "20")
            self.assertEqual(record.table_rows[18][5], "20")
            self.assertEqual(record.table_rows[18][6], "20")
            self.assertEqual(record.table_rows[18][7], "0.95")
            self.assertEqual(record.table_rows[18][8], "6.42")
            self.assertEqual(record.table_rows[18][9], "2")
            self.assertEqual(record.table_rows[18][10], "0")
            self.assertEqual(record.table_rows[18][11], "Undef")
            self.assertEqual(record.table_rows[18][12], "10787")
            self.assertEqual(record.table_rows[18][13], "P")
            self.assertEqual(len(record.table_rows[19]), 14)
            self.assertEqual(record.table_rows[19][0], "AFFX-BioB-5_st")
            self.assertEqual(record.table_rows[19][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[19][2], "5")
            self.assertEqual(record.table_rows[19][3], "3")
            self.assertEqual(record.table_rows[19][4], "20")
            self.assertEqual(record.table_rows[19][5], "20")
            self.assertEqual(record.table_rows[19][6], "19")
            self.assertEqual(record.table_rows[19][7], "0.25")
            self.assertEqual(record.table_rows[19][8], "0.48")
            self.assertEqual(record.table_rows[19][9], "0")
            self.assertEqual(record.table_rows[19][10], "0")
            self.assertEqual(record.table_rows[19][11], "1.7")
            self.assertEqual(record.table_rows[19][12], "80")
            self.assertEqual(record.table_rows[19][13], "A")
            self.assertEqual(len(record.table_rows[20]), 14)
            self.assertEqual(record.table_rows[20][0], "AFFX-BioB-M_st")
            self.assertEqual(record.table_rows[20][1], "RHMu8LarB")
            self.assertEqual(record.table_rows[20][2], "2")
            self.assertEqual(record.table_rows[20][3], "3")
            self.assertEqual(record.table_rows[20][4], "20")
            self.assertEqual(record.table_rows[20][5], "20")
            self.assertEqual(record.table_rows[20][6], "17")
            self.assertEqual(record.table_rows[20][7], "0.10")
            self.assertEqual(record.table_rows[20][8], "0.16")
            self.assertEqual(record.table_rows[20][9], "0")
            self.assertEqual(record.table_rows[20][10], "0")
            self.assertEqual(record.table_rows[20][11], "0.7")
            self.assertEqual(record.table_rows[20][12], "24")
            self.assertEqual(record.table_rows[20][13], "A")

    def test_soft_ex_series(self):
        path = "Geo/soft_ex_series.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SERIES")
            self.assertEqual(record.entity_id, "Bone_marrow_stromal_cells")
            self.assertEqual(len(record.entity_attributes), 13)
            self.assertEqual(
                record.entity_attributes["Series_variable_description_1"], "HS-5"
            )
            self.assertEqual(
                record.entity_attributes["Series_variable_sample_list_1"],
                "GSM10001, GSM10002",
            )
            self.assertEqual(len(record.entity_attributes["Series_sample_id"]), 4)
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][0], "GSM10001"
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][1], "GSM10002"
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][2], "GSM10003"
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][3], "GSM10004"
            )
            self.assertEqual(record.entity_attributes["Series_variable_2"], "cell line")
            self.assertEqual(record.entity_attributes["Series_pubmed_id"], "123456789")
            self.assertEqual(len(record.entity_attributes["Series_contributor"]), 5)
            self.assertEqual(
                record.entity_attributes["Series_contributor"][0], "Jane,Doe"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][1], "John,A,Smith"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][2], "Hans,van Elton"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][3], "John,Smithers Jr"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][4], "Jie,D,Chen"
            )
            self.assertEqual(
                record.entity_attributes["Series_summary"],
                "Two human stromal cell lines, HS-5 and HS-27a, represent functionally distinct components of the bone marrow microenvironment.1,2 HS-27a supports cobblestone area formation by early hematopoietic progenitors, whereas HS-5 secretes multiple cytokines that support the proliferation of committed progenitors. These cell lines have been distributed to research groups worldwide for use as a tool to understand interactions between hematopoietic cells and their microenvironment. We have used DNA microarray technology to characterize and compare the expression of over 17 000 genes in these cell lines. Gene expression differences in cytokines/chemokines, G-protein signaling molecules, and multiple extracellular matrix proteins add to the known protein and functional characterization of the lines, leading to new insight into the differences in their support function for hematopoietic progenitors.",
            )
            self.assertEqual(
                record.entity_attributes["Series_type"], "Cell Line Comparison"
            )
            self.assertEqual(record.entity_attributes["Series_variable_1"], "cell line")
            self.assertEqual(
                record.entity_attributes["Series_variable_description_2"], "HS-27a"
            )
            self.assertEqual(
                record.entity_attributes["Series_title"],
                "Profiling of the functionally distinct human bone marrow stromal cell lines HS-5 and HS-27a.",
            )
            self.assertEqual(
                record.entity_attributes["Series_variable_sample_list_2"],
                "GSM10003, GSM10004",
            )
            self.assertEqual(
                record.entity_attributes["Series_overall_design"],
                "We analyzed 2 arrays for HS-5 cell line and 2 arrays for HS-27a cell line",
            )
            self.assertEqual(len(record.col_defs), 0)
            self.assertEqual(len(record.table_rows), 0)

    def test_GSM691(self):
        path = "Geo/GSM691.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "GSM691")
            self.assertEqual(len(record.entity_attributes), 20)
            self.assertEqual(
                record.entity_attributes["Sample_submitter_institute"],
                "National Cancer Institute",
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_department"],
                "Cancer Genome Anatomy Project",
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_web_link"],
                "http://cgap.nci.nih.gov/",
            )
            self.assertEqual(len(record.entity_attributes["Sample_description"]), 12)
            self.assertEqual(
                record.entity_attributes["Sample_description"][0],
                "This library represents a Cancer Genome Anatomy Project library, which was either produced through CGAP funding, or donated to CGAP.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][1],
                "The Cancer Genome Anatomy Project (CGAP: http://cgap.nci.nih.gov) is an interdisciplinary program established and administered by the National Cancer Institute (NCI: http://www.nci.nih.gov) to generate the information and technological tools needed to decipher the molecular anatomy of the cancer cell.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][2],
                "Library constructed by Riggins laboratory Tissue supplied by Jeffrey Marks, Ph.D.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][3], "Organ: Breast"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][4],
                "Tissue_type: normal epithelial organoids",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][5],
                "Library treatment: non-normalized",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][6],
                "Tissue description: Breast, Isolated normal epithelial organoids. Derived from a reduction mammoplasty.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][7],
                "Tissue supplier: Jeffrey Marks, Ph.D.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][8], "Sample type: Bulk"
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][9],
                "Producer: Riggins Laboratory",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][10],
                "Clones generated to date: 768",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][11],
                "Sequences generated to date: 572",
            )
            self.assertEqual(len(record.entity_attributes["Sample_author"]), 3)
            self.assertEqual(
                record.entity_attributes["Sample_author"][0], "Jeffrey,,Marks"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][1], "Gregory,J,Riggins"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][2], "Robert,L,Strausberg"
            )
            self.assertEqual(
                record.entity_attributes["Sample_web_link"], "http://cgap.nci.nih.gov"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_phone"], "301-496-1550"
            )
            self.assertEqual(record.entity_attributes["Sample_series_id"], "GSE14")
            self.assertEqual(record.entity_attributes["Sample_tag_count"], "7165")
            self.assertEqual(record.entity_attributes["Sample_type"], "sage")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_name"], "Robert,L,Strausberg"
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL4")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_city"],
                "Bethesda,MD,20892,USA",
            )
            self.assertEqual(
                record.entity_attributes["Sample_status"], "Public on Nov 28 2001"
            )
            self.assertEqual(record.entity_attributes["Sample_anchor"], "NlaIII")
            self.assertEqual(record.entity_attributes["Sample_title"], "SAGE_Duke_40N")
            self.assertEqual(
                record.entity_attributes["Sample_organism"], "Homo sapiens"
            )
            self.assertEqual(
                record.entity_attributes["Sample_target_source"],
                "Breast, isolated normal epithelial organoids",
            )
            self.assertEqual(
                record.entity_attributes["Sample_submission_date"], "Nov 28 2001"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_email"],
                "cgapbs-r@mail.nih.gov",
            )
            self.assertEqual(len(record.col_defs), 3)
            self.assertEqual(record.col_defs["COUNT"], "Absolute tag count")
            self.assertEqual(
                record.col_defs["TPM"],
                "Tags per million, or (1000000*COUNT)/(Total tags)",
            )
            self.assertEqual(
                record.col_defs["TAG"],
                'Ten base SAGE tag, LINK_PRE:"http://www.ncbi.nlm.nih.gov/SAGE/SAGEtag.cgi?tag="',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 3)
            self.assertEqual(record.table_rows[0][0], "TAG")
            self.assertEqual(record.table_rows[0][1], "COUNT")
            self.assertEqual(record.table_rows[0][2], "TPM")
            self.assertEqual(len(record.table_rows[1]), 3)
            self.assertEqual(record.table_rows[1][0], "TTGGGGTTTC")
            self.assertEqual(record.table_rows[1][1], "202")
            self.assertEqual(record.table_rows[1][2], "28192.6")
            self.assertEqual(len(record.table_rows[2]), 3)
            self.assertEqual(record.table_rows[2][0], "TAGGTTGTCT")
            self.assertEqual(record.table_rows[2][1], "129")
            self.assertEqual(record.table_rows[2][2], "18004.2")
            self.assertEqual(len(record.table_rows[3]), 3)
            self.assertEqual(record.table_rows[3][0], "GAGGGAGTTT")
            self.assertEqual(record.table_rows[3][1], "109")
            self.assertEqual(record.table_rows[3][2], "15212.8")
            self.assertEqual(len(record.table_rows[4]), 3)
            self.assertEqual(record.table_rows[4][0], "TGCACGTTTT")
            self.assertEqual(record.table_rows[4][1], "92")
            self.assertEqual(record.table_rows[4][2], "12840.2")
            self.assertEqual(len(record.table_rows[5]), 3)
            self.assertEqual(record.table_rows[5][0], "CTGGGTTAAT")
            self.assertEqual(record.table_rows[5][1], "83")
            self.assertEqual(record.table_rows[5][2], "11584.1")
            self.assertEqual(len(record.table_rows[6]), 3)
            self.assertEqual(record.table_rows[6][0], "GTTGTGGTTA")
            self.assertEqual(record.table_rows[6][1], "82")
            self.assertEqual(record.table_rows[6][2], "11444.5")
            self.assertEqual(len(record.table_rows[7]), 3)
            self.assertEqual(record.table_rows[7][0], "GATCCCAACT")
            self.assertEqual(record.table_rows[7][1], "63")
            self.assertEqual(record.table_rows[7][2], "8792.74")
            self.assertEqual(len(record.table_rows[8]), 3)
            self.assertEqual(record.table_rows[8][0], "TGCAGTCACT")
            self.assertEqual(record.table_rows[8][1], "59")
            self.assertEqual(record.table_rows[8][2], "8234.47")
            self.assertEqual(len(record.table_rows[9]), 3)
            self.assertEqual(record.table_rows[9][0], "GGATTTGGCC")
            self.assertEqual(record.table_rows[9][1], "58")
            self.assertEqual(record.table_rows[9][2], "8094.91")
            self.assertEqual(len(record.table_rows[10]), 3)
            self.assertEqual(record.table_rows[10][0], "GGGCTGGGGT")
            self.assertEqual(record.table_rows[10][1], "56")
            self.assertEqual(record.table_rows[10][2], "7815.77")
            self.assertEqual(len(record.table_rows[11]), 3)
            self.assertEqual(record.table_rows[11][0], "ATAATTCTTT")
            self.assertEqual(record.table_rows[11][1], "44")
            self.assertEqual(record.table_rows[11][2], "6140.96")
            self.assertEqual(len(record.table_rows[12]), 3)
            self.assertEqual(record.table_rows[12][0], "CTTCCTTGCC")
            self.assertEqual(record.table_rows[12][1], "42")
            self.assertEqual(record.table_rows[12][2], "5861.83")
            self.assertEqual(len(record.table_rows[13]), 3)
            self.assertEqual(record.table_rows[13][0], "TTGGTCCTCT")
            self.assertEqual(record.table_rows[13][1], "40")
            self.assertEqual(record.table_rows[13][2], "5582.69")
            self.assertEqual(len(record.table_rows[14]), 3)
            self.assertEqual(record.table_rows[14][0], "GGCAAGCCCC")
            self.assertEqual(record.table_rows[14][1], "36")
            self.assertEqual(record.table_rows[14][2], "5024.42")
            self.assertEqual(len(record.table_rows[15]), 3)
            self.assertEqual(record.table_rows[15][0], "AACTAAAAAA")
            self.assertEqual(record.table_rows[15][1], "34")
            self.assertEqual(record.table_rows[15][2], "4745.29")
            self.assertEqual(len(record.table_rows[16]), 3)
            self.assertEqual(record.table_rows[16][0], "AGGGCTTCCA")
            self.assertEqual(record.table_rows[16][1], "34")
            self.assertEqual(record.table_rows[16][2], "4745.29")
            self.assertEqual(len(record.table_rows[17]), 3)
            self.assertEqual(record.table_rows[17][0], "AGGCTACGGA")
            self.assertEqual(record.table_rows[17][1], "33")
            self.assertEqual(record.table_rows[17][2], "4605.72")
            self.assertEqual(len(record.table_rows[18]), 3)
            self.assertEqual(record.table_rows[18][0], "GTGAAACCCC")
            self.assertEqual(record.table_rows[18][1], "32")
            self.assertEqual(record.table_rows[18][2], "4466.15")
            self.assertEqual(len(record.table_rows[19]), 3)
            self.assertEqual(record.table_rows[19][0], "AACTAACAAA")
            self.assertEqual(record.table_rows[19][1], "31")
            self.assertEqual(record.table_rows[19][2], "4326.59")
            self.assertEqual(len(record.table_rows[20]), 3)
            self.assertEqual(record.table_rows[20][0], "GAAAAATGGT")
            self.assertEqual(record.table_rows[20][1], "30")
            self.assertEqual(record.table_rows[20][2], "4187.02")

    def test_soft_ex_family(self):
        path = "Geo/soft_ex_family.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "PLATFORM")
            self.assertEqual(
                record.entity_id, "Murine 15K long oligo array version 2.0"
            )
            self.assertEqual(len(record.entity_attributes), 12)
            self.assertEqual(
                record.entity_attributes["Platform_title"],
                "Murine 15K long oligo array version 2.0",
            )
            self.assertEqual(
                record.entity_attributes["Platform_web_link"],
                "http://www.microarray.protocols.html",
            )
            self.assertEqual(record.entity_attributes["platform_table_end"], "")
            self.assertEqual(record.entity_attributes["Platform_support"], "glass")
            self.assertEqual(
                record.entity_attributes["Platform_manufacturer"],
                "Un. London microarray facility",
            )
            self.assertEqual(record.entity_attributes["Platform_coating"], "polysine")
            self.assertEqual(
                record.entity_attributes["Platform_technology"],
                "spotted oligonucleotide",
            )
            self.assertEqual(record.entity_attributes["platform_table_begin"], "")
            self.assertEqual(
                len(record.entity_attributes["Platform_manufacture_protocol"]), 12
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][0],
                "1.  Oligos are arrayed in Greiner 384-well flat-bottom plates. Each well contains 600 pmol of 70-mer oligo.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][1],
                "2. Resuspend oligos in water to 20 uM and rearray 5 \xb5L into 384-well, Genetix polystyrene V-bottom plates (cat# X6004).",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][2],
                "3. Allow Genetix plates to dry through passive water evaporation in a protected environment (e.g., chemical hood).",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][3],
                "4. Before printing, add 5 \xb5L of 1X Printing Buffer to each well. This can be done the night before a print run is started.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][4],
                "5. Seal plates with Corning seals.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][5],
                "6. Incubate at 37\xb0C for 30 minutes to aid resuspension of DNA.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][6],
                "7. Shake plates near maximum rotational speed on flat-bed shaker for 1 minute.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][7],
                "8. Centrifuge plates at 2000 rpm for 3 minutes.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][8],
                "9. Remove seals and cover with plate lids. Place in appropriate location of plate cassette. This should be done with first plates just before print run is started to minimize evaporation time before printing. For second and third cassettes, wait until 30 minutes before next cassette is needed to begin centrifugation.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][9],
                "10. Make sure plates rest behind both holding clips in the cassettes. Push plates back into the cassettes as far as they will go, putting them in the proper position for the server arm.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][10],
                "11. After the print run is completed, allow plates to dry through passive evaporation in a protected environment.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_manufacture_protocol"][11],
                "12. For each subsequent preparation of these plates for a print run, add water to the wells instead of sodium phosphate buffer. The amount of water should be decreased by 0.25 \xb5L per print run, as this is the amount drawn up by the pin capillary during each dip.",
            )
            self.assertEqual(
                record.entity_attributes["Platform_organism"], "Mus musculus"
            )
            self.assertEqual(len(record.entity_attributes["Platform_contributor"]), 5)
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][0], "Jane,Doe"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][1], "John,A,Smith"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][2], "Hans,van Elton"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][3], "John,Smithers Jr"
            )
            self.assertEqual(
                record.entity_attributes["Platform_contributor"][4], "Jie,D,Chen"
            )
            self.assertEqual(
                record.entity_attributes["Platform_distribution"], "non-commercial"
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["Gene_Desc"], "Gene description")
            self.assertEqual(record.col_defs["SEQUENCE"], "Probe sequence information")
            self.assertEqual(record.col_defs["Gene_Sym"], "Gene symbols")
            self.assertEqual(
                record.col_defs["GB_ACC"],
                'GenBank accession number of sequence used to design oligonucleotide probe   LINK_PRE:"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&term="',
            )
            self.assertEqual(record.col_defs["SPOT_ID"], "alternative identifier")
            self.assertEqual(record.col_defs["ID"], "")
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID")
            self.assertEqual(record.table_rows[0][1], "GB_ACC")
            self.assertEqual(record.table_rows[0][2], "Gene_Desc")
            self.assertEqual(record.table_rows[0][3], "Gene_Sym")
            self.assertEqual(record.table_rows[0][4], "SPOT_ID")
            self.assertEqual(record.table_rows[0][5], "SEQUENCE")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "U02079")
            self.assertEqual(
                record.table_rows[1][2],
                "nuclear factor of activated T-cells, cytoplasmic 2",
            )
            self.assertEqual(record.table_rows[1][3], "Nfatc2")
            self.assertEqual(record.table_rows[1][4], "")
            self.assertEqual(
                record.table_rows[1][5],
                "ACCTGGATGACGCAGCCACTTCAGAAAGCTGGGTTGGGACAGAAAGGTATATAGAGAGAAAATTTTGGAA",
            )
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "NM_008154")
            self.assertEqual(record.table_rows[2][2], "G-protein coupled receptor 3")
            self.assertEqual(record.table_rows[2][3], "Gpr3")
            self.assertEqual(record.table_rows[2][4], "")
            self.assertEqual(
                record.table_rows[2][5],
                "CTGTACAATGCTCTCACTTACTACTCAGAGACAACGGTAACTCGGACTTATGTGATGCTGGCCTTGGTGT",
            )
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "AK015719")
            self.assertEqual(record.table_rows[3][2], "tropomodulin 2")
            self.assertEqual(record.table_rows[3][3], "Tmod2")
            self.assertEqual(record.table_rows[3][4], "")
            self.assertEqual(
                record.table_rows[3][5],
                "CACCAGGCTCAGTGCCTAGTATCGGCTTCACCTAGTGTGGTTACTCAGGGCACGCAGAGCTACAGAACAC",
            )
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "AK003367")
            self.assertEqual(
                record.table_rows[4][2], "mitochondrial ribosomal protein L15"
            )
            self.assertEqual(record.table_rows[4][3], "Mrpl15")
            self.assertEqual(record.table_rows[4][4], "")
            self.assertEqual(
                record.table_rows[4][5],
                "CAAGAAGTCTAGAAATTCTGTGCAAGCCTATTCCATTCTTTCTGCGGGGACAACCAATTCCGAAAAGAAT",
            )
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "BC003333")
            self.assertEqual(record.table_rows[5][2], "RIKEN cDNA 0610033I05 gene")
            self.assertEqual(record.table_rows[5][3], "0610033I05Rik")
            self.assertEqual(record.table_rows[5][4], "")
            self.assertEqual(
                record.table_rows[5][5],
                "AGAACTGGGTGGCAGATATCCTAGAGTTTTGACCAACGTTCACAGCACACATATTGATCTTATAGGACCT",
            )
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "NM_008462")
            self.assertEqual(
                record.table_rows[6][2],
                "killer cell lectin-like receptor, subfamily A, member 2",
            )
            self.assertEqual(record.table_rows[6][3], "Klra2")
            self.assertEqual(record.table_rows[6][4], "")
            self.assertEqual(
                record.table_rows[6][5],
                "TGAATTGAAGTTCCTTAAATCCCAACTTCAAAGAAACACATACTGGATTTCACTGACACATCATAAAAGC",
            )
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "NM_008029")
            self.assertEqual(record.table_rows[7][2], "FMS-like tyrosine kinase 4")
            self.assertEqual(record.table_rows[7][3], "Flt4")
            self.assertEqual(record.table_rows[7][4], "")
            self.assertEqual(
                record.table_rows[7][5],
                "GAGGTGCTGTGGGATGACCGCCGGGGCATGCGGGTGCCCACTCAACTGTTGCGCGATGCCCTGTACCTGC",
            )
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "NM_054088")
            self.assertEqual(record.table_rows[8][2], "adiponutrin")
            self.assertEqual(record.table_rows[8][3], "Adpn")
            self.assertEqual(record.table_rows[8][4], "")
            self.assertEqual(
                record.table_rows[8][5],
                "GTCTGAGTTCCATTCCAAAGACGAAGTCGTGGATGCCCTGGTGTGTTCCTGCTTCATTCCCCTCTTCTCT",
            )
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "NM_009750")
            self.assertEqual(
                record.table_rows[9][2],
                "nerve growth factor receptor (TNFRSF16) associated protein 1",
            )
            self.assertEqual(record.table_rows[9][3], "Ngfrap1")
            self.assertEqual(record.table_rows[9][4], "")
            self.assertEqual(
                record.table_rows[9][5],
                "TACAGCTGAGAAATTGTCTACGCATCCTTATGGGGGAGCTGTCTAACCACCACGATCACCATGATGAATT",
            )
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "AB045323")
            self.assertEqual(
                record.table_rows[10][2], "DNA segment, Chr 8, ERATO Doi 594, expressed"
            )
            self.assertEqual(record.table_rows[10][3], "D8Ertd594e")
            self.assertEqual(record.table_rows[10][4], "")
            self.assertEqual(
                record.table_rows[10][5],
                "GATTCAGACTCGGGAGGAGCATCCCAACCTCTCCTTGAGGATAAAGGCCTGAGCGATTGCCCTGGGGAGC",
            )
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "AK005789")
            self.assertEqual(
                record.table_rows[11][2], "dynein, cytoplasmic, light chain 2B"
            )
            self.assertEqual(record.table_rows[11][3], "Dncl2b")
            self.assertEqual(record.table_rows[11][4], "")
            self.assertEqual(
                record.table_rows[11][5],
                "TGCAGAAGGCATTCCAATCCGAACAACCCTGGACAACTCCACAACGGTTCAGTATGCGGGTCTTCTCCAC",
            )
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "NM_010517")
            self.assertEqual(
                record.table_rows[12][2], "insulin-like growth factor binding protein 4"
            )
            self.assertEqual(record.table_rows[12][3], "Igfbp4")
            self.assertEqual(record.table_rows[12][4], "")
            self.assertEqual(
                record.table_rows[12][5],
                "GGAGAAGCTGGCGCGCTGCCGCCCCCCCGTGGGTTGCGAGGAGTTGGTGCGGGAGCCAGGCTGCGGTTGT",
            )
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "AK010722")
            self.assertEqual(record.table_rows[13][2], "RIKEN cDNA 2410075D05 gene")
            self.assertEqual(record.table_rows[13][3], "2410075D05Rik")
            self.assertEqual(record.table_rows[13][4], "")
            self.assertEqual(
                record.table_rows[13][5],
                "GGAGCATCTGGAGTTCCGCTTACCGGAAATAAAGTCTTTACTATCGGTGATTGGAGGGCAGTTCACTAAC",
            )
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "AK003755")
            self.assertEqual(
                record.table_rows[14][2], "DNA segment, Chr 4, ERATO Doi 421, expressed"
            )
            self.assertEqual(record.table_rows[14][3], "D4Ertd421e")
            self.assertEqual(record.table_rows[14][4], "")
            self.assertEqual(
                record.table_rows[14][5],
                "AGCAAAGAGATCTCCCTCAGTGTGCCCATAGGTGGCGGTGCGAGCTTGCGGTTATTGGCCAGTGACTTGC",
            )
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "BC003241")
            self.assertEqual(
                record.table_rows[15][2],
                "cleavage stimulation factor, 3' pre-RNA, subunit 3",
            )
            self.assertEqual(record.table_rows[15][3], "Cstf3")
            self.assertEqual(record.table_rows[15][4], "")
            self.assertEqual(
                record.table_rows[15][5],
                "AAATTAGAAGAAAATCCATATGACCTTGATGCTTGGAGCATTCTCATTCGAGAGGCACAGAATCAACCTA",
            )
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "AK004937")
            self.assertEqual(record.table_rows[16][2], "RIKEN cDNA 1300007O09 gene")
            self.assertEqual(record.table_rows[16][3], "1300007O09Rik")
            self.assertEqual(record.table_rows[16][4], "")
            self.assertEqual(
                record.table_rows[16][5],
                "CAGACACAAACCCTAGGTTGTATTGTAGACCGGAGTTTAAGCAGGCACTACCTGTCTGTCTTTTCTTCAT",
            )
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "AK004524")
            self.assertEqual(
                record.table_rows[17][2],
                "unnamed protein product; hypothetical SOCS domain",
            )
            self.assertEqual(record.table_rows[17][3], "")
            self.assertEqual(record.table_rows[17][4], "")
            self.assertEqual(
                record.table_rows[17][5],
                "CGGAGCCCTGCGCGCCCAGAGCCCCCTCCCACCCGCTTCCACCAAGTGCATGGAGCCAACATCCGCATGG",
            )
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "NM_025999")
            self.assertEqual(record.table_rows[18][2], "RIKEN cDNA 2610110L04 gene")
            self.assertEqual(record.table_rows[18][3], "2610110L04Rik")
            self.assertEqual(record.table_rows[18][4], "")
            self.assertEqual(
                record.table_rows[18][5],
                "TGCATTGATAAATGGAGTGATCGACACAGGAACTGCCCCATTTGTCGCCTACAGATGACTGGAGCAAATG",
            )
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "")
            self.assertEqual(record.table_rows[19][2], "")
            self.assertEqual(record.table_rows[19][3], "")
            self.assertEqual(record.table_rows[19][4], "-- CONTROL")
            self.assertEqual(record.table_rows[19][5], "")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "NM_023120")
            self.assertEqual(
                record.table_rows[20][2],
                "guanine nucleotide binding protein (G protein), beta polypeptide 1-like",
            )
            self.assertEqual(record.table_rows[20][3], "Gnb1l")
            self.assertEqual(record.table_rows[20][4], "")
            self.assertEqual(
                record.table_rows[20][5],
                "ACCGCCTGGTCCCAGATTTGTCCTCCGAGGCACACAGTCGGCTGTGAACACGCTCCATTTCTGCCCACCA",
            )
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(
                record.entity_id, "Control Embyronic Stem Cell Replicate 1"
            )
            self.assertEqual(len(record.entity_attributes), 24)
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch2"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_platform_id"],
                "Murine 15K long oligo array version 2.0",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "Control Embyronic Stem Cell Replicate 1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"], "file1.gpr"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch2"], "Mus musculus"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"], "Mus musculus"
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "Cy5")
            self.assertEqual(len(record.entity_attributes["Sample_scan_protocol"]), 2)
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][0],
                "Scanned on an Agilent G2565AA scanner.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][1],
                "Images were quantified using Agilent Feature Extraction Software (version A.7.5).",
            )
            self.assertEqual(record.entity_attributes["sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch2"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal.",
            )
            self.assertEqual(record.entity_attributes["sample_table_end"], "")
            self.assertEqual(record.entity_attributes["Sample_label_ch2"], "Cy3")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Biological replicate 1 of 4. Control embryonic stem cells, untreated, harvested after several passages.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch2"],
                "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch2"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 4
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "ES-D3 cell line (CRL-1934)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1], "Age: day 4"
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][2],
                "Tissue: blastocytes",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][3],
                "Strain: 129/Sv mice",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch2"]), 3
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][0],
                "Strain: C57BL/6",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][1],
                "Age: e17.5 d",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][2],
                "Tissue: whole embryo",
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"],
                "log(REDsignal/GREENsignal) per feature (processed signals used).",
            )
            self.assertEqual(
                record.col_defs["gProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.',
            )
            self.assertEqual(
                record.col_defs["LogRatioError"],
                "error of the log ratio calculated according to the error model chosen.",
            )
            self.assertEqual(
                record.col_defs["PValueLogRatio"],
                "Significance level of the Log Ratio computed for a feature.",
            )
            self.assertEqual(
                record.col_defs["rProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LogRatioError")
            self.assertEqual(record.table_rows[0][3], "PValueLogRatio")
            self.assertEqual(record.table_rows[0][4], "gProcessedSignal")
            self.assertEqual(record.table_rows[0][5], "rProcessedSignal")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "-1.6274758")
            self.assertEqual(record.table_rows[1][2], "1.36E-01")
            self.assertEqual(record.table_rows[1][3], "6.41E-33")
            self.assertEqual(record.table_rows[1][4], "9.13E+03")
            self.assertEqual(record.table_rows[1][5], "2.15E+02")
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "0.1412248")
            self.assertEqual(record.table_rows[2][2], "1.34E+00")
            self.assertEqual(record.table_rows[2][3], "1.00E+00")
            self.assertEqual(record.table_rows[2][4], "4.14E+01")
            self.assertEqual(record.table_rows[2][5], "5.72E+01")
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.1827684")
            self.assertEqual(record.table_rows[3][2], "5.19E-02")
            self.assertEqual(record.table_rows[3][3], "4.33E-04")
            self.assertEqual(record.table_rows[3][4], "5.13E+03")
            self.assertEqual(record.table_rows[3][5], "7.81E+03")
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.3932267")
            self.assertEqual(record.table_rows[4][2], "6.08E-02")
            self.assertEqual(record.table_rows[4][3], "1.02E-10")
            self.assertEqual(record.table_rows[4][4], "4.65E+03")
            self.assertEqual(record.table_rows[4][5], "1.88E+03")
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "-0.9865994")
            self.assertEqual(record.table_rows[5][2], "1.05E-01")
            self.assertEqual(record.table_rows[5][3], "6.32E-21")
            self.assertEqual(record.table_rows[5][4], "2.91E+03")
            self.assertEqual(record.table_rows[5][5], "3.01E+02")
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "0.0238812")
            self.assertEqual(record.table_rows[6][2], "1.02E-01")
            self.assertEqual(record.table_rows[6][3], "8.15E-01")
            self.assertEqual(record.table_rows[6][4], "7.08E+02")
            self.assertEqual(record.table_rows[6][5], "7.48E+02")
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-1.4841822")
            self.assertEqual(record.table_rows[7][2], "1.25E-01")
            self.assertEqual(record.table_rows[7][3], "1.42E-32")
            self.assertEqual(record.table_rows[7][4], "1.02E+04")
            self.assertEqual(record.table_rows[7][5], "3.36E+02")
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "-1.8261356")
            self.assertEqual(record.table_rows[8][2], "4.15E-01")
            self.assertEqual(record.table_rows[8][3], "1.10E-05")
            self.assertEqual(record.table_rows[8][4], "7.19E+02")
            self.assertEqual(record.table_rows[8][5], "1.07E+01")
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "-1.0344779")
            self.assertEqual(record.table_rows[9][2], "1.78E+00")
            self.assertEqual(record.table_rows[9][3], "1.00E+00")
            self.assertEqual(record.table_rows[9][4], "9.62E+01")
            self.assertEqual(record.table_rows[9][5], "8.89E+00")
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "0.2405891")
            self.assertEqual(record.table_rows[10][2], "3.09E-01")
            self.assertEqual(record.table_rows[10][3], "4.36E-01")
            self.assertEqual(record.table_rows[10][4], "1.61E+02")
            self.assertEqual(record.table_rows[10][5], "2.80E+02")
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "0.3209366")
            self.assertEqual(record.table_rows[11][2], "3.59E-01")
            self.assertEqual(record.table_rows[11][3], "3.71E-01")
            self.assertEqual(record.table_rows[11][4], "1.25E+02")
            self.assertEqual(record.table_rows[11][5], "2.61E+02")
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "0.358304")
            self.assertEqual(record.table_rows[12][2], "2.06E+00")
            self.assertEqual(record.table_rows[12][3], "1.00E+00")
            self.assertEqual(record.table_rows[12][4], "2.04E+01")
            self.assertEqual(record.table_rows[12][5], "4.66E+01")
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "-0.0122072")
            self.assertEqual(record.table_rows[13][2], "3.64E-01")
            self.assertEqual(record.table_rows[13][3], "9.73E-01")
            self.assertEqual(record.table_rows[13][4], "1.84E+02")
            self.assertEqual(record.table_rows[13][5], "1.79E+02")
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "-1.5480396")
            self.assertEqual(record.table_rows[14][2], "1.30E-01")
            self.assertEqual(record.table_rows[14][3], "7.21E-33")
            self.assertEqual(record.table_rows[14][4], "1.02E+04")
            self.assertEqual(record.table_rows[14][5], "2.90E+02")
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "0.0073419")
            self.assertEqual(record.table_rows[15][2], "2.98E-01")
            self.assertEqual(record.table_rows[15][3], "9.80E-01")
            self.assertEqual(record.table_rows[15][4], "2.21E+02")
            self.assertEqual(record.table_rows[15][5], "2.25E+02")
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "-0.2267015")
            self.assertEqual(record.table_rows[16][2], "9.44E-01")
            self.assertEqual(record.table_rows[16][3], "8.10E-01")
            self.assertEqual(record.table_rows[16][4], "8.90E+01")
            self.assertEqual(record.table_rows[16][5], "5.28E+01")
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "-0.1484023")
            self.assertEqual(record.table_rows[17][2], "8.01E-01")
            self.assertEqual(record.table_rows[17][3], "8.53E-01")
            self.assertEqual(record.table_rows[17][4], "9.65E+01")
            self.assertEqual(record.table_rows[17][5], "6.86E+01")
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.6122195")
            self.assertEqual(record.table_rows[18][2], "1.28E-01")
            self.assertEqual(record.table_rows[18][3], "1.69E-06")
            self.assertEqual(record.table_rows[18][4], "1.12E+03")
            self.assertEqual(record.table_rows[18][5], "2.73E+02")
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.0796905")
            self.assertEqual(record.table_rows[19][2], "8.78E-02")
            self.assertEqual(record.table_rows[19][3], "3.64E-01")
            self.assertEqual(record.table_rows[19][4], "8.21E+02")
            self.assertEqual(record.table_rows[19][5], "9.87E+02")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "-0.084895")
            self.assertEqual(record.table_rows[20][2], "9.38E-01")
            self.assertEqual(record.table_rows[20][3], "9.28E-01")
            self.assertEqual(record.table_rows[20][4], "7.68E+01")
            self.assertEqual(record.table_rows[20][5], "6.32E+01")
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(
                record.entity_id, "Control Embyronic Stem Cell Replicate 2"
            )
            self.assertEqual(len(record.entity_attributes), 24)
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch2"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_platform_id"],
                "Murine 15K long oligo array version 2.0",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "Control Embyronic Stem Cell Replicate 2",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"], "file2.gpr"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch2"], "Mus musculus"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"], "Mus musculus"
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "Cy5")
            self.assertEqual(len(record.entity_attributes["Sample_scan_protocol"]), 2)
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][0],
                "Scanned on an Agilent G2565AA scanner.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"][1],
                "Images were quantified using Agilent Feature Extraction Software (version A.7.5).",
            )
            self.assertEqual(record.entity_attributes["sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch2"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal.",
            )
            self.assertEqual(record.entity_attributes["sample_table_end"], "")
            self.assertEqual(record.entity_attributes["Sample_label_ch2"], "Cy3")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Biological replicate 2 of 4. Control embryonic stem cells, untreated, harvested after several passages.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Total RNA from murine ES-D3 embryonic stem cells labeled with Cyanine-5 (red).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch2"],
                "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch2"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 4
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "ES-D3 cell line (CRL-1934)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1], "Age: day 4"
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][2],
                "Tissue: blastocytes",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][3],
                "Strain: 129/Sv mice",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch2"]), 3
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][0],
                "Strain: C57BL/6",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][1],
                "Age: e17.5 d",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][2],
                "Tissue: whole embryo",
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"],
                "log(REDsignal/GREENsignal) per feature (processed signals used).",
            )
            self.assertEqual(
                record.col_defs["gProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.',
            )
            self.assertEqual(
                record.col_defs["LogRatioError"],
                "error of the log ratio calculated according to the error model chosen.",
            )
            self.assertEqual(
                record.col_defs["PValueLogRatio"],
                "Significance level of the Log Ratio computed for a feature.",
            )
            self.assertEqual(
                record.col_defs["rProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LogRatioError")
            self.assertEqual(record.table_rows[0][3], "PValueLogRatio")
            self.assertEqual(record.table_rows[0][4], "gProcessedSignal")
            self.assertEqual(record.table_rows[0][5], "rProcessedSignal")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "-1.1697263")
            self.assertEqual(record.table_rows[1][2], "1.23E-01")
            self.assertEqual(record.table_rows[1][3], "2.14E-21")
            self.assertEqual(record.table_rows[1][4], "3.17E+03")
            self.assertEqual(record.table_rows[1][5], "2.14E+02")
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "-0.1111353")
            self.assertEqual(record.table_rows[2][2], "1.63E+00")
            self.assertEqual(record.table_rows[2][3], "9.46E-01")
            self.assertEqual(record.table_rows[2][4], "5.43E+01")
            self.assertEqual(record.table_rows[2][5], "4.20E+01")
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.1400597")
            self.assertEqual(record.table_rows[3][2], "5.11E-02")
            self.assertEqual(record.table_rows[3][3], "6.17E-03")
            self.assertEqual(record.table_rows[3][4], "6.72E+03")
            self.assertEqual(record.table_rows[3][5], "9.28E+03")
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.4820633")
            self.assertEqual(record.table_rows[4][2], "6.38E-02")
            self.assertEqual(record.table_rows[4][3], "4.06E-14")
            self.assertEqual(record.table_rows[4][4], "6.46E+03")
            self.assertEqual(record.table_rows[4][5], "2.13E+03")
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "-1.2116196")
            self.assertEqual(record.table_rows[5][2], "1.22E-01")
            self.assertEqual(record.table_rows[5][3], "2.31E-23")
            self.assertEqual(record.table_rows[5][4], "3.62E+03")
            self.assertEqual(record.table_rows[5][5], "2.22E+02")
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "-0.0230528")
            self.assertEqual(record.table_rows[6][2], "1.04E-01")
            self.assertEqual(record.table_rows[6][3], "8.24E-01")
            self.assertEqual(record.table_rows[6][4], "8.76E+02")
            self.assertEqual(record.table_rows[6][5], "8.31E+02")
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-1.1380152")
            self.assertEqual(record.table_rows[7][2], "1.13E-01")
            self.assertEqual(record.table_rows[7][3], "9.23E-24")
            self.assertEqual(record.table_rows[7][4], "3.94E+03")
            self.assertEqual(record.table_rows[7][5], "2.86E+02")
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "-1.834596")
            self.assertEqual(record.table_rows[8][2], "5.40E-01")
            self.assertEqual(record.table_rows[8][3], "6.74E-04")
            self.assertEqual(record.table_rows[8][4], "6.44E+02")
            self.assertEqual(record.table_rows[8][5], "9.43E+00")
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "-0.9747637")
            self.assertEqual(record.table_rows[9][2], "2.14E+00")
            self.assertEqual(record.table_rows[9][3], "1.00E+00")
            self.assertEqual(record.table_rows[9][4], "9.17E+01")
            self.assertEqual(record.table_rows[9][5], "9.72E+00")
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "0.3874005")
            self.assertEqual(record.table_rows[10][2], "2.92E-01")
            self.assertEqual(record.table_rows[10][3], "1.85E-01")
            self.assertEqual(record.table_rows[10][4], "1.69E+02")
            self.assertEqual(record.table_rows[10][5], "4.11E+02")
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "0.5340442")
            self.assertEqual(record.table_rows[11][2], "3.29E-01")
            self.assertEqual(record.table_rows[11][3], "1.04E-01")
            self.assertEqual(record.table_rows[11][4], "1.23E+02")
            self.assertEqual(record.table_rows[11][5], "4.20E+02")
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "0.3260696")
            self.assertEqual(record.table_rows[12][2], "1.92E+00")
            self.assertEqual(record.table_rows[12][3], "8.65E-01")
            self.assertEqual(record.table_rows[12][4], "2.73E+01")
            self.assertEqual(record.table_rows[12][5], "5.77E+01")
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "0.3010618")
            self.assertEqual(record.table_rows[13][2], "2.84E-01")
            self.assertEqual(record.table_rows[13][3], "2.90E-01")
            self.assertEqual(record.table_rows[13][4], "1.93E+02")
            self.assertEqual(record.table_rows[13][5], "3.87E+02")
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "-1.0760413")
            self.assertEqual(record.table_rows[14][2], "1.08E-01")
            self.assertEqual(record.table_rows[14][3], "1.63E-23")
            self.assertEqual(record.table_rows[14][4], "4.06E+03")
            self.assertEqual(record.table_rows[14][5], "3.41E+02")
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "-0.1167371")
            self.assertEqual(record.table_rows[15][2], "3.87E-01")
            self.assertEqual(record.table_rows[15][3], "7.63E-01")
            self.assertEqual(record.table_rows[15][4], "2.32E+02")
            self.assertEqual(record.table_rows[15][5], "1.77E+02")
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "-0.1936322")
            self.assertEqual(record.table_rows[16][2], "9.44E-01")
            self.assertEqual(record.table_rows[16][3], "8.38E-01")
            self.assertEqual(record.table_rows[16][4], "1.02E+02")
            self.assertEqual(record.table_rows[16][5], "6.56E+01")
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "-0.3275898")
            self.assertEqual(record.table_rows[17][2], "7.87E-01")
            self.assertEqual(record.table_rows[17][3], "6.77E-01")
            self.assertEqual(record.table_rows[17][4], "1.41E+02")
            self.assertEqual(record.table_rows[17][5], "6.65E+01")
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.4805853")
            self.assertEqual(record.table_rows[18][2], "1.14E-01")
            self.assertEqual(record.table_rows[18][3], "2.41E-05")
            self.assertEqual(record.table_rows[18][4], "1.34E+03")
            self.assertEqual(record.table_rows[18][5], "4.42E+02")
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.1109524")
            self.assertEqual(record.table_rows[19][2], "9.56E-02")
            self.assertEqual(record.table_rows[19][3], "2.46E-01")
            self.assertEqual(record.table_rows[19][4], "8.38E+02")
            self.assertEqual(record.table_rows[19][5], "1.08E+03")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "0.1677912")
            self.assertEqual(record.table_rows[20][2], "6.51E-01")
            self.assertEqual(record.table_rows[20][3], "7.97E-01")
            self.assertEqual(record.table_rows[20][4], "9.84E+01")
            self.assertEqual(record.table_rows[20][5], "1.45E+02")
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(
                record.entity_id,
                "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1",
            )
            self.assertEqual(len(record.entity_attributes), 25)
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch2"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Oligoarray control targets and hybridization buffer (Agilent In Situ Hybridization Kit Plus) were added, and samples were applied to microarrays enclosed in Agilent SureHyb-enabled hybridization chambers. After hybridization, slides were washed sequentially with 6x SSC/0.005% Triton X-102 and 0.1x SSC/0.005% Triton X-102 before scanning. Slides were hybridized for 17 h at 60\xb0C in a rotating oven, and washed.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "TriZol procedure",
            )
            self.assertEqual(
                record.entity_attributes["Sample_platform_id"],
                "Murine 15K long oligo array version 2.0",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"], "file3.gpr"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch2"], "Mus musculus"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"], "Mus musculus"
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "Cy5")
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "Scanned on an Agilent G2565AA scanner.",
            )
            self.assertEqual(record.entity_attributes["sample_table_begin"], "")
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch2"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP, with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "10 \xb5g of total RNA were primed with 2 \xb5l of 100 \xb5M T16N2 DNA primer at 70\xb0C for 10 min, then reversed transcribed at 42\xb0C for 1 h in the presence of 400 U SuperScript II RTase (Invitrogen), and 100 \xb5M each dATP, dTTP, dGTP with 25 \xb5M dCTP, 25 \xb5M Cy5-labeled dCTP (NEN Life Science, Boston, MA), and RNase inhibitor (Invitrogen). RNA was then degraded with RNase A, and labeled cDNAs were purified using QIAquick PCR columns (Qiagen).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "LOWESS normalized, background subtracted VALUE data obtained from log of processed Red signal/processed Green signal.",
            )
            self.assertEqual(record.entity_attributes["sample_table_end"], "")
            self.assertEqual(record.entity_attributes["Sample_label_ch2"], "Cy3")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Biological replicate 1 of 3. Stable triple-fusion-reporter-gene transfected embryonic stem cells, harvested after several passages.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Total RNA from murine ES-D3 triple-transfected embryonic stem cells labeled with Cyanine-5 (red).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch2"],
                "Total RNA from pooled whole mouse embryos e17.5, labeled with Cyanine-3 (green).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch2"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "ES cells were kept in an undifferentiated, pluripotent state by using 1000 IU/ml leukemia inhibitory factor (LIF; Chemicon, ESGRO, ESG1107), and grown on top of murine embryonic fibroblasts feeder layer inactivated by 10 ug/ml of mitomycin C (Sigma, St. Louis). ES cells were cultured on 0.1% gelatin-coated plastic dishes in ES medium containing Dulbecco modified Eagle medium supplemented with 15% fetal calf serum, 0.1 mM beta-mercaptoethanol, 2 mM glutamine, and 0.1 mN non-essential amino acids.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 5
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "ES-D3 cell line (CRL-1934)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Transfected with pUb-fluc-mrfp-ttk triple fusion reporter gene.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][2], "Age: day 4"
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][3],
                "Tissue: blastocytes",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][4],
                "Strain: 129/Sv mice",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch2"]), 3
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][0],
                "Strain: C57BL/6",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][1],
                "Age: e17.5 d",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch2"][2],
                "Tissue: whole embryo",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "PCR amplification and standard cloning techniques were used to insert fluc and mrfp genes from plasmids pCDNA 3.1-CMV-fluc (Promega, Madison, WI) and pCDNA3.1-CMV-mrfp in frame with the ttk gene into the pCDNA3.1-truncated	sr39tk. This triple fusion (TF) reporter gene fragment (3.3 kbp) was released from the plasmid with Not1 and BamH1 restriction enzymes before blunt-end ligation into the multiple cloning site of lentiviral transfer vector, FUG, driven by the human ubiquitin-C promoter. Self-inactivating (SIN) lentivirus was prepared by transient transfection of 293T cells. Briefly, pFUG-TF containing the triple fusion reporter gene was co-transfected into 293T cells with HIV-1 packaging vector (?8.9) and vesicular stomatitis virus G glycoprotein-pseudotyped envelop vector (pVSVG). Lentivirus supernatant was concentrated by sediment centrifugation using a SW29 rotor at 50,000 x g for two hours. Concentrated virus was titered on 293T cells. Murine ES cells were transfected with LV-pUb-fluc-mrfp-ttk at a multiplicity of infection (MOI) of 10.",
            )
            self.assertEqual(len(record.col_defs), 6)
            self.assertEqual(record.col_defs["ID_REF"], "")
            self.assertEqual(
                record.col_defs["VALUE"],
                "log(REDsignal/GREENsignal) per feature (processed signals used).",
            )
            self.assertEqual(
                record.col_defs["gProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," green "channel," used for computation of log ratio.',
            )
            self.assertEqual(
                record.col_defs["LogRatioError"],
                "error of the log ratio calculated according to the error model chosen.",
            )
            self.assertEqual(
                record.col_defs["PValueLogRatio"],
                "Significance level of the Log Ratio computed for a feature.",
            )
            self.assertEqual(
                record.col_defs["rProcessedSignal"],
                'Dye-normalized signal after surrogate "algorithm," red "channel," used for computation of log ratio.',
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 6)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LogRatioError")
            self.assertEqual(record.table_rows[0][3], "PValueLogRatio")
            self.assertEqual(record.table_rows[0][4], "gProcessedSignal")
            self.assertEqual(record.table_rows[0][5], "rProcessedSignal")
            self.assertEqual(len(record.table_rows[1]), 6)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "-0.7837546")
            self.assertEqual(record.table_rows[1][2], "1.30E-01")
            self.assertEqual(record.table_rows[1][3], "1.70E-09")
            self.assertEqual(record.table_rows[1][4], "2.10E+03")
            self.assertEqual(record.table_rows[1][5], "3.46E+02")
            self.assertEqual(len(record.table_rows[2]), 6)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "0.3797837")
            self.assertEqual(record.table_rows[2][2], "1.15E+00")
            self.assertEqual(record.table_rows[2][3], "7.41E-01")
            self.assertEqual(record.table_rows[2][4], "5.59E+01")
            self.assertEqual(record.table_rows[2][5], "1.34E+02")
            self.assertEqual(len(record.table_rows[3]), 6)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.2079269")
            self.assertEqual(record.table_rows[3][2], "5.38E-02")
            self.assertEqual(record.table_rows[3][3], "1.12E-04")
            self.assertEqual(record.table_rows[3][4], "5.04E+03")
            self.assertEqual(record.table_rows[3][5], "8.14E+03")
            self.assertEqual(len(record.table_rows[4]), 6)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.4730291")
            self.assertEqual(record.table_rows[4][2], "6.71E-02")
            self.assertEqual(record.table_rows[4][3], "1.86E-12")
            self.assertEqual(record.table_rows[4][4], "5.66E+03")
            self.assertEqual(record.table_rows[4][5], "1.91E+03")
            self.assertEqual(len(record.table_rows[5]), 6)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "-0.9481128")
            self.assertEqual(record.table_rows[5][2], "1.19E-01")
            self.assertEqual(record.table_rows[5][3], "1.30E-15")
            self.assertEqual(record.table_rows[5][4], "3.10E+03")
            self.assertEqual(record.table_rows[5][5], "3.49E+02")
            self.assertEqual(len(record.table_rows[6]), 6)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "-0.0159867")
            self.assertEqual(record.table_rows[6][2], "1.33E-01")
            self.assertEqual(record.table_rows[6][3], "9.05E-01")
            self.assertEqual(record.table_rows[6][4], "8.45E+02")
            self.assertEqual(record.table_rows[6][5], "8.14E+02")
            self.assertEqual(len(record.table_rows[7]), 6)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-0.819922")
            self.assertEqual(record.table_rows[7][2], "1.14E-01")
            self.assertEqual(record.table_rows[7][3], "7.01E-13")
            self.assertEqual(record.table_rows[7][4], "2.75E+03")
            self.assertEqual(record.table_rows[7][5], "4.16E+02")
            self.assertEqual(len(record.table_rows[8]), 6)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "-0.1559774")
            self.assertEqual(record.table_rows[8][2], "9.16E-01")
            self.assertEqual(record.table_rows[8][3], "8.65E-01")
            self.assertEqual(record.table_rows[8][4], "1.34E+02")
            self.assertEqual(record.table_rows[8][5], "9.34E+01")
            self.assertEqual(len(record.table_rows[9]), 6)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "0.145267")
            self.assertEqual(record.table_rows[9][2], "3.90E+00")
            self.assertEqual(record.table_rows[9][3], "1.00E+00")
            self.assertEqual(record.table_rows[9][4], "2.22E+01")
            self.assertEqual(record.table_rows[9][5], "3.10E+01")
            self.assertEqual(len(record.table_rows[10]), 6)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "0.3611211")
            self.assertEqual(record.table_rows[10][2], "3.40E-01")
            self.assertEqual(record.table_rows[10][3], "2.88E-01")
            self.assertEqual(record.table_rows[10][4], "1.97E+02")
            self.assertEqual(record.table_rows[10][5], "4.52E+02")
            self.assertEqual(len(record.table_rows[11]), 6)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "0.5092089")
            self.assertEqual(record.table_rows[11][2], "4.39E-01")
            self.assertEqual(record.table_rows[11][3], "2.46E-01")
            self.assertEqual(record.table_rows[11][4], "1.24E+02")
            self.assertEqual(record.table_rows[11][5], "4.01E+02")
            self.assertEqual(len(record.table_rows[12]), 6)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "0.3715387")
            self.assertEqual(record.table_rows[12][2], "1.69E+00")
            self.assertEqual(record.table_rows[12][3], "8.26E-01")
            self.assertEqual(record.table_rows[12][4], "3.84E+01")
            self.assertEqual(record.table_rows[12][5], "9.04E+01")
            self.assertEqual(len(record.table_rows[13]), 6)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "0.1734934")
            self.assertEqual(record.table_rows[13][2], "3.57E-01")
            self.assertEqual(record.table_rows[13][3], "6.27E-01")
            self.assertEqual(record.table_rows[13][4], "2.37E+02")
            self.assertEqual(record.table_rows[13][5], "3.53E+02")
            self.assertEqual(len(record.table_rows[14]), 6)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "-0.9340707")
            self.assertEqual(record.table_rows[14][2], "1.20E-01")
            self.assertEqual(record.table_rows[14][3], "6.90E-15")
            self.assertEqual(record.table_rows[14][4], "2.96E+03")
            self.assertEqual(record.table_rows[14][5], "3.45E+02")
            self.assertEqual(len(record.table_rows[15]), 6)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "-0.2956317")
            self.assertEqual(record.table_rows[15][2], "5.78E-01")
            self.assertEqual(record.table_rows[15][3], "6.09E-01")
            self.assertEqual(record.table_rows[15][4], "2.46E+02")
            self.assertEqual(record.table_rows[15][5], "1.25E+02")
            self.assertEqual(len(record.table_rows[16]), 6)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "-0.2321102")
            self.assertEqual(record.table_rows[16][2], "1.22E+00")
            self.assertEqual(record.table_rows[16][3], "8.49E-01")
            self.assertEqual(record.table_rows[16][4], "1.09E+02")
            self.assertEqual(record.table_rows[16][5], "6.37E+01")
            self.assertEqual(len(record.table_rows[17]), 6)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "-0.1603561")
            self.assertEqual(record.table_rows[17][2], "1.16E+00")
            self.assertEqual(record.table_rows[17][3], "8.90E-01")
            self.assertEqual(record.table_rows[17][4], "1.06E+02")
            self.assertEqual(record.table_rows[17][5], "7.34E+01")
            self.assertEqual(len(record.table_rows[18]), 6)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.5063897")
            self.assertEqual(record.table_rows[18][2], "1.63E-01")
            self.assertEqual(record.table_rows[18][3], "1.95E-03")
            self.assertEqual(record.table_rows[18][4], "1.15E+03")
            self.assertEqual(record.table_rows[18][5], "3.58E+02")
            self.assertEqual(len(record.table_rows[19]), 6)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.1990761")
            self.assertEqual(record.table_rows[19][2], "1.32E-01")
            self.assertEqual(record.table_rows[19][3], "1.32E-01")
            self.assertEqual(record.table_rows[19][4], "6.65E+02")
            self.assertEqual(record.table_rows[19][5], "1.05E+03")
            self.assertEqual(len(record.table_rows[20]), 6)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "0.2985912")
            self.assertEqual(record.table_rows[20][2], "8.89E-01")
            self.assertEqual(record.table_rows[20][3], "7.37E-01")
            self.assertEqual(record.table_rows[20][4], "8.06E+01")
            self.assertEqual(record.table_rows[20][5], "1.60E+02")
            record = next(records)
            self.assertEqual(record.entity_type, "SERIES")
            self.assertEqual(record.entity_id, "Murine ES Cells")
            self.assertEqual(len(record.entity_attributes), 7)
            self.assertEqual(len(record.entity_attributes["Series_sample_id"]), 3)
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][0],
                "Control Embyronic Stem Cell Replicate 1",
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][1],
                "Control Embyronic Stem Cell Replicate 2",
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][2],
                "Triple-Fusion Transfected Embryonic Stem Cells Replicate 1",
            )
            self.assertEqual(record.entity_attributes["Series_pubmed_id"], "16390873")
            self.assertEqual(len(record.entity_attributes["Series_contributor"]), 9)
            self.assertEqual(
                record.entity_attributes["Series_contributor"][0], "Joseph,C,Wu"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][1], "Joshua,M,Spin"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][2], "Feng,,Cao"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][3], "Shaun,,Lin"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][4], "Olivier,,Gheysens"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][5], "Ian,Y,Chen"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][6], "Anya,,Tsalenko"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][7], "Sanjiv,S,Ghambhir"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][8], "Thomas,,Quertermous"
            )
            self.assertEqual(
                record.entity_attributes["Series_summary"],
                "Transcriptional profiling of mouse embryonic stem cells comparing control untreated ES cells with ES cells transfected with a pUb-fluc-mrfp-ttk triple fusion reporter gene. The latter makes ES visualization possible by FACS and single ce",
            )
            self.assertEqual(
                record.entity_attributes["Series_type"], "Genetic modification"
            )
            self.assertEqual(
                record.entity_attributes["Series_title"],
                "Murine ES Cells: Control vs. Triple-Fusion Transfected",
            )
            self.assertEqual(
                record.entity_attributes["Series_overall_design"],
                "Two-condition experiment, ES vs. TF-ES cells. Biological replicates: 4 control, 3 transfected, independently grown and harvested. One replicate per array.",
            )
            self.assertEqual(len(record.col_defs), 0)
            self.assertEqual(len(record.table_rows), 0)

    def test_GSM804(self):
        path = "Geo/GSM804.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "GSM804")
            self.assertEqual(len(record.entity_attributes), 18)
            self.assertEqual(record.entity_attributes["Sample_pubmed_id"], "11687795")
            self.assertEqual(
                record.entity_attributes["Sample_submitter_institute"],
                "University of California San Francisco",
            )
            self.assertEqual(len(record.entity_attributes["Sample_author"]), 19)
            self.assertEqual(
                record.entity_attributes["Sample_author"][0], "Antoine,M,Snijders"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][1], "Norma,,Nowak"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][2], "Richard,,Segraves"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][3], "Stephanie,,Blackwood"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][4], "Nils,,Brown"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][5], "Jeffery,,Conroy"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][6], "Greg,,Hamilton"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][7], "Anna,K,Hindle"
            )
            self.assertEqual(record.entity_attributes["Sample_author"][8], "Bing,,Huey")
            self.assertEqual(
                record.entity_attributes["Sample_author"][9], "Karen,,Kimura"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][10], "Sindy,,Law"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][11], "Ken,,Myambo"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][12], "Joel,,Palmer"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][13], "Bauke,,Ylstra"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][14], "Jingzhu,P,Yue"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][15], "Joe,W,Gray"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][16], "Ajay,N,Jain"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][17], "Daniel,,Pinkel"
            )
            self.assertEqual(
                record.entity_attributes["Sample_author"][18], "Donna,G,Albertson"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_phone"], "415 502-8463"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_department"],
                "Comprehensive Cancer Center",
            )
            self.assertEqual(len(record.entity_attributes["Sample_description"]), 4)
            self.assertEqual(
                record.entity_attributes["Sample_description"][0],
                'Coriell Cell Repositories cell line <a href="http://locus.umdnj.edu/nigms/nigms_cgi/display.cgi?GM05296">GM05296</a>.',
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][1],
                "Fibroblast cell line derived from a 1 month old female with multiple congenital malformations, dysmorphic features, intrauterine growth retardation, heart murmur, cleft palate, equinovarus deformity, microcephaly, coloboma of right iris, clinodactyly, reduced RBC catalase activity, and 1 copy of catalase gene.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][2],
                "Chromosome abnormalities are present.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_description"][3],
                "Karyotype is 46,XX,-11,+der(11)inv ins(11;10)(11pter> 11p13::10q21>10q24::11p13>11qter)mat",
            )
            self.assertEqual(
                record.entity_attributes["Sample_target_source2"],
                "normal male reference genomic DNA",
            )
            self.assertEqual(
                record.entity_attributes["Sample_target_source1"], "Cell line GM05296"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_name"], "Donna,G,Albertson"
            )
            self.assertEqual(record.entity_attributes["Sample_platform_id"], "GPL28")
            self.assertEqual(
                record.entity_attributes["Sample_type"], "dual channel genomic"
            )
            self.assertEqual(
                record.entity_attributes["Sample_status"], "Public on Feb 12 2002"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_email"],
                "albertson@cc.ucsf.edu",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"], "CGH_Albertson_GM05296-001218"
            )
            self.assertEqual(
                record.entity_attributes["Sample_organism"], "Homo sapiens"
            )
            self.assertEqual(record.entity_attributes["Sample_series_id"], "GSE16")
            self.assertEqual(
                record.entity_attributes["Sample_submission_date"], "Jan 17 2002"
            )
            self.assertEqual(
                record.entity_attributes["Sample_submitter_city"],
                "San Francisco,CA,94143,USA",
            )
            self.assertEqual(len(record.col_defs), 5)
            self.assertEqual(
                record.col_defs["NO_REPLICATES"],
                "Number of replicate spot measurements",
            )
            self.assertEqual(
                record.col_defs["LOG2STDDEV"], "Standard deviation of VALUE"
            )
            self.assertEqual(
                record.col_defs["ID_REF"],
                "Unique row identifier, genome position order",
            )
            self.assertEqual(
                record.col_defs["VALUE"],
                "aka LOG2RATIO, mean of log base 2 of LINEAR_RATIO",
            )
            self.assertEqual(
                record.col_defs["LINEAR_RATIO"], "Mean of replicate Cy3/Cy5 ratios"
            )
            self.assertEqual(len(record.table_rows), 21)
            self.assertEqual(len(record.table_rows[0]), 5)
            self.assertEqual(record.table_rows[0][0], "ID_REF")
            self.assertEqual(record.table_rows[0][1], "VALUE")
            self.assertEqual(record.table_rows[0][2], "LINEAR_RATIO")
            self.assertEqual(record.table_rows[0][3], "LOG2STDDEV")
            self.assertEqual(record.table_rows[0][4], "NO_REPLICATES")
            self.assertEqual(len(record.table_rows[1]), 5)
            self.assertEqual(record.table_rows[1][0], "1")
            self.assertEqual(record.table_rows[1][1], "")
            self.assertEqual(record.table_rows[1][2], "1.047765")
            self.assertEqual(record.table_rows[1][3], "0.011853")
            self.assertEqual(record.table_rows[1][4], "3")
            self.assertEqual(len(record.table_rows[2]), 5)
            self.assertEqual(record.table_rows[2][0], "2")
            self.assertEqual(record.table_rows[2][1], "")
            self.assertEqual(record.table_rows[2][2], "")
            self.assertEqual(record.table_rows[2][3], "")
            self.assertEqual(record.table_rows[2][4], "0")
            self.assertEqual(len(record.table_rows[3]), 5)
            self.assertEqual(record.table_rows[3][0], "3")
            self.assertEqual(record.table_rows[3][1], "0.008824")
            self.assertEqual(record.table_rows[3][2], "1.006135")
            self.assertEqual(record.table_rows[3][3], "0.00143")
            self.assertEqual(record.table_rows[3][4], "3")
            self.assertEqual(len(record.table_rows[4]), 5)
            self.assertEqual(record.table_rows[4][0], "4")
            self.assertEqual(record.table_rows[4][1], "-0.000894")
            self.assertEqual(record.table_rows[4][2], "0.99938")
            self.assertEqual(record.table_rows[4][3], "0.001454")
            self.assertEqual(record.table_rows[4][4], "3")
            self.assertEqual(len(record.table_rows[5]), 5)
            self.assertEqual(record.table_rows[5][0], "5")
            self.assertEqual(record.table_rows[5][1], "0.075875")
            self.assertEqual(record.table_rows[5][2], "1.054")
            self.assertEqual(record.table_rows[5][3], "0.003077")
            self.assertEqual(record.table_rows[5][4], "3")
            self.assertEqual(len(record.table_rows[6]), 5)
            self.assertEqual(record.table_rows[6][0], "6")
            self.assertEqual(record.table_rows[6][1], "0.017303")
            self.assertEqual(record.table_rows[6][2], "1.012066")
            self.assertEqual(record.table_rows[6][3], "0.005876")
            self.assertEqual(record.table_rows[6][4], "2")
            self.assertEqual(len(record.table_rows[7]), 5)
            self.assertEqual(record.table_rows[7][0], "7")
            self.assertEqual(record.table_rows[7][1], "-0.006766")
            self.assertEqual(record.table_rows[7][2], "0.995321")
            self.assertEqual(record.table_rows[7][3], "0.013881")
            self.assertEqual(record.table_rows[7][4], "3")
            self.assertEqual(len(record.table_rows[8]), 5)
            self.assertEqual(record.table_rows[8][0], "8")
            self.assertEqual(record.table_rows[8][1], "0.020755")
            self.assertEqual(record.table_rows[8][2], "1.014491")
            self.assertEqual(record.table_rows[8][3], "0.005506")
            self.assertEqual(record.table_rows[8][4], "3")
            self.assertEqual(len(record.table_rows[9]), 5)
            self.assertEqual(record.table_rows[9][0], "9")
            self.assertEqual(record.table_rows[9][1], "-0.094938")
            self.assertEqual(record.table_rows[9][2], "0.936313")
            self.assertEqual(record.table_rows[9][3], "0.012662")
            self.assertEqual(record.table_rows[9][4], "3")
            self.assertEqual(len(record.table_rows[10]), 5)
            self.assertEqual(record.table_rows[10][0], "10")
            self.assertEqual(record.table_rows[10][1], "-0.054527")
            self.assertEqual(record.table_rows[10][2], "0.96291")
            self.assertEqual(record.table_rows[10][3], "0.01073")
            self.assertEqual(record.table_rows[10][4], "3")
            self.assertEqual(len(record.table_rows[11]), 5)
            self.assertEqual(record.table_rows[11][0], "11")
            self.assertEqual(record.table_rows[11][1], "-0.025057")
            self.assertEqual(record.table_rows[11][2], "0.982782")
            self.assertEqual(record.table_rows[11][3], "0.003855")
            self.assertEqual(record.table_rows[11][4], "3")
            self.assertEqual(len(record.table_rows[12]), 5)
            self.assertEqual(record.table_rows[12][0], "12")
            self.assertEqual(record.table_rows[12][1], "")
            self.assertEqual(record.table_rows[12][2], "")
            self.assertEqual(record.table_rows[12][3], "")
            self.assertEqual(record.table_rows[12][4], "0")
            self.assertEqual(len(record.table_rows[13]), 5)
            self.assertEqual(record.table_rows[13][0], "13")
            self.assertEqual(record.table_rows[13][1], "0.108454")
            self.assertEqual(record.table_rows[13][2], "1.078072")
            self.assertEqual(record.table_rows[13][3], "0.005196")
            self.assertEqual(record.table_rows[13][4], "3")
            self.assertEqual(len(record.table_rows[14]), 5)
            self.assertEqual(record.table_rows[14][0], "14")
            self.assertEqual(record.table_rows[14][1], "0.078633")
            self.assertEqual(record.table_rows[14][2], "1.056017")
            self.assertEqual(record.table_rows[14][3], "0.009165")
            self.assertEqual(record.table_rows[14][4], "3")
            self.assertEqual(len(record.table_rows[15]), 5)
            self.assertEqual(record.table_rows[15][0], "15")
            self.assertEqual(record.table_rows[15][1], "0.098571")
            self.assertEqual(record.table_rows[15][2], "1.070712")
            self.assertEqual(record.table_rows[15][3], "0.007834")
            self.assertEqual(record.table_rows[15][4], "3")
            self.assertEqual(len(record.table_rows[16]), 5)
            self.assertEqual(record.table_rows[16][0], "16")
            self.assertEqual(record.table_rows[16][1], "0.044048")
            self.assertEqual(record.table_rows[16][2], "1.031003")
            self.assertEqual(record.table_rows[16][3], "0.013651")
            self.assertEqual(record.table_rows[16][4], "3")
            self.assertEqual(len(record.table_rows[17]), 5)
            self.assertEqual(record.table_rows[17][0], "17")
            self.assertEqual(record.table_rows[17][1], "0.018039")
            self.assertEqual(record.table_rows[17][2], "1.012582")
            self.assertEqual(record.table_rows[17][3], "0.005471")
            self.assertEqual(record.table_rows[17][4], "3")
            self.assertEqual(len(record.table_rows[18]), 5)
            self.assertEqual(record.table_rows[18][0], "18")
            self.assertEqual(record.table_rows[18][1], "-0.088807")
            self.assertEqual(record.table_rows[18][2], "0.9403")
            self.assertEqual(record.table_rows[18][3], "0.010571")
            self.assertEqual(record.table_rows[18][4], "3")
            self.assertEqual(len(record.table_rows[19]), 5)
            self.assertEqual(record.table_rows[19][0], "19")
            self.assertEqual(record.table_rows[19][1], "0.016349")
            self.assertEqual(record.table_rows[19][2], "1.011397")
            self.assertEqual(record.table_rows[19][3], "0.007113")
            self.assertEqual(record.table_rows[19][4], "3")
            self.assertEqual(len(record.table_rows[20]), 5)
            self.assertEqual(record.table_rows[20][0], "20")
            self.assertEqual(record.table_rows[20][1], "0.030977")
            self.assertEqual(record.table_rows[20][2], "1.021704")
            self.assertEqual(record.table_rows[20][3], "0.016798")
            self.assertEqual(record.table_rows[20][4], "3")

    def test_soft_ex_affy_chp(self):
        path = "Geo/soft_ex_affy_chp.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "Drosophila_T0-1")
            self.assertEqual(len(record.entity_attributes), 16)
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"],
                "Drosophila melanogaster",
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "biotin")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 2
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "Genotype: yellow white and Oregon R parents",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "Trizol extraction of total RNA was performed according to the manufacturer's instructions.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Drosophila embryos before nuclear cycle 9 (maternal transcripts)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_table"], "Drosophila_T0-1.CHP"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "embryo at T0, biological rep1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"],
                "Drosophila_T0-1.CEL",
            )
            self.assertEqual(len(record.col_defs), 0)
            self.assertEqual(len(record.table_rows), 0)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "Drosophila_T0-2")
            self.assertEqual(len(record.entity_attributes), 16)
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"],
                "Drosophila melanogaster",
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "biotin")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Gene expression data from embryos younger than nuclear cycle 9, i.e. before zygotic genome activation.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 2
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "Genotype: yellow white and Oregon R parents",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Age: embryos younger than nuclear cycle 9, i.e. before pole cells budding",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "Trizol extraction of total RNA was performed according to the manufacturer's instructions.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Drosophila embryos before nuclear cycle 9 (maternal transcripts)",
            )
            self.assertEqual(
                record.entity_attributes["Sample_table"], "Drosophila_T0-2.CHP"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "embryo at T0, biological rep2",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"],
                "Drosophila_T0-2.CEL",
            )
            self.assertEqual(len(record.col_defs), 0)
            self.assertEqual(len(record.table_rows), 0)
            record = next(records)
            self.assertEqual(record.entity_type, "SAMPLE")
            self.assertEqual(record.entity_id, "Drosophila_T1-1")
            self.assertEqual(len(record.entity_attributes), 16)
            self.assertEqual(
                record.entity_attributes["Sample_organism_ch1"],
                "Drosophila melanogaster",
            )
            self.assertEqual(record.entity_attributes["Sample_label_ch1"], "biotin")
            self.assertEqual(
                record.entity_attributes["Sample_description"],
                "Gene expression data from embryos in slow phase of cellularisation.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_growth_protocol_ch1"],
                "30 min egg collections of OreR and yw flies at 25C were aged at room temperature (RT) according to the different temporal classes T0-T4.",
            )
            self.assertEqual(
                len(record.entity_attributes["Sample_characteristics_ch1"]), 2
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][0],
                "Genotype: yellow white and Oregon R parents",
            )
            self.assertEqual(
                record.entity_attributes["Sample_characteristics_ch1"][1],
                "Age: embryos in slow phase of cellularisation",
            )
            self.assertEqual(
                record.entity_attributes["Sample_scan_protocol"],
                "GeneChips were scanned using the Hewlett-Packard GeneArray Scanner G2500A.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_hyb_protocol"],
                "Following fragmentation, 10 microg of cRNA were hybridized for 16 hr at 45C on GeneChip Drosophila Genome Array. GeneChips were washed and stained in the Affymetrix Fluidics Station 400.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_extract_protocol_ch1"],
                "Trizol extraction of total RNA was performed according to the manufacturer's instructions.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_source_name_ch1"],
                "Drosophila embryos in slow phase of cellularisation",
            )
            self.assertEqual(
                record.entity_attributes["Sample_table"], "Drosophila_T1-1.CHP"
            )
            self.assertEqual(
                record.entity_attributes["Sample_molecule_ch1"], "total RNA"
            )
            self.assertEqual(
                record.entity_attributes["Sample_label_protocol_ch1"],
                "Biotinylated cRNA were prepared according to the standard Affymetrix protocol from 6 microg total RNA (Expression Analysis Technical Manual, 2001, Affymetrix).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_data_processing"],
                "The data were analyzed with Microarray Suite version 5.0 (MAS 5.0) using Affymetrix default analysis settings and global scaling as normalization method. The trimmed mean target intensity of each array was arbitrarily set to 100.",
            )
            self.assertEqual(
                record.entity_attributes["Sample_treatment_protocol_ch1"],
                "Embryos were dechorionated with 50% bleach, put on a cover slip and covered with Halocarbon oil 27 (Sigma). Embryos of the appropriate stage were manually selected under the dissecting scope. Selected embryos were transferred to a basket, rinsed with PBS with 0,7% NaCl, 0,04% triton-X100 and placed on ice in the Trizol solution (GibcoBRL).",
            )
            self.assertEqual(
                record.entity_attributes["Sample_title"],
                "embryo at T1, biological rep1",
            )
            self.assertEqual(
                record.entity_attributes["Sample_supplementary_file"],
                "Drosophila_T1-1.CEL",
            )
            self.assertEqual(len(record.col_defs), 0)
            self.assertEqual(len(record.table_rows), 0)
            record = next(records)
            self.assertEqual(record.entity_type, "SERIES")
            self.assertEqual(record.entity_id, "Dros_embryo_timecourse")
            self.assertEqual(len(record.entity_attributes), 6)
            self.assertEqual(len(record.entity_attributes["Series_sample_id"]), 3)
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][0], "Drosophila_T0-1"
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][1], "Drosophila_T0-2"
            )
            self.assertEqual(
                record.entity_attributes["Series_sample_id"][2], "Drosophila_T1-1"
            )
            self.assertEqual(len(record.entity_attributes["Series_contributor"]), 5)
            self.assertEqual(
                record.entity_attributes["Series_contributor"][0], "Jane,Doe"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][1], "John,A,Smith"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][2], "Hans,van Elton"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][3], "John,Smithers Jr"
            )
            self.assertEqual(
                record.entity_attributes["Series_contributor"][4], "Jie,D,Chen"
            )
            self.assertEqual(len(record.entity_attributes["Series_summary"]), 2)
            self.assertEqual(
                record.entity_attributes["Series_summary"][0],
                "Morphogenesis of epithelial tissues relies on the precise developmental control of cell polarity and architecture. In the early Drosophila embryo, the primary epithelium forms during cellularisation, following a tightly controlled genetic programme where specific sets of genes are up-regulated. Some of them, for instance, control membrane invagination between the nuclei anchored at the apical surface of the syncytium.",
            )
            self.assertEqual(
                record.entity_attributes["Series_summary"][1],
                "We used microarrays to detail the global programme of gene expression underlying cellularisation and identified distinct classes of up-regulated genes during this process.",
            )
            self.assertEqual(record.entity_attributes["Series_type"], "time course")
            self.assertEqual(
                record.entity_attributes["Series_title"],
                "Expression data from early Drosophila embryo",
            )
            self.assertEqual(
                record.entity_attributes["Series_overall_design"],
                "Drosophila embryos were selected at successive stages of early development for RNA extraction and hybridization on Affymetrix microarrays. We sought to obtain homogeneous populations of embryos at each developmental stage in order to increase the temporal resolution of expression profiles. To that end, we hand-selected embryos according to morphological criteria at five time-points: before pole cell formation, i.e. before zygotic transcription (T0), during the slow phase (T1) and the fast phase (T2) of cellularisation and at the beginning (T3) and the end (T4) of gastrulation.",
            )
            self.assertEqual(len(record.col_defs), 0)
            self.assertEqual(len(record.table_rows), 0)

    def test_record_str(self):
        path = "Geo/GSM804.txt"
        with open(path, encoding="latin") as handle:
            records = Geo.parse(handle)
            record = next(records)
            self.assertEqual(
                str(record),
                """\
GEO Type: SAMPLE
GEO Id: GSM804
Sample_author: Antoine,M,Snijders

Sample_author: Norma,,Nowak

Sample_author: Richard,,Segraves

Sample_author: Stephanie,,Blackwood

Sample_author: Nils,,Brown

Sample_author: Jeffery,,Conroy

Sample_author: Greg,,Hamilton

Sample_author: Anna,K,Hindle

Sample_author: Bing,,Huey

Sample_author: Karen,,Kimura

Sample_author: Sindy,,Law

Sample_author: Ken,,Myambo

Sample_author: Joel,,Palmer

Sample_author: Bauke,,Ylstra

Sample_author: Jingzhu,P,Yue

Sample_author: Joe,W,Gray

Sample_author: Ajay,N,Jain

Sample_author: Daniel,,Pinkel

Sample_author: Donna,G,Albertson

Sample_description: Coriell Cell Repositories cell line <a h
ref="http://locus.umdnj.edu/nigms/nigms_cgi/display.cgi?GM05296">GM05296</a>.

Sample_description: Fibroblast cell line derived from a 1 mo
nth old female with multiple congenital malformations, dysmorphic features, intr
auterine growth retardation, heart murmur, cleft palate, equinovarus deformity, \
\nmicrocephaly, coloboma of right iris, clinodactyly, reduced RBC catalase activit
y, and 1 copy of catalase gene.

Sample_description: Chromosome abnormalities are present.

Sample_description: Karyotype is 46,XX,-11,+der(11)inv ins(1
1;10)(11pter> 11p13::10q21>10q24::11p13>11qter)mat

Sample_organism: Homo sapiens

Sample_platform_id: GPL28

Sample_pubmed_id: 11687795

Sample_series_id: GSE16

Sample_status: Public on Feb 12 2002

Sample_submission_date: Jan 17 2002

Sample_submitter_city: San Francisco,CA,94143,USA

Sample_submitter_department: Comprehensive Cancer Center

Sample_submitter_email: albertson@cc.ucsf.edu

Sample_submitter_institute: University of California San Francisco

Sample_submitter_name: Donna,G,Albertson

Sample_submitter_phone: 415 502-8463

Sample_target_source1: Cell line GM05296

Sample_target_source2: normal male reference genomic DNA

Sample_title: CGH_Albertson_GM05296-001218

Sample_type: dual channel genomic

Column Header Definitions
    ID_REF: Unique row identifier, genome position o
    rder

    LINEAR_RATIO: Mean of replicate Cy3/Cy5 ratios

    LOG2STDDEV: Standard deviation of VALUE

    NO_REPLICATES: Number of replicate spot measurements

    VALUE: aka LOG2RATIO, mean of log base 2 of LIN
    EAR_RATIO

0: ID_REF	VALUE	LINEAR_RATIO	LOG2STDDEV	NO_REPLICATES\t
1: 1		1.047765	0.011853	3\t
2: 2				0\t
3: 3	0.008824	1.006135	0.00143	3\t
4: 4	-0.000894	0.99938	0.001454	3\t
5: 5	0.075875	1.054	0.003077	3\t
6: 6	0.017303	1.012066	0.005876	2\t
7: 7	-0.006766	0.995321	0.013881	3\t
8: 8	0.020755	1.014491	0.005506	3\t
9: 9	-0.094938	0.936313	0.012662	3\t
10: 10	-0.054527	0.96291	0.01073	3\t
11: 11	-0.025057	0.982782	0.003855	3\t
12: 12				0\t
13: 13	0.108454	1.078072	0.005196	3\t
14: 14	0.078633	1.056017	0.009165	3\t
15: 15	0.098571	1.070712	0.007834	3\t
16: 16	0.044048	1.031003	0.013651	3\t
17: 17	0.018039	1.012582	0.005471	3\t
18: 18	-0.088807	0.9403	0.010571	3\t
19: 19	0.016349	1.011397	0.007113	3\t
20: 20	0.030977	1.021704	0.016798	3\t
""",
            )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
