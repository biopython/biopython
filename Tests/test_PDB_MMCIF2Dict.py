# Copyright 2017 by Francesco Gastaldello. All rights reserved.
# Revisions copyright 2017 by Peter Cock.  All rights reserved.
#
# Converted by Francesco Gastaldello from an older unit test copyright 2002
# by Thomas Hamelryck.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Unit tests for the Bio.PDB.MMCIF2Dict module."""

import unittest

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.PDB.")

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

import io
import textwrap


class MMCIF2dictTests(unittest.TestCase):

    def test_MMCIF2dict(self):
        filename = "PDB/1A8O.cif"
        mmcif = MMCIF2Dict(filename)
        self.assertEqual(len(mmcif.keys()), 575)
        self.assertEqual(mmcif['_entity_poly_seq.mon_id'], ['MSE', 'ASP', 'ILE', 'ARG', 'GLN', 'GLY', 'PRO', 'LYS', 'GLU', 'PRO', 'PHE', 'ARG', 'ASP', 'TYR', 'VAL', 'ASP', 'ARG', 'PHE', 'TYR', 'LYS', 'THR', 'LEU', 'ARG', 'ALA', 'GLU', 'GLN', 'ALA', 'SER', 'GLN', 'GLU', 'VAL', 'LYS', 'ASN', 'TRP', 'MSE', 'THR', 'GLU', 'THR', 'LEU', 'LEU', 'VAL', 'GLN', 'ASN', 'ALA', 'ASN', 'PRO', 'ASP', 'CYS', 'LYS', 'THR', 'ILE', 'LEU', 'LYS', 'ALA', 'LEU', 'GLY', 'PRO', 'GLY', 'ALA', 'THR', 'LEU', 'GLU', 'GLU', 'MSE', 'MSE', 'THR', 'ALA', 'CYS', 'GLN', 'GLY'])
        self.assertEqual(mmcif['_atom_site.Cartn_x'], ['19.594', '20.255', '20.351', '19.362', '19.457', '20.022', '21.718', '21.424', '21.554', '21.835', '21.947', '21.678', '23.126', '23.098', '23.433', '22.749', '22.322', '22.498', '21.220', '20.214', '23.062', '24.282', '23.423', '25.429', '21.280', '20.173', '20.766', '21.804', '19.444', '18.724', '18.011', '17.416', '16.221', '15.459', '15.824', '20.116', '20.613', '20.546', '19.488', '19.837', '20.385', '19.526', '18.365', '20.090', '21.675', '21.698', '20.859', '20.729', '20.260', '19.435', '20.158', '19.512', '18.993', '20.056', '20.300', '21.486', '22.285', '23.286', '24.155', '23.025', '22.117', '21.236', '20.159', '19.231', '23.152', '24.037', '23.563', '22.398', '24.086', '25.003', '24.858', '23.861', '25.748', '24.459', '24.089', '23.580', '24.111', '25.415', '26.116', '25.852', '22.544', '21.960', '22.965', '22.928', '20.793', '19.999', '19.234', '20.019', '18.495', '19.286', '18.523', '23.861', '24.870', '25.788', '26.158', '25.684', '26.777', '26.215', '27.235', '28.136', '28.155', '29.030', '26.137', '26.994', '26.279', '26.880', '27.408', '28.345', '28.814', '28.620', '24.992', '24.151', '24.025', '24.139', '22.787', '21.629', '21.657', '20.489', '20.571', '19.408', '19.450', '18.365', '23.839', '23.720', '24.962', '24.853', '23.502', '23.661', '22.120', '26.137', '27.387', '27.511', '27.925', '28.595', '28.723', '28.016', '29.545', '27.136', '27.202', '26.238', '26.585', '26.850', '27.835', '27.667', '26.352', '25.494', '25.797', '24.325', '25.037', '23.984', '24.456', '24.305', '22.761', '21.538', '21.301', '20.586', '20.130', '19.415', '19.186', '25.033', '25.526', '26.755', '27.015', '25.771', '24.608', '23.508', '24.583', '22.406', '23.490', '22.406', '21.326', '27.508', '28.691', '28.183', '28.705', '29.455', '30.787', '31.428', '32.618', '33.153', '27.116', '26.508', '25.826', '25.827', '25.475', '26.150', '24.741', '25.264', '24.587', '25.587', '25.302', '23.789', '22.707', '21.787', '21.910', '26.767', '27.806', '28.299', '28.656', '29.006', '28.944', '30.295', '30.744', '30.326', '29.441', '30.787', '28.332', '28.789', '27.943', '28.374', '28.803', '26.740', '25.833', '25.775', '24.998', '24.425', '24.354', '24.816', '24.535', '25.454', '26.601', '26.645', '25.240', '24.885', '27.391', '28.884', '29.200', '28.729', '29.998', '24.438', '23.066', '23.001', '23.824', '22.370', '22.035', '21.831', '21.174', '20.852', '20.917', '19.638', '20.949', '20.315', '18.908', '18.539', '20.262', '19.688', '20.414', '21.592', '19.714', '18.136', '16.775', '16.738', '15.875', '16.101', '15.478', '14.341', '13.247', '14.542', '17.668', '17.730', '18.064', '17.491', '18.754', '18.932', '18.279', '18.971', '19.343', '18.126', '17.905', '20.444', '21.777', '22.756', '24.069', '24.913', '17.344', '16.136', '15.146', '14.599', '15.468', '16.242', '17.164', '15.865', '14.932', '14.017', '14.495', '13.700', '13.904', '13.254', '12.332', '13.484', '11.975', '12.666', '14.303', '12.641', '14.280', '13.452', '15.793', '16.368', '16.285', '16.053', '17.815', '17.939', '17.221', '18.427', '16.438', '16.375', '14.950', '14.778', '16.869', '18.228', '16.791', '13.947', '12.529', '12.045', '11.151', '11.625', '11.950', '11.054', '11.086', '10.326', '12.589', '12.177', '13.076', '12.888', '11.978', '13.202', '10.883', '14.054', '14.963', '15.702', '15.846', '15.935', '15.286', '16.327', '14.580', '16.162', '16.876', '15.961', '16.391', '17.402', '18.238', '19.553', '18.506', '14.695', '13.703', '13.270', '13.262', '12.460', '11.372', '12.854', '12.954', '12.503', '13.541', '13.184', '12.008', '10.830', '10.505', '10.626', '10.093', '14.820', '15.887', '16.443', '17.416', '17.014', '16.627', '15.451', '17.619', '15.830', '16.248', '15.758', '14.809', '15.689', '16.404', '16.005', '14.639', '14.122', '17.109', '17.396', '16.559', '18.588', '14.018', '12.706', '12.516', '11.536', '12.617', '13.288', '14.522', '13.454', '13.383', '13.351', '12.406', '14.564', '14.482', '13.353', '15.552', '14.378', '14.488', '13.443', '12.968', '15.902', '16.144', '13.061', '12.087', '10.746', '10.157', '11.879', '11.014', '11.003', '10.171', '10.269', '10.273', '9.002', '9.101', '8.227', '8.612', '8.611', '7.224', '10.191', '10.458', '10.518', '9.916', '11.791', '11.677', '12.184', '12.967', '11.222', '11.377', '10.082', '9.885', '12.416', '13.824', '14.764', '14.287', '9.214', '7.937', '7.048', '6.294', '7.230', '7.828', '7.618', '8.090', '7.916', '7.189', '6.419', '6.871', '6.391', '6.449', '7.815', '8.305', '7.481', '7.371', '9.788', '10.832', '12.217', '10.789', '6.886', '6.080', '6.922', '8.149', '6.294', '7.024', '7.912', '7.680', '5.901', '4.734', '4.839', '8.952', '9.861', '10.886', '11.642', '10.910', '11.884', '13.285', '13.524', '11.599', '14.199', '15.563', '16.391', '16.022', '16.290', '16.498', '15.473', '17.509', '18.426', '18.875', '19.012', '19.645', '20.773', '20.264', '21.920', '19.082', '19.510', '18.471', '18.816', '19.784', '21.035', '20.954', '19.902', '21.955', '17.199', '16.109', '16.001', '15.690', '14.787', '14.776', '13.539', '13.220', '12.888', '16.301', '16.274', '17.413', '17.209', '16.429', '15.284', '15.332', '13.844', '18.606', '19.764', '19.548', '19.922', '21.047', '21.507', '23.105', '22.645', '18.915', '18.636', '17.640', '17.807', '18.050', '18.998', '17.730', '16.631', '15.593', '16.104', '15.685', '14.486', '17.033', '17.572', '18.985', '19.634', '17.525', '15.855', '19.451', '20.802', '21.001', '20.066', '21.152', '20.421', '20.725', '21.768', '19.817', '22.226', '22.536', '23.683', '24.328', '23.949', '15.165', '19.774', '22.152', '12.938', '23.499', '17.568', '13.544', '15.524', '31.249', '11.999', '14.511', '7.439', '19.303', '17.114', '21.867', '17.573', '26.151', '20.974', '20.796', '28.370', '29.565', '21.248', '25.744', '8.691', '30.789', '30.905', '28.623', '24.935', '23.462', '9.924', '28.729', '13.579', '23.652', '25.631', '17.799', '23.547', '16.363', '24.125', '33.063', '29.209', '10.391', '12.221', '18.997', '16.360', '27.915', '28.158', '21.975', '27.069', '30.148', '21.196', '8.864', '13.228', '18.577', '20.526', '25.758', '7.838', '20.569', '13.009', '19.229', '17.655', '30.445', '9.014', '3.398', '31.603', '16.543', '12.037', '7.261', '5.607', '23.532', '30.701', '32.300', '34.351', '9.450', '29.476', '13.681', '26.728', '10.004', '30.553', '23.569', '10.927', '17.983', '8.191', '32.095', '11.520', '13.249', '15.919', '11.187', '16.743'])
        self.assertEqual(
            mmcif['_struct_ref.pdbx_seq_one_letter_code'],
            textwrap.dedent('''\
            GARASVLSGGELDKWEKIRLRPGGKKQYKLKHIVWASRELERFAVNPGLLETSEGCRQILGQLQPSLQTGSEELRSLYNT
            IAVLYCVHQRIDVKDTKEALDKIEEEQNKSKKKAQQAAADTGNNSQVSQNYPIVQNLQGQMVHQAISPRTLNAWVKVVEE
            KAFSPEVIPMFSALSEGATPQDLNTMLNTVGGHQAAMQMLKETINEEAAEWDRLHPVHAGPIAPGQMREPRGSDIAGTTS
            TLQEQIGWMTHNPPIPVGEIYKRWIILGLNKIVRMYSPTSILDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETL
            LVQNANPDCKTILKALGPGATLEEMMTACQGVGGPGHKARVLAEAMSQVTNPATIMIQKGNFRNQRKTVKCFNCGKEGHI
            AKNCRAPRKKGCWKCGKEGHQMKDCTERQANFLGKIWPSHKGRPGNFLQSRPEPTAPPEESFRFGEETTTPSQKQEPIDK
            ELYPLASLRSLFGSDPSSQ'''))

    def test_underscores(self):
        # Test values starting with an underscore are not treated as keys
        filename = "PDB/4Q9R_min.cif"
        mmcif = MMCIF2Dict(filename)
        self.assertEqual(len(mmcif.keys()), 5)
        self.assertEqual(mmcif['_pdbx_audit_revision_item.item'], ['_atom_site.B_iso_or_equiv', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z'])

    def test_quotefix(self):
        # Test quote characters parse correctly
        filename = "PDB/1MOM_min.cif"
        mmcif = MMCIF2Dict(filename)
        self.assertEqual(len(mmcif.keys()), 21)
        self.assertEqual(mmcif['_struct_conf.pdbx_PDB_helix_id'], ['A', 'A\'', 'B', 'C', 'B\'', 'D', 'E', 'C\'', 'F', 'G', 'H', 'D\'', 'E\'', 'A\'"', 'BC', 'CD', 'DE'])

    def test_splitline(self):
        filename = "PDB/4Q9R_min.cif"
        mmcif = MMCIF2Dict(filename)
        self.assertEqual(list(mmcif._splitline("foo bar")), ["foo", "bar"])
        self.assertEqual(list(mmcif._splitline("  foo bar  ")), ["foo", "bar"])
        self.assertEqual(list(mmcif._splitline("'foo' bar")), ["foo", "bar"])
        self.assertEqual(list(mmcif._splitline("foo \"bar\"")), ["foo", "bar"])
        self.assertEqual(list(mmcif._splitline("foo 'bar a' b")), ["foo", "bar a", "b"])
        self.assertEqual(list(mmcif._splitline("foo 'bar'a' b")), ["foo", "bar'a", "b"])
        self.assertEqual(list(mmcif._splitline("foo \"bar' a\" b")), ["foo", "bar' a", "b"])
        self.assertEqual(list(mmcif._splitline("foo '' b")), ["foo", "", "b"])

        # A hash (#) starts a comment iff it is preceded by whitespace or is at
        # the beginning of a line:
        # https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax#lex
        self.assertEqual(list(mmcif._splitline("foo#bar")), ["foo#bar"])
        self.assertEqual(list(mmcif._splitline("foo #bar")), ["foo"])
        self.assertEqual(list(mmcif._splitline("foo# bar")), ["foo#", "bar"])
        self.assertEqual(list(mmcif._splitline("#foo bar")), [])

        self.assertRaises(ValueError, list, mmcif._splitline("foo 'bar"))
        self.assertRaises(ValueError, list, mmcif._splitline("foo 'ba'r  "))
        self.assertRaises(ValueError, list, mmcif._splitline("foo \"bar'"))
        self.assertRaises(ValueError, list, mmcif._splitline("foo b'ar'"))

    def test_verbatim_block(self):
        """Verbatim bocks parsed correctly.

        Verbatim blocks delimited by ";...;" should have the final newline
        stripped. Whitespace may be stripped from the end of the line but not
        the beginning.
        """
        mmcif_dict = MMCIF2Dict(io.StringIO(textwrap.dedent(u"""\
            data_verbatim_test
            _test_value
            ;First line
                Second line
            Third line
            ;
        """)))
        self.assertEqual(
            mmcif_dict["_test_value"],
            "First line\n    Second line\nThird line")

    def test_inline_comments(self):
        """Comments may begin outside of column 1 if preceded by whitespace."""
        mmcif_dict = MMCIF2Dict(io.StringIO(textwrap.dedent(u"""\
            data_verbatim_test
            _test_key_value_1 foo # Ignore this comment
            _test_key_value_2 foo#NotIgnored
            loop_
            _test_loop
            a b c d # Ignore this comment
            e f g

        """)))
        self.assertEqual(mmcif_dict["_test_key_value_1"], "foo")
        self.assertEqual(mmcif_dict["_test_key_value_2"], "foo#NotIgnored")
        self.assertEqual(mmcif_dict["_test_loop"], list("abcdefg"))

    def test_loop_keyword_case_insensitive(self):
        """Comments may begin outside of column 1."""
        test_data = u"""\
            data_verbatim_test
            _test_key_value foo # Ignore this comment
            loop_
            _test_loop
            a b c d # Ignore this comment
            e f g

        """
        mmcif_dict = MMCIF2Dict(io.StringIO(textwrap.dedent(test_data)))

        mmcif_dict2 = MMCIF2Dict(io.StringIO(textwrap.dedent(test_data.replace("loop_", "LOOP_"))))
        self.assertDictEqual(mmcif_dict,
                             mmcif_dict2)

        mmcif_dict2 = MMCIF2Dict(io.StringIO(textwrap.dedent(test_data.replace("loop_", "looP_"))))
        self.assertDictEqual(mmcif_dict,
                             mmcif_dict2)

        mmcif_dict2 = MMCIF2Dict(io.StringIO(textwrap.dedent(test_data.replace("_loop", "_LOOP"))))
        self.assertNotEqual(mmcif_dict,
                            mmcif_dict2)


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
