# Copyright 2014-2016 Marco Galardini.  All rights reserved.
# Adapted from test_Mymodule.py by Jeff Chang
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Tests for the Bio.phenotype module."""

try:
    import numpy

    del numpy
except ImportError:
    from Bio import MissingExternalDependencyError

    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.phenotype."
    ) from None

import json
import unittest

from io import StringIO

from Bio import BiopythonExperimentalWarning

import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio import phenotype

# Example plate files
SMALL_JSON_PLATE = "phenotype/SmallPlate.json"
SMALL_JSON_PLATE_2 = "phenotype/SmallPlate_2.json"
JSON_PLATE = "phenotype/Plate.json"
JSON_PLATE_2 = "phenotype/Plate_2.json"
JSON_PLATE_3 = "phenotype/Plate_3.json"
JSON_BAD = "phenotype/BadPlate.json"

SMALL_CSV_PLATES = "phenotype/SmallPlates.csv"
CSV_PLATES = "phenotype/Plates.csv"


class TestPhenoMicro(unittest.TestCase):
    """Tests for phenotype support."""

    def test_phenotype_IO_errors(self):
        """Test bad arguments to phenotype IO methods."""
        self.assertRaises(ValueError, phenotype.read, CSV_PLATES, "pm-csv")
        self.assertRaises(ValueError, phenotype.read, CSV_PLATES, "pm-json")
        self.assertRaises(ValueError, phenotype.read, CSV_PLATES, "pm-noformat")
        self.assertRaises(ValueError, phenotype.read, CSV_PLATES, "PM-CSV")
        self.assertRaises(TypeError, phenotype.read, CSV_PLATES, 1)
        self.assertRaises(KeyError, phenotype.read, JSON_BAD, "pm-json")

    def test_phenotype_IO(self):
        """Test basic functionalities of phenotype IO methods."""
        p1 = phenotype.read(SMALL_JSON_PLATE, "pm-json")
        p2 = next(phenotype.parse(SMALL_CSV_PLATES, "pm-csv"))

        handle = StringIO()

        c = phenotype.write([p1, p2], handle, "pm-json")
        self.assertEqual(c, 2)

        handle.flush()
        handle.seek(0)
        # Now ready to read back from the handle...
        try:
            records = list(phenotype.parse(handle, "pm-json"))
        except ValueError as e:
            # This is BAD.  We can't read our own output.
            # I want to see the output when called from the test harness,
            # run_tests.py (which can be funny about new lines on Windows)
            handle.seek(0)
            self.fail(f"{e}\n\n{handle.read()!r}\n\n{records!r}")

        self.assertEqual(p1, records[0])

        handle.close()
        handle = StringIO()
        self.assertRaises(TypeError, phenotype.write, p1, handle, 1)
        self.assertRaises(ValueError, phenotype.write, p1, handle, "PM-JSON")
        self.assertRaises(ValueError, phenotype.write, p1, handle, "pm-csv")
        handle.close()

    def test_PlateRecord_errors(self):
        """Test bad arguments with PlateRecord objects."""
        self.assertRaises(
            ValueError, phenotype.phen_micro.PlateRecord, "test", [1, 2, 3]
        )
        self.assertRaises(TypeError, phenotype.phen_micro.PlateRecord, "test", 1)

    def test_PlateRecord(self):
        """Test basic functionalities of PlateRecord objects."""
        with open(SMALL_JSON_PLATE) as handle:
            j = json.load(handle)

        p = phenotype.phen_micro.PlateRecord(j["csv_data"]["Plate Type"])

        times = j["measurements"]["Hour"]
        for k in j["measurements"]:
            if k == "Hour":
                continue
            p[k] = phenotype.phen_micro.WellRecord(
                k,
                signals={times[i]: j["measurements"][k][i] for i in range(len(times))},
            )

        del j["measurements"]
        p.qualifiers = j

        self.assertEqual(p.id, "PM01")
        self.assertEqual(len(p), 24)
        self.assertEqual(p.qualifiers, j)
        self.assertRaises(ValueError, p._is_well, "a")
        self.assertEqual(p["A01"].id, "A01")
        self.assertRaises(KeyError, p.__getitem__, "test")
        self.assertEqual(len(p[1]), 12)
        self.assertEqual(len(p[1:2:2]), 12)
        self.assertEqual(p[1, 2], p["B03"])
        self.assertEqual(len(p[:, 1]), 2)
        self.assertEqual(len(p[:, 1:4:2]), 4)
        self.assertRaises(TypeError, p.__getitem__, 1, 2, 3)
        self.assertRaises(IndexError, p.__getitem__, 13)
        self.assertRaises(ValueError, p.__setitem__, "A02", p["A01"])
        self.assertRaises(ValueError, p.__setitem__, "A02", "a")
        p["A02"] = p["A02"]
        for w in p:
            pass
        self.assertIn("A01", p)
        self.assertNotIn("test", p)
        self.assertRaises(ValueError, next, p.get_row("test"))
        self.assertEqual(next(p.get_row("A")), p["A01"])
        self.assertRaises(ValueError, next, p.get_column("test"))
        self.assertEqual(next(p.get_column("12")), p["A12"])
        self.assertEqual(next(p.get_column("1")), p["A01"])
        self.assertRaises(ValueError, p.subtract_control, "A121")
        self.assertRaises(ValueError, p.subtract_control, wells=["A121"])
        p2 = p.subtract_control()
        self.assertEqual(p2.id, p.id)
        self.assertEqual(p2["A02"], p["A02"] - p["A01"])
        self.assertEqual(
            repr(p),
            "PlateRecord('WellRecord['A01'], WellRecord['A02'], "
            "WellRecord['A03'], ..., WellRecord['B12']')",
        )
        self.assertEqual(
            str(p),
            "Plate ID: PM01\nWell: 24\nRows: 2\nColumns: 12\n"
            "PlateRecord('WellRecord['A01'], WellRecord['A02'], "
            "WellRecord['A03'], ..., WellRecord['B12']')",
        )

        with open(SMALL_JSON_PLATE_2) as handle:
            j = json.load(handle)

        p1 = phenotype.phen_micro.PlateRecord(j["csv_data"]["Plate Type"])

        times = j["measurements"]["Hour"]
        for k in j["measurements"]:
            if k == "Hour":
                continue
            p1[k] = phenotype.phen_micro.WellRecord(
                k,
                signals={times[i]: j["measurements"][k][i] for i in range(len(times))},
            )

        del j["measurements"]
        p1.qualifiers = j

        self.assertRaises(TypeError, p.__add__, "a")
        self.assertRaises(TypeError, p.__sub__, "a")

        p3 = p + p1
        self.assertEqual(p3["A02"], p["A02"] + p1["A02"])

        p3 = p - p1
        self.assertEqual(p3["A02"], p["A02"] - p1["A02"])

        del p["A02"]
        self.assertRaises(ValueError, p.__add__, p1)
        self.assertRaises(ValueError, p.__sub__, p1)

    def test_bad_fit_args(self):
        """Test error handling of the fit method."""
        with open(JSON_PLATE) as handle:
            p = json.load(handle)

        times = p["measurements"]["Hour"]
        w = phenotype.phen_micro.WellRecord(
            "A10",
            signals={times[i]: p["measurements"]["A10"][i] for i in range(len(times))},
        )

        self.assertRaises(ValueError, w.fit, "wibble")
        self.assertRaises(ValueError, w.fit, ["wibble"])
        self.assertRaises(ValueError, w.fit, ("logistic", "wibble"))
        self.assertRaises(ValueError, w.fit, ("wibble", "logistic"))
        self.assertRaises(ValueError, w.fit, "logistic")  # should be a list/tuple!

    def test_WellRecord(self):
        """Test basic functionalities of WellRecord objects."""
        with open(JSON_PLATE) as handle:
            p = json.load(handle)

        times = p["measurements"]["Hour"]
        w = phenotype.phen_micro.WellRecord(
            "A10",
            signals={times[i]: p["measurements"]["A10"][i] for i in range(len(times))},
        )

        w1 = phenotype.phen_micro.WellRecord(
            "H12",
            signals={times[i]: p["measurements"]["H12"][i] for i in range(len(times))},
        )

        # self.assertIsInstance(w.plate,
        #                       phenotype.phen_micro.PlateRecord)
        self.assertIsInstance(w.plate, phenotype.phen_micro.PlateRecord)
        self.assertEqual(w.id, "A10")
        self.assertEqual(len(w), len(times))
        self.assertEqual(len(w), 384)
        self.assertEqual(max(w), (95.75, 217.0))
        self.assertEqual(min(w), (0.0, 37.0))
        self.assertEqual(max(w, key=lambda x: x[1]), (16.75, 313.0))  # noqa: E731
        self.assertEqual(min(w, key=lambda x: x[1]), (0.25, 29.0))  # noqa: E731
        self.assertEqual(len(w[:]), 96)
        self.assertEqual(w[1], 29.0)
        self.assertEqual(w[12], 272.0)
        self.assertEqual(w[1:5], [29.0, 35.0, 39.0, 43.0])
        self.assertRaises(ValueError, w.__getitem__, "a")
        self.assertAlmostEqual(w[1:2:0.25][0], 29.0)
        self.assertAlmostEqual(w[1.3567], 33.7196)
        self.assertEqual(w.get_raw()[0], (0.0, 37.0))
        self.assertEqual(w.get_raw()[-1], (95.75, 217.0))
        self.assertEqual(w.get_times()[0], 0.0)
        self.assertEqual(w.get_times()[-1], 95.75)
        self.assertEqual(w.get_signals()[0], 37.0)
        self.assertEqual(w.get_signals()[-1], 217.0)
        self.assertEqual(
            repr(w),
            "WellRecord('(0.0, 37.0), (0.25, 29.0), (0.5, 32.0),"
            " (0.75, 30.0), (1.0, 29.0), ..., (95.75, 217.0)')",
        )
        self.assertEqual(
            str(w),
            "Well ID: A10\nTime points: 384\nMinum signal 0.25 at time 29.00\n"
            "Maximum signal 16.75 at time 313.00\n"
            "WellRecord('(0.0, 37.0), (0.25, 29.0), (0.5, 32.0), (0.75, 30.0), "
            "(1.0, 29.0), ..., (95.75, 217.0)')",
        )

        w.fit(None)
        self.assertIsNone(w.area)
        self.assertIsNone(w.model)
        self.assertIsNone(w.lag)
        self.assertIsNone(w.plateau)
        self.assertIsNone(w.slope)
        self.assertIsNone(w.v)
        self.assertIsNone(w.y0)
        self.assertEqual(w.max, 313.0)
        self.assertEqual(w.min, 29.0)
        self.assertEqual(w.average_height, 217.82552083333334)

        self.assertRaises(TypeError, w.__add__, "a")

        w2 = w + w1
        self.assertEqual(w2.id, "A10")
        self.assertEqual(len(w2), len(times))
        self.assertEqual(len(w2), 384)
        self.assertEqual(max(w2), (95.75, 327.0))
        self.assertEqual(min(w2), (0.0, 63.0))
        self.assertEqual(max(w2, key=lambda x: x[1]), (18.25, 357.0))  # noqa: E731
        self.assertEqual(min(w2, key=lambda x: x[1]), (0.25, 55.0))  # noqa: E731
        self.assertEqual(w2[1], 71.0)
        self.assertEqual(w2[12], 316.0)
        self.assertEqual(w2[1:5], [71.0, 88.0, 94.0, 94.0])
        self.assertAlmostEqual(w2[1:2:0.25][0], 71.0)
        self.assertAlmostEqual(w2[1.3567], 77.7196)
        self.assertEqual(w2.get_raw()[0], (0.0, 63.0))
        self.assertEqual(w2.get_raw()[-1], (95.75, 327.0))
        self.assertEqual(w2.get_times()[0], 0.0)
        self.assertEqual(w2.get_times()[-1], 95.75)
        self.assertEqual(w2.get_signals()[0], 63.0)
        self.assertEqual(w2.get_signals()[-1], 327.0)

        self.assertRaises(TypeError, w.__sub__, "a")

        w2 = w - w1
        self.assertEqual(w2.id, "A10")
        self.assertEqual(len(w2), len(times))
        self.assertEqual(len(w2), 384)
        self.assertEqual(max(w2), (95.75, 107.0))
        self.assertEqual(min(w2), (0.0, 11.0))
        self.assertEqual(max(w2, key=lambda x: x[1]), (15.75, 274.0))  # noqa: E731
        self.assertEqual(min(w2, key=lambda x: x[1]), (3.25, -20.0))  # noqa: E731
        self.assertEqual(w2[1], -13.0)
        self.assertEqual(w2[12], 228.0)
        self.assertEqual(w2[1:5], [-13.0, -18.0, -16.0, -8.0])
        self.assertAlmostEqual(w2[1:2:0.25][0], -13.0)
        self.assertAlmostEqual(w2[1.3567], -10.2804)
        self.assertEqual(w2.get_raw()[0], (0.0, 11.0))
        self.assertEqual(w2.get_raw()[-1], (95.75, 107.0))
        self.assertEqual(w2.get_times()[0], 0.0)
        self.assertEqual(w2.get_times()[-1], 95.75)
        self.assertEqual(w2.get_signals()[0], 11.0)
        self.assertEqual(w2.get_signals()[-1], 107.0)

        w[1] = 1

    def test_JsonIterator(self):
        """Test basic functionalities of JsonIterator file parser."""
        # Parse file content big enough to trigger issue #3783
        handle = StringIO(
            '{"csv_data": {"Plate Type": "PM-999"}, "measurements": {"Hour": 9}}'
        )
        with self.assertWarnsRegex(UserWarning, "PM-999"):
            for w in phenotype.phen_micro.JsonIterator(handle):
                self.assertEqual(w.id, "PM999")

    def test_CsvIterator(self):
        """Test basic functionalities of CsvIterator file parser."""
        # Parse file content big enough to trigger issue #3783
        handle = StringIO('"Data File",3\n"Plate Type",PM-33\n')
        with self.assertWarnsRegex(UserWarning, "PM-33"):
            for w in phenotype.phen_micro.CsvIterator(handle):
                self.assertEqual(w.id, "PM33")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
