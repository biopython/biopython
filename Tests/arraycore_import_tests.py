# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Isolated tests for ArrayCore extension initialization."""

import gc
import importlib.machinery
import importlib.util
import sys
import types
import unittest
import weakref


class _InvalidArray:
    pass


def _load_arraycore(extension):
    module_name = "Bio.Align.substitution_matrices._arraycore"
    loader = importlib.machinery.ExtensionFileLoader(module_name, extension)
    spec = importlib.util.spec_from_file_location(module_name, extension, loader=loader)
    if spec is None:
        raise AssertionError("could not create an arraycore import specification")
    return importlib.util.module_from_spec(spec)


def _check_successful_arraycore_import(extension):
    import numpy

    base_references = sys.getrefcount(numpy.ndarray)
    module = _load_arraycore(extension)
    reference_delta = sys.getrefcount(numpy.ndarray) - base_references
    if reference_delta not in (0, 2):
        raise AssertionError(
            "successful import retained an unexpected number of "
            f"numpy.ndarray references: {reference_delta}"
        )
    if module.Array.__base__ is not numpy.ndarray:
        raise AssertionError("Array does not retain numpy.ndarray as its base")
    del module
    gc.collect()


def _check_retried_arraycore_import(extension):
    fake_numpy = types.ModuleType("numpy")
    fake_numpy.ndarray = _InvalidArray
    previous_numpy = sys.modules.get("numpy")
    sys.modules["numpy"] = fake_numpy
    try:
        try:
            _load_arraycore(extension)
        except Exception:
            pass
        else:
            raise AssertionError("arraycore accepted a dynamically allocated base type")

        fake_numpy.ndarray = object
        try:
            _load_arraycore(extension)
        except RuntimeError as exception:
            if str(exception) != "Array type initialization previously failed":
                raise
        else:
            raise AssertionError("arraycore retried a poisoned static type")
    finally:
        if previous_numpy is None:
            del sys.modules["numpy"]
        else:
            sys.modules["numpy"] = previous_numpy


def _check_failed_arraycore_import(mode, extension):
    fake_numpy = types.ModuleType("numpy")
    invalid_array = _InvalidArray() if mode == "invalid" else None
    if invalid_array is not None:
        fake_numpy.ndarray = invalid_array

    numpy_reference = weakref.ref(fake_numpy)
    invalid_reference = (
        weakref.ref(invalid_array) if invalid_array is not None else None
    )
    previous_numpy = sys.modules.get("numpy")
    sys.modules["numpy"] = fake_numpy
    try:
        try:
            _load_arraycore(extension)
        except RuntimeError:
            pass
        else:
            raise AssertionError("arraycore accepted an invalid numpy.ndarray")
    finally:
        if previous_numpy is None:
            del sys.modules["numpy"]
        else:
            sys.modules["numpy"] = previous_numpy

    del fake_numpy
    del invalid_array
    gc.collect()
    if numpy_reference() is not None:
        raise AssertionError("failed import retained the NumPy module")
    if invalid_reference is not None and invalid_reference() is not None:
        raise AssertionError("failed import retained the invalid base object")


class ArrayCoreImportTests(unittest.TestCase):
    def test_successful_import(self):
        _check_successful_arraycore_import(sys.argv[2])

    def test_missing_ndarray(self):
        _check_failed_arraycore_import("missing", sys.argv[2])

    def test_invalid_ndarray(self):
        _check_failed_arraycore_import("invalid", sys.argv[2])

    def test_poisoned_retry(self):
        _check_retried_arraycore_import(sys.argv[2])


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0], sys.argv[1]])
