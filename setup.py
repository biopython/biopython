#!/usr/bin/env python
"""Freely available tools for computational molecular biology."""

from setuptools import setup

import sys
import warnings

_DEPRECATION_MESSAGE = (
    "Invoking setup.py is deprecated and will be removed in a future release of Biopython.\n"
    "Please use `pip install` or `python -m build` instead of `python setup.py` commands.\n"
    "For further information see https://packaging.python.org/en/latest/discussions/setup-py-deprecated/\n"
)

setup()

warnings.simplefilter("always", DeprecationWarning)
warnings.warn(_DEPRECATION_MESSAGE, category=DeprecationWarning, stacklevel=1)
print("DEPRECATION WARNING: " + _DEPRECATION_MESSAGE, file=sys.stderr)
