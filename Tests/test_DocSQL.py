#!/usr/bin/env python
"""Test the Bio.DocSQL module
"""

from __future__ import print_function

import warnings
from Bio import BiopythonDeprecationWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonDeprecationWarning)
    import Bio.DocSQL

print("Skipping Bio.DocSQL doctests.")
# print("Running Bio.DocSQL doctests...")
# Bio.DocSQL._test()
# print("Bio.DocSQL doctests complete.")
