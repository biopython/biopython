"""This module is deprecated and will eventually be removed.

Bio.sequtils and Bio.SeqUtils duplicate a lot of code and Bio.SeqUtils
includes a directory where other utility code can be deposited. As a result
sequtils is being deprecated for now and will eventually be removed.

To update your code change:
    from Bio import sequtils
to:
    from Bio import SeqUtils
"""
import warnings
warnings.warn("Bio.sequtils is deprecated. Please use Bio.SeqUtils instead.",
              DeprecationWarning)
del warnings

# import the implementations of the code from SeqUtils
# temporary way to avoid duplication of code
from Bio.SeqUtils import *
