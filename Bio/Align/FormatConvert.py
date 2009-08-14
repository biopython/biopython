"""
For conversion between different formats for representing alignments (DEPRECATED).

This module is considered obsolete and has been deprecated.  Please use
Bio.AlignIO instead for reading and writing alignments in different file
formats.

classes:
FormatConverter
"""

# biopython
from Bio.Fasta.FastaAlign import FastaAlignment
from Bio.Clustalw import ClustalAlignment

import warnings
warnings.warn("Bio.Align.FormatConvert is deprecated. Please use Bio.AlignIO "
              "or the Alignment object's format method instead.",
              DeprecationWarning)

class FormatConverter:
    """Convert between different alignment representation formats.

    The basic idea behind the converter is that it takes a given format,
    converts it into the base Alignment class, and then can return an object
    in any other supported format with the info from the basal alignment.

    Supported formats are:
    o Clustal format (*.aln)
    o Fasta format (*.fasta)
    """
    def __init__(self, to_convert):
        """Initialize a converter with a given object.

        Arguments:
        o to_convert - The class which we are going to be converting.
        """
        self.orig_format = to_convert

        self._base_alphabet, self._base_records = self._get_base_info()

    def _get_base_info(self):
        """Retrieve all of the basal (ie Generic.Alignment) info.

        The idea is that this info is present in all of the classes and
        this is the information that will be retained in a conversion.
        Format specific information will be lost.
        """
        return self.orig_format._alphabet, self.orig_format._records

    def to_fasta(self):
        """Convert the current info into a FastaAlignment object.
        """
        return_format = FastaAlignment(self._base_alphabet)
        return_format._records = self._base_records

        return return_format

    def to_clustal(self):
        """Convert the current info into a ClustalAlignment object.
        """
        return_format = ClustalAlignment(self._base_alphabet)
        return_format._records = self._base_records

        return return_format
