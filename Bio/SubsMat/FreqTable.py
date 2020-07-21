# Copyright 2000 by Iddo Friedberg idoerg@cc.huji.ac.il
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

r"""A class to handle frequency tables or letter count files.

Example files for a DNA alphabet:

A count file (whitespace separated)::

 A  50
 C  37
 G  23
 T  58

The same info as a frequency file::

 A 0.2976
 C 0.2202
 G 0.1369
 T 0.3452

Functions:
  :read_count(f): read a count file from stream f. Then convert to
                  frequencies.
  :read_freq(f): read a frequency data file from stream f. Of course, we then
                 don't have the counts, but it is usually the letter frequencies
                 which are interesting.

Methods:
  (all internal)

Attributes:
  :alphabet: The letters you are using as indices into the table.
  :data: Frequency dictionary.
  :count: Count dictionary. Empty if no counts are provided.

Example of use:
    >>> import io
    >>> from Bio.SubsMat import FreqTable
    >>> f_count = io.StringIO(u"A  50\nC  37\nG  23\nT  58")
    >>> ftab = FreqTable.read_count(f_count)
    >>> for nb in sorted(ftab):
    ...     print("%s %0.4f" %(nb, ftab[nb]))
    ...
    A 0.2976
    C 0.2202
    G 0.1369
    T 0.3452

"""


COUNT = 1
FREQ = 2


class FreqTable(dict):
    """Define class to handle frequency tables or letter count files."""

    def _freq_from_count(self):
        """Calculate frequency from count values (PRIVATE)."""
        total = float(sum(self.count.values()))
        for i, v in self.count.items():
            self[i] = v / total

    def _alphabet_from_input(self):
        """Order the alphabet (PRIVATE)."""
        s = ""
        for i in sorted(self):
            s += i
        return s

    def __init__(self, in_dict, dict_type, alphabet=None):
        """Initialize the class."""
        self.alphabet = alphabet
        if dict_type == COUNT:
            self.count = in_dict
            self._freq_from_count()
        elif dict_type == FREQ:
            self.count = {}
            self.update(in_dict)
        else:
            raise ValueError("bad dict_type")
        if not alphabet:
            self.alphabet = self._alphabet_from_input()


def read_count(f):
    """Read a count file f and load values to the Frequency Table."""
    count = {}
    for line in f:
        key, value = line.strip().split()
        count[key] = int(value)
    return FreqTable(count, COUNT)


def read_freq(f):
    """Read a frequency data file f and load values to the Frequency Table."""
    freq_dict = {}
    for line in f:
        key, value = line.strip().split()
        freq_dict[key] = float(value)
    return FreqTable(freq_dict, FREQ)
