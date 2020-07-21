# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# Copyright 2012-2013 by Michiel JL de Hoon.  All rights reserved.
# Revisions copyright 2019 by Victor Lin.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tools for sequence motif analysis.

Bio.motifs contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.
It also includes functionality for parsing output from the AlignACE, MEME,
and MAST programs, as well as files in the TRANSFAC format.
"""

import warnings

from urllib.parse import urlencode
from urllib.request import urlopen, Request

from Bio import BiopythonDeprecationWarning


def create(instances, alphabet="ACGT"):
    """Create a Motif object."""
    instances = Instances(instances, alphabet)
    return Motif(instances=instances, alphabet=alphabet)


def parse(handle, fmt, strict=True):
    """Parse an output file from a motif finding program.

    Currently supported formats (case is ignored):
     - AlignAce:         AlignAce output file format
     - ClusterBuster:    Cluster Buster position frequency matrix format
     - XMS:              XMS matrix format
     - MEME:             MEME output file motif
     - MINIMAL:          MINIMAL MEME output file motif
     - MAST:             MAST output file motif
     - TRANSFAC:         TRANSFAC database file format
     - pfm-four-columns: Generic position-frequency matrix format with four columns. (cisbp, homer, hocomoco, neph, tiffin)
     - pfm-four-rows:    Generic position-frequency matrix format with four row. (scertf, yetfasco, hdpi, idmmpmm, flyfactor survey)
     - pfm:              JASPAR-style position-frequency matrix
     - jaspar:           JASPAR-style multiple PFM format
     - sites:            JASPAR-style sites file

    As files in the pfm and sites formats contain only a single motif,
    it is easier to use Bio.motifs.read() instead of Bio.motifs.parse()
    for those.

    For example:

    >>> from Bio import motifs
    >>> with open("motifs/alignace.out") as handle:
    ...     for m in motifs.parse(handle, "AlignAce"):
    ...         print(m.consensus)
    ...
    TCTACGATTGAG
    CTGCACCTAGCTACGAGTGAG
    GTGCCCTAAGCATACTAGGCG
    GCCACTAGCAGAGCAGGGGGC
    CGACTCAGAGGTT
    CCACGCTAAGAGAAGTGCCGGAG
    GCACGTCCCTGAGCA
    GTCCATCGCAAAGCGTGGGGC
    GAGATCAGAGGGCCG
    TGGACGCGGGG
    GACCAGAGCCTCGCATGGGGG
    AGCGCGCGTG
    GCCGGTTGCTGTTCATTAGG
    ACCGACGGCAGCTAAAAGGG
    GACGCCGGGGAT
    CGACTCGCGCTTACAAGG

    If strict is True (default), the parser will raise a ValueError if the
    file contents does not strictly comply with the specified file format.
    """
    fmt = fmt.lower()
    if fmt == "alignace":
        from Bio.motifs import alignace

        return alignace.read(handle)
    elif fmt == "meme":
        from Bio.motifs import meme

        return meme.read(handle)
    elif fmt == "minimal":
        from Bio.motifs import minimal

        return minimal.read(handle)
    elif fmt == "clusterbuster":
        from Bio.motifs import clusterbuster

        return clusterbuster.read(handle)
    elif fmt in ("pfm-four-columns", "pfm-four-rows"):
        from Bio.motifs import pfm

        return pfm.read(handle, fmt)
    elif fmt == "xms":
        from Bio.motifs import xms

        return xms.read(handle)
    elif fmt == "mast":
        from Bio.motifs import mast

        return mast.read(handle)
    elif fmt == "transfac":
        from Bio.motifs import transfac

        return transfac.read(handle, strict)
    elif fmt in ("pfm", "sites", "jaspar"):
        from Bio.motifs import jaspar

        return jaspar.read(handle, fmt)
    else:
        raise ValueError("Unknown format %s" % fmt)


def read(handle, fmt, strict=True):
    """Read a motif from a handle using the specified file-format.

    This supports the same formats as Bio.motifs.parse(), but
    only for files containing exactly one motif.  For example,
    reading a JASPAR-style pfm file:

    >>> from Bio import motifs
    >>> with open("motifs/SRF.pfm") as handle:
    ...     m = motifs.read(handle, "pfm")
    >>> m.consensus
    Seq('GCCCATATATGG')

    Or a single-motif MEME file,

    >>> from Bio import motifs
    >>> with open("motifs/meme.psp_test.classic.zoops.xml") as handle:
    ...     m = motifs.read(handle, "meme")
    >>> m.consensus
    Seq('GCTTATGTAA')

    If the handle contains no records, or more than one record,
    an exception is raised:

    >>> from Bio import motifs
    >>> with open("motifs/alignace.out") as handle:
    ...     motif = motifs.read(handle, "AlignAce")
    Traceback (most recent call last):
        ...
    ValueError: More than one motif found in handle

    If however you want the first motif from a file containing
    multiple motifs this function would raise an exception (as
    shown in the example above).  Instead use:

    >>> from Bio import motifs
    >>> with open("motifs/alignace.out") as handle:
    ...     record = motifs.parse(handle, "alignace")
    >>> motif = record[0]
    >>> motif.consensus
    Seq('TCTACGATTGAG')

    Use the Bio.motifs.parse(handle, fmt) function if you want
    to read multiple records from the handle.

    If strict is True (default), the parser will raise a ValueError if the
    file contents does not strictly comply with the specified file format.
    """
    fmt = fmt.lower()
    motifs = parse(handle, fmt, strict)
    if len(motifs) == 0:
        raise ValueError("No motifs found in handle")
    if len(motifs) > 1:
        raise ValueError("More than one motif found in handle")
    motif = motifs[0]
    return motif


class Instances(list):
    """Class containing a list of sequences that made the motifs."""

    def __init__(self, instances=None, alphabet="ACGT"):
        """Initialize the class."""
        from Bio.Seq import Seq

        if instances is None:
            instances = []
        self.length = None
        for instance in instances:
            if self.length is None:
                self.length = len(instance)
            elif self.length != len(instance):
                message = (
                    "All instances should have the same length (%d found, %d expected)"
                    % (len(instance), self.length)
                )
                raise ValueError(message)
        for instance in instances:
            if not isinstance(instance, Seq):
                sequence = str(instance)
                instance = Seq(sequence)
            self.append(instance)
        self.alphabet = alphabet

    def __str__(self):
        """Return a string containing the sequences of the motif."""
        text = ""
        for instance in self:
            text += str(instance) + "\n"
        return text

    def count(self):
        """Count nucleotides in a position."""
        counts = {}
        for letter in self.alphabet:
            counts[letter] = [0] * self.length
        for instance in self:
            for position, letter in enumerate(instance):
                counts[letter][position] += 1
        return counts

    def search(self, sequence):
        """Find positions of motifs in a given sequence.

        This is a generator function, returning found positions of motif
        instances in a given sequence.
        """
        for pos in range(0, len(sequence) - self.length + 1):
            for instance in self:
                if str(instance) == str(sequence[pos : pos + self.length]):
                    yield (pos, instance)
                    break  # no other instance will fit (we don't want to return multiple hits)

    def reverse_complement(self):
        """Compute reverse complement of sequences."""
        instances = Instances(alphabet=self.alphabet)
        instances.length = self.length
        for instance in self:
            instance = instance.reverse_complement()
            instances.append(instance)
        return instances


class Motif:
    """A class representing sequence motifs."""

    def __init__(self, alphabet="ACGT", instances=None, counts=None):
        """Initialize the class."""
        from . import matrix

        self.name = ""
        if counts is not None and instances is not None:
            raise Exception(
                ValueError, "Specify either instances or counts, don't specify both"
            )
        elif counts is not None:
            self.instances = None
            self.counts = matrix.FrequencyPositionMatrix(alphabet, counts)
            self.length = self.counts.length
        elif instances is not None:
            self.instances = instances
            alphabet = self.instances.alphabet
            counts = self.instances.count()
            self.counts = matrix.FrequencyPositionMatrix(alphabet, counts)
            self.length = self.counts.length
        else:
            self.counts = None
            self.instances = None
            self.length = None
        self.alphabet = alphabet
        self.pseudocounts = None
        self.background = None
        self.mask = None

    def __get_mask(self):
        return self.__mask

    def __set_mask(self, mask):
        if self.length is None:
            self.__mask = ()
        elif mask is None:
            self.__mask = (1,) * self.length
        elif len(mask) != self.length:
            raise ValueError(
                "The length (%d) of the mask is inconsistent with the length (%d) of the motif",
                (len(mask), self.length),
            )
        elif isinstance(mask, str):
            self.__mask = []
            for char in mask:
                if char == "*":
                    self.__mask.append(1)
                elif char == " ":
                    self.__mask.append(0)
                else:
                    raise ValueError(
                        "Mask should contain only '*' or ' ' and not a '%s'" % char
                    )
            self.__mask = tuple(self.__mask)
        else:
            self.__mask = tuple(int(bool(c)) for c in mask)

    mask = property(__get_mask, __set_mask)
    del __get_mask
    del __set_mask

    def __get_pseudocounts(self):
        return self._pseudocounts

    def __set_pseudocounts(self, value):
        self._pseudocounts = {}
        if isinstance(value, dict):
            self._pseudocounts = {letter: value[letter] for letter in self.alphabet}
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = dict.fromkeys(self.alphabet, value)

    pseudocounts = property(__get_pseudocounts, __set_pseudocounts)
    del __get_pseudocounts
    del __set_pseudocounts

    def __get_background(self):
        return self._background

    def __set_background(self, value):
        if isinstance(value, dict):
            self._background = {letter: value[letter] for letter in self.alphabet}
        elif value is None:
            self._background = dict.fromkeys(self.alphabet, 1.0)
        else:
            if sorted(self.alphabet) != ["A", "C", "G", "T"]:
                raise ValueError(
                    "Setting the background to a single value only works for DNA motifs"
                    " (in which case the value is interpreted as the GC content)"
                )
            self._background["A"] = (1.0 - value) / 2.0
            self._background["C"] = value / 2.0
            self._background["G"] = value / 2.0
            self._background["T"] = (1.0 - value) / 2.0
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total

    background = property(__get_background, __set_background)
    del __get_background
    del __set_background

    @property
    def pwm(self):
        """Compute position weight matrices."""
        return self.counts.normalize(self._pseudocounts)

    @property
    def pssm(self):
        """Compute position specific scoring matrices."""
        return self.pwm.log_odds(self._background)

    def __str__(self, masked=False):
        """Return string representation of a motif."""
        text = ""
        if self.instances is not None:
            text += str(self.instances)

        if masked:
            for i in range(self.length):
                if self.__mask[i]:
                    text += "*"
                else:
                    text += " "
            text += "\n"
        return text

    def __len__(self):
        """Return the length of a motif.

        Please use this method (i.e. invoke len(m)) instead of referring to m.length directly.
        """
        if self.length is None:
            return 0
        else:
            return self.length

    def reverse_complement(self):
        """Return the reverse complement of the motif as a new motif."""
        alphabet = self.alphabet
        if self.instances is not None:
            instances = self.instances.reverse_complement()
            res = Motif(alphabet=alphabet, instances=instances)
        else:  # has counts
            counts = {
                "A": self.counts["T"][::-1],
                "C": self.counts["G"][::-1],
                "G": self.counts["C"][::-1],
                "T": self.counts["A"][::-1],
            }
            res = Motif(alphabet=alphabet, counts=counts)
        res.__mask = self.__mask[::-1]
        res.background = {
            "A": self.background["T"],
            "C": self.background["G"],
            "G": self.background["C"],
            "T": self.background["A"],
        }
        res.pseudocounts = {
            "A": self.pseudocounts["T"],
            "C": self.pseudocounts["G"],
            "G": self.pseudocounts["C"],
            "T": self.pseudocounts["A"],
        }
        return res

    @property
    def consensus(self):
        """Return the consensus sequence."""
        return self.counts.consensus

    @property
    def anticonsensus(self):
        """Return the least probable pattern to be generated from this motif."""
        return self.counts.anticonsensus

    @property
    def degenerate_consensus(self):
        """Return the degenerate consensus sequence.

        Following the rules adapted from
        D. R. Cavener: "Comparison of the consensus sequence flanking
        translational start sites in Drosophila and vertebrates."
        Nucleic Acids Research 15(4): 1353-1361. (1987).

        The same rules are used by TRANSFAC.
        """
        return self.counts.degenerate_consensus

    def weblogo(self, fname, fmt="PNG", version="2.8.2", **kwds):
        """Download and save a weblogo using the Berkeley weblogo service.

        Requires an internet connection.

        The parameters from ``**kwds`` are passed directly to the weblogo server.

        Currently, this method uses WebLogo version 3.3.
        These are the arguments and their default values passed to
        WebLogo 3.3; see their website at http://weblogo.threeplusone.com
        for more information::

            'stack_width' : 'medium',
            'stacks_per_line' : '40',
            'alphabet' : 'alphabet_dna',
            'ignore_lower_case' : True,
            'unit_name' : "bits",
            'first_index' : '1',
            'logo_start' : '1',
            'logo_end': str(self.length),
            'composition' : "comp_auto",
            'percentCG' : '',
            'scale_width' : True,
            'show_errorbars' : True,
            'logo_title' : '',
            'logo_label' : '',
            'show_xaxis': True,
            'xaxis_label': '',
            'show_yaxis': True,
            'yaxis_label': '',
            'yaxis_scale': 'auto',
            'yaxis_tic_interval' : '1.0',
            'show_ends' : True,
            'show_fineprint' : True,
            'color_scheme': 'color_auto',
            'symbols0': '',
            'symbols1': '',
            'symbols2': '',
            'symbols3': '',
            'symbols4': '',
            'color0': '',
            'color1': '',
            'color2': '',
            'color3': '',
            'color4': '',

        """
        if set(self.alphabet) == set("ACDEFGHIKLMNPQRSTVWY"):
            alpha = "alphabet_protein"
        elif set(self.alphabet) == set("ACGU"):
            alpha = "alphabet_rna"
        elif set(self.alphabet) == set("ACGT"):
            alpha = "alphabet_dna"
        else:
            alpha = "auto"

        frequencies = format(self, "transfac")
        url = "http://weblogo.threeplusone.com/create.cgi"
        values = {
            "sequences": frequencies,
            "format": fmt.lower(),
            "stack_width": "medium",
            "stacks_per_line": "40",
            "alphabet": alpha,
            "ignore_lower_case": True,
            "unit_name": "bits",
            "first_index": "1",
            "logo_start": "1",
            "logo_end": str(self.length),
            "composition": "comp_auto",
            "percentCG": "",
            "scale_width": True,
            "show_errorbars": True,
            "logo_title": "",
            "logo_label": "",
            "show_xaxis": True,
            "xaxis_label": "",
            "show_yaxis": True,
            "yaxis_label": "",
            "yaxis_scale": "auto",
            "yaxis_tic_interval": "1.0",
            "show_ends": True,
            "show_fineprint": True,
            "color_scheme": "color_auto",
            "symbols0": "",
            "symbols1": "",
            "symbols2": "",
            "symbols3": "",
            "symbols4": "",
            "color0": "",
            "color1": "",
            "color2": "",
            "color3": "",
            "color4": "",
        }

        values.update({k: "" if v is False else str(v) for k, v in kwds.items()})
        data = urlencode(values).encode("utf-8")
        req = Request(url, data)
        response = urlopen(req)
        with open(fname, "wb") as f:
            im = response.read()
            f.write(im)

    def __format__(self, format_spec):
        """Return a string representation of the Motif in the given format.

        Currently supported formats:
         - clusterbuster: Cluster Buster position frequency matrix format
         - pfm : JASPAR single Position Frequency Matrix
         - jaspar : JASPAR multiple Position Frequency Matrix
         - transfac : TRANSFAC like files

        """
        if format_spec in ("pfm", "jaspar"):
            from Bio.motifs import jaspar

            motifs = [self]
            return jaspar.write(motifs, format_spec)
        elif format_spec == "transfac":
            from Bio.motifs import transfac

            motifs = [self]
            return transfac.write(motifs)
        elif format_spec == "clusterbuster":
            from Bio.motifs import clusterbuster

            motifs = [self]
            return clusterbuster.write(motifs)
        else:
            raise ValueError("Unknown format type %s" % format_spec)

    def format(self, format_spec):
        """Return a string representation of the Motif in the given format [DEPRECATED].

        This method is deprecated; instead of motif.format(format_spec),
        please use format(motif, format_spec).
        """
        warnings.warn(
            "Motif.format has been deprecated, and we intend to remove it in a "
            "future release of Biopython. Instead of motif.format(format_spec), "
            "please use format(motif, format_spec).",
            BiopythonDeprecationWarning,
        )
        return self.__format__(format_spec)


def write(motifs, fmt):
    """Return a string representation of motifs in the given format.

    Currently supported formats (case is ignored):
     - clusterbuster: Cluster Buster position frequency matrix format
     - pfm : JASPAR simple single Position Frequency Matrix
     - jaspar : JASPAR multiple PFM format
     - transfac : TRANSFAC like files

    """
    fmt = fmt.lower()
    if fmt in ("pfm", "jaspar"):
        from Bio.motifs import jaspar

        return jaspar.write(motifs, fmt)
    elif fmt == "transfac":
        from Bio.motifs import transfac

        return transfac.write(motifs)
    elif fmt == "clusterbuster":
        from Bio.motifs import clusterbuster

        return clusterbuster.write(motifs)
    else:
        raise ValueError("Unknown format type %s" % fmt)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
