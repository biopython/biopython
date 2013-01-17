# Copyright 2003-2009 by Bartek Wilczynski.  All rights reserved.
# Copyright 2012-2013 by Michiel JL de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Module containing different tools for sequence motif analysis.

It contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.
It also includes functionality for parsing AlignACE and MEME programs.
"""

import warnings
from Bio import BiopythonExperimentalWarning


def create(instances, alphabet=None):
    from Bio.Alphabet import IUPAC
    from Bio.Seq import Seq
    for instance in instances:
        try:
            a = instance.alphabet
        except AttributeError:
            # The instance is a plain string
            continue
        if alphabet is None:
            alphabet = a
        elif alphabet != a:
            raise ValueError("Alphabets are inconsistent")
    if alphabet is None or alphabet.letters is None:
        # If we didn't get a meaningful alphabet from the instances,
        # assume it is DNA.
        alphabet = IUPAC.unambiguous_dna
    seqs = []
    for instance in instances:
        seq = Seq(str(instance), alphabet=alphabet)
        seqs.append(seq)
    return NewMotif(instances=seqs, alphabet=alphabet)


def parse(handle, format):
    """Parses an output file of motif finding programs.

    Currently supported formats:
     - AlignAce:      AlignAce output file format
     - MEME:          MEME output file motif
     - TRANSFAC:      TRANSFAC database file format
     - pfm:           JASPAR-style position-frequency matrix
     - sites:         JASPAR-style sites file
     - jaspar-pfm:    JASPAR-style position-frequency matrix [DEPRECATED]
     - jaspar-sites:  JASPAR-style sites file [DEPRECATED]
    As files in the pfm and sites formats contain only a single motif,
    it is easier to use Bio.Motif.read() instead of Bio.Motif.parse()
    for those.

    For example:

    >>> from Bio import Motif
    >>> for motif in Motif.parse(open("Motif/alignace.out"),"AlignAce"):
    ...     print motif.consensus()
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
    >>> for motif in Motif.parse(open("Motif/alignace.out"),"alignace"):
    ...     print motif.consensus
    TCTACGATTGAG
    CTGCAGCTAGCTACGAGTGAG
    GTGCTCTAAGCATAGTAGGCG
    GCCACTAGCAGAGCAGGGGGC
    CGACTCAGAGGTT
    CCACGCTAAGAGAGGTGCCGGAG
    GCGCGTCGCTGAGCA
    GTCCATCGCAAAGCGTGGGGC
    GGGATCAGAGGGCCG
    TGGAGGCGGGG
    GACCAGAGCTTCGCATGGGGG
    GGCGTGCGTG
    GCTGGTTGCTGTTCATTAGG
    GCCGGCGGCAGCTAAAAGGG
    GAGGCCGGGGAT
    CGACTCGTGCTTAGAAGG
    """
    if format=="AlignAce":
        # Old Motif code
        from Bio.Motif.Parsers import AlignAce
        record = AlignAce.read(handle)
        return iter(record.motifs)
    elif format=="alignace":
        # Old Motif code
        from Bio.Motif import AlignAce
        record = AlignAce.read(handle)
        return record
    elif format=="MEME":
        from Bio.Motif.Parsers import MEME
        record = MEME.read(handle)
        return iter(record.motifs)
    elif format=="meme":
        from Bio.Motif import MEME
        record = MEME.read(handle)
        return record
    elif format=="TRANSFAC":
        from Bio.Motif import TRANSFAC
        record = TRANSFAC.read(handle)
        return record
    elif format in ('pfm', 'sites'):
        from Bio.Motif import Jaspar
        motif = Jaspar.read(handle, format)
    elif format=="jaspar-pfm":
        motif = OldMotif()._from_jaspar_pfm(handle)
        return iter([motif])
    elif format=="jaspar-sites":
        motif = OldMotif()._from_jaspar_sites(handle)
        return iter([motif])
    else:
        raise ValueError("Unknown format %s" % format)
    # Treat the single-motif formats
    motifs = [motif]
    return motifs


def read(handle, format):
    """Reads a motif from a handle using a specified file-format.

    This supports the same formats as Bio.Motif.parse(), but
    only for files containing exactly one motif.  For example,
    reading a JASPAR-style pfm file:

    (old Motif code):
    >>> from Bio import Motif
    >>> motif = Motif.read(open("Motif/SRF.pfm"), "jaspar-pfm")
    >>> motif.consensus()
    Seq('GCCCATATATGG', IUPACUnambiguousDNA())

    (new Motif code):
    >>> from Bio import Motif
    >>> motif = Motif.read(open("Motif/SRF.pfm"), "pfm")
    >>> motif.consensus
    Seq('GCCCATATATGG', IUPACUnambiguousDNA())

    Or a single-motif MEME file,

    (old Motif code)
    >>> from Bio import Motif
    >>> motif =  Motif.read(open("Motif/meme.out"),"MEME")
    >>> motif.consensus()
    Seq('CTCAATCGTA', IUPACUnambiguousDNA())

    (new Motif code)
    >>> from Bio import Motif
    >>> motif =  Motif.read(open("Motif/meme.out"),"meme")
    >>> motif.consensus
    Seq('CTCAATCGTA', IUPACUnambiguousDNA())

    If the handle contains no records, or more than one record,
    an exception is raised:

    >>> from Bio import Motif
    >>> motif = Motif.read(open("Motif/alignace.out"),"AlignAce")
    Traceback (most recent call last):
        ...
    ValueError: More than one motif found in handle

    If however you want the first motif from a file containing
    multiple motifs this function would raise an exception (as
    shown in the example above).  Instead use:

    (old Motif code)
    >>> from Bio import Motif
    >>> motif = Motif.parse(open("Motif/alignace.out"),"AlignAce").next()
    >>> motif.consensus()
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

    (new Motif code)
    >>> from Bio import Motif
    >>> motifs = Motif.parse(open("Motif/alignace.out"),"alignace")
    >>> motif = motifs[0]
    >>> motif.consensus
    Seq('TCTACGATTGAG', IUPACUnambiguousDNA())

    Use the Bio.Motif.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    if format in ("AlignAce", "MEME", "jaspar-pfm", "jaspar-sites"):
        iterator = parse(handle, format)
        try:
            first = iterator.next()
        except StopIteration:
            first = None
        if first is None:
            raise ValueError("No motifs found in handle")
        try:
            second = iterator.next()
        except StopIteration:
            second = None
        if second is not None:
            raise ValueError("More than one motif found in handle")
        return first
    motifs = parse(handle, format)
    if len(motifs)==0:
        raise ValueError("No motifs found in handle")
    if len(motifs) > 1:
        raise ValueError("More than one motif found in handle")
    motif = motifs[0]
    return motif


class Motif(object):
    """
    A class representing sequence motifs.
    """
    def __init__(self, alphabet=None, instances=None, counts=None):
        import Matrix
        from Bio.Alphabet import IUPAC
        self.length=None
        self.name=""
        if counts is not None and instances is not None:
            raise Exception(ValueError,
                "Specify either instances or counts, don't specify both")
        elif counts is not None:
            warnings.warn("This is experimental code, and may change in future versions", BiopythonExperimentalWarning)
            if alphabet is None:
                alphabet = IUPAC.unambiguous_dna
            for letter in counts:
                length = len(counts[letter])
                if self.length is None:
                    self.length = length
                elif self.length!=length:
                    raise Exception("counts matrix has inconsistent lengths")
            self.instances = None
            self.counts = Matrix.FrequencyPositionMatrix(alphabet, counts)
        elif instances is not None:
            warnings.warn("This is experimental code, and may change in future versions", BiopythonExperimentalWarning)
            self.instances = []
            for instance in instances:
                if alphabet != instance.alphabet:
                    raise ValueError("Alphabets are inconsistent")
                if self.length is None:
                    self.length = len(instance)
                elif self.length != len(instance):
                    message = "All instances should have the same length (%d found, %d expected)" % (len(instance), self.length)
                    raise ValueError(message)
                self.instances.append(instance)
            counts = {}
            for letter in alphabet.letters:
                counts[letter] = [0] * self.length
            for instance in self.instances:
                for position, letter in enumerate(instance):
                    counts[letter][position] += 1
            self.counts = Matrix.FrequencyPositionMatrix(alphabet, counts)
        else:
            self.counts = None
            self.instances = None
            if alphabet is None:
                alphabet = IUPAC.unambiguous_dna
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
        elif len(mask)!=self.length:
            raise ValueError("The length (%d) of the mask is inconsistent with the length (%d) of the motif", (len(mask), self.length))
        elif isinstance(mask, str):
            self.__mask=[]
            for char in mask:
                if char=="*":
                    self.__mask.append(1)
                elif char==" ":
                    self.__mask.append(0)
                else:
                    raise ValueError("Mask should contain only '*' or ' ' and not a '%s'"%char)
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
            self._pseudocounts = dict((letter, value[letter]) for letter in self.alphabet.letters)
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = dict.fromkeys(self.alphabet.letters, value)

    pseudocounts = property(__get_pseudocounts, __set_pseudocounts)
    del __get_pseudocounts
    del __set_pseudocounts

    def __get_background(self):
        return self._background

    def __set_background(self, value):
        if isinstance(value, dict):
            self._background = dict((letter, value[letter]) for letter in self.alphabet.letters)
        elif value is None:
            self._background = dict.fromkeys(self.alphabet.letters, 1.0)
        else:
            if sorted(self.alphabet.letters)!=["A", "C", "G", "T"]:
                raise Exception("Setting the background to a single value only works for DNA motifs (in which case the value is interpreted as the GC content")
            self._background['A'] = (1.0-value)/2.0
            self._background['C'] = value/2.0
            self._background['G'] = value/2.0
            self._background['T'] = (1.0-value)/2.0
        total = sum(self._background.values())
        for letter in self.alphabet.letters:
            self._background[letter] /= total

    background = property(__get_background, __set_background)
    del __get_background
    del __set_background

    @property
    def pwm(self):
        return self.counts.normalize(self._pseudocounts)

    @property
    def pssm(self):
        return self.pwm.log_odds(self._background)

    def search_instances(self,sequence):
        """
        a generator function, returning found positions of instances of the motif in a given sequence
        """
        if self.instances is None:
            raise ValueError("This motif has no instances")
        for pos in xrange(0,len(sequence)-self.length+1):
            for instance in self.instances:
                if str(instance) == str(sequence[pos:pos+self.length]):
                    yield(pos,instance)
                    break # no other instance will fit (we don't want to return multiple hits)

    def __str__(self,masked=False):
        """ string representation of a motif.
        """
        string = ""
        if self.instances is not None:
            for inst in self.instances:
                string += str(inst) + "\n"

        if masked:
            for i in xrange(self.length):
                if self.__mask[i]:
                    string += "*"
                else:
                    string += " "
            string += "\n"
        return string

    def __len__(self):
        """return the length of a motif

        Please use this method (i.e. invoke len(m)) instead of referring to m.length directly.
        """
        if self.length is None:
            return 0
        else:
            return self.length

    def reverse_complement(self):
        """
        Gives the reverse complement of the motif
        """
        alphabet = self.alphabet
        if self.instances is not None:
            instances = []
            for instance in self.instances:
                instance = instance.reverse_complement()
                instances.append(instance)
            res = NewMotif(instances=instances, alphabet=alphabet)
        else:  # has counts
            res = NewMotif(alphabet)
            res.counts={}
            res.counts["A"]=self.counts["T"][::-1]
            res.counts["T"]=self.counts["A"][::-1]
            res.counts["G"]=self.counts["C"][::-1]
            res.counts["C"]=self.counts["G"][::-1]
            res.length=self.length
        res.__mask = self.__mask[::-1]
        return res

    @property
    def consensus(self):
        """Returns the consensus sequence.
        """
        return self.counts.consensus

    @property
    def anticonsensus(self):
        """returns the least probable pattern to be generated from this motif.
        """
        return self.counts.anticonsensus

    @property
    def degenerate_consensus(self):
        """Following the rules adapted from
D. R. Cavener: "Comparison of the consensus sequence flanking
translational start sites in Drosophila and vertebrates."
Nucleic Acids Research 15(4): 1353-1361. (1987).
The same rules are used by TRANSFAC."""
        return self.counts.degenerate_consensus

    def weblogo(self,fname,format="PNG",version="2.8.2", **kwds):
        """
        uses the Berkeley weblogo service to download and save a weblogo of
        itself

        requires an internet connection.
        The parameters from **kwds are passed directly to the weblogo server.

        Currently, this method uses WebLogo version 3.3.
        These are the arguments and their default values passed to
        WebLogo 3.3; see their website at http://weblogo.threeplusone.com
        for more information:

            'stack_width' : 'medium',
            'stack_per_line' : '40',
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
        import urllib
        import urllib2
        frequencies= self._to_transfac()
        url = 'http://weblogo.threeplusone.com/create.cgi'
        values = {'sequences' : frequencies,
                  'format' : format.lower(),
                  'stack_width' : 'medium',
                  'stack_per_line' : '40',
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
                  }
        for k,v in kwds.iteritems():
            if type(values[k])==bool:
                if not v:
                    v = ""
            values[k]=str(v)

        data = urllib.urlencode(values)
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req)
        f=open(fname,"w")
        im=response.read()

        f.write(im)
        f.close()

    def format(self,format):
        """Returns a string representation of the Motif in a given format

        Currently supported fromats:
         - pfm : JASPAR Position Frequency Matrix
         - transfac : TRANSFAC like files
        """

        if format=="pfm":
            from Bio.Motif import Jaspar
            return Jaspar.write(self)
        elif format=="transfac":
            from Bio.Motif import TRANSFAC
            motifs = [self]
            return TRANSFAC.write(motifs)
        else:
            raise ValueError("Unknown format type %s" % format)


def write(motifs, format):
    """Returns a string representation of motifs in a given format

    Currently supported fromats:
     - pfm : JASPAR Position Frequency Matrix
             [only if len(motifs)==1]
     - transfac : TRANSFAC like files
    """

    if format=="pfm":
        from Bio.Motif import Jaspar
        if len(motifs)!=1:
            raise Exception("Only a single motif can be written in the JASPAR Position Frequency Matrix (pfm) format")
        motif = motifs[0]
        return Jaspar.write(motif)
    elif format=="transfac":
        from Bio.Motif import TRANSFAC
        return TRANSFAC.write(motifs)
    else:
        raise ValueError("Unknown format type %s" % format)


NewMotif = Motif
from Bio.Motif._Motif import Motif
OldMotif = Motif

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
