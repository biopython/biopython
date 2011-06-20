# Copyright 2000-2002 Brad Chapman.
# Copyright 2004-2005 by M de Hoon.
# Copyright 2007-2010 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Provides objects to represent biological sequences with alphabets.

See also U{http://biopython.org/wiki/Seq} and the chapter in our tutorial:
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.html}
 - U{http://biopython.org/DIST/docs/tutorial/Tutorial.pdf}
"""
__docformat__ ="epytext en" #Don't just use plain text in epydoc API pages!

import string #for maketrans only
import array
import sys

from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data.IUPACData import ambiguous_dna_complement, ambiguous_rna_complement
from Bio.Data import CodonTable

def _maketrans(complement_mapping):
    """Makes a python string translation table (PRIVATE).

    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.

    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.
    
    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    before = ''.join(complement_mapping.keys())
    after  = ''.join(complement_mapping.values())
    before = before + before.lower()
    after  = after + after.lower()
    if sys.version_info[0] == 3 :
        return str.maketrans(before, after)
    else:
        return string.maketrans(before, after)

_dna_complement_table = _maketrans(ambiguous_dna_complement)
_rna_complement_table = _maketrans(ambiguous_rna_complement)

class Seq(object):
    """A read-only sequence object (essentially a string with an alphabet).

    Like normal python strings, our basic sequence object is immutable.
    This prevents you from doing my_seq[5] = "A" for example, but does allow
    Seq objects to be used as dictionary keys.

    The Seq object provides a number of string like methods (such as count,
    find, split and strip), which are alphabet aware where appropriate.

    In addition to the string like sequence, the Seq object has an alphabet
    property. This is an instance of an Alphabet class from Bio.Alphabet,
    for example generic DNA, or IUPAC DNA. This describes the type of molecule
    (e.g. RNA, DNA, protein) and may also indicate the expected symbols
    (letters).

    The Seq object also provides some biological methods, such as complement,
    reverse_complement, transcribe, back_transcribe and translate (which are
    not applicable to sequences with a protein alphabet).
    """
    def __init__(self, data, alphabet = Alphabet.generic_alphabet):
        """Create a Seq object.

        Arguments:
         - seq      - Sequence, required (string)
         - alphabet - Optional argument, an Alphabet object from Bio.Alphabet
        
        You will typically use Bio.SeqIO to read in sequences from files as
        SeqRecord objects, whose sequence will be exposed as a Seq object via
        the seq property.

        However, will often want to create your own Seq objects directly:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
        ...              IUPAC.protein)
        >>> my_seq
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())
        >>> print my_seq
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
        >>> my_seq.alphabet
        IUPACProtein()

        """
        # Enforce string storage
        if not isinstance(data, basestring):
            raise TypeError("The sequence data given to a Seq object should "
                            "be a string (not another Seq object etc)")
        self._data = data
        self.alphabet = alphabet  # Seq API requirement
 
    # A data property is/was a Seq API requirement
    # Note this is read only since the Seq object is meant to be imutable
    @property
    def data(self) :
        """Sequence as a string (DEPRECATED).

        This is a read only property provided for backwards compatility with
        older versions of Biopython (as is the tostring() method). We now
        encourage you to use str(my_seq) instead of my_seq.data or the method
        my_seq.tostring().

        In recent releases of Biopython it was possible to change a Seq object
        by updating its data property, but this triggered a deprecation warning.
        Now the data property is read only, since Seq objects are meant to be
        immutable:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_dna
        >>> my_seq = Seq("ACGT", generic_dna)
        >>> str(my_seq) == my_seq.tostring() == "ACGT"
        True
        >>> my_seq.data = "AAAA"
        Traceback (most recent call last):
           ...
        AttributeError: can't set attribute
        """
        import warnings
        import Bio
        warnings.warn("Accessing the .data attribute is deprecated. Please "
                      "use str(my_seq) or my_seq.tostring() instead of "
                      "my_seq.data.", Bio.BiopythonDeprecationWarning)
        return str(self)

    def __repr__(self):
        """Returns a (truncated) representation of the sequence for debugging."""
        if len(self) > 60:
            #Shows the last three letters as it is often useful to see if there
            #is a stop codon at the end of a sequence.
            #Note total length is 54+3+3=60
            return "%s('%s...%s', %s)" % (self.__class__.__name__,
                                   str(self)[:54], str(self)[-3:],
                                   repr(self.alphabet))
        else:
            return "%s(%s, %s)" % (self.__class__.__name__,
                                  repr(self._data),
                                   repr(self.alphabet))
    def __str__(self):
        """Returns the full sequence as a python string, use str(my_seq).

        Note that Biopython 1.44 and earlier would give a truncated
        version of repr(my_seq) for str(my_seq).  If you are writing code
        which need to be backwards compatible with old Biopython, you
        should continue to use my_seq.tostring() rather than str(my_seq).
        """
        return self._data

    def __hash__(self):
        """Hash for comparison.

        See the __cmp__ documentation - we plan to change this!
        """
        return id(self) #Currently use object identity for equality testing
    
    def __cmp__(self, other):
        """Compare the sequence to another sequence or a string (README).

        Historically comparing Seq objects has done Python object comparison.
        After considerable discussion (keeping in mind constraints of the
        Python language, hashes and dictionary support) a future release of
        Biopython will change this to use simple string comparison. The plan is
        that comparing incompatible alphabets (e.g. DNA to RNA) will trigger a
        warning.

        This version of Biopython still does Python object comparison, but with
        a warning about this future change. During this transition period,
        please just do explicit comparisons:

        >>> seq1 = Seq("ACGT")
        >>> seq2 = Seq("ACGT")
        >>> id(seq1) == id(seq2)
        False
        >>> str(seq1) == str(seq2)
        True

        Note - This method indirectly supports ==, < , etc.
        """
        if hasattr(other, "alphabet"):
            #other should be a Seq or a MutableSeq
            import warnings
            warnings.warn("In future comparing Seq objects will use string "
                          "comparison (not object comparison). Incompatible "
                          "alphabets will trigger a warning (not an exception). "
                          "In the interim please use id(seq1)==id(seq2) or "
                          "str(seq1)==str(seq2) to make your code explicit "
                          "and to avoid this warning.", FutureWarning)
        return cmp(id(self), id(other))

    def __len__(self):
        """Returns the length of the sequence, use len(my_seq)."""
        return len(self._data)       # Seq API requirement

    def __getitem__(self, index) :                 # Seq API requirement
        """Returns a subsequence of single letter, use my_seq[index]."""
        #Note since Python 2.0, __getslice__ is deprecated
        #and __getitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int):
            #Return a single letter as a string
            return self._data[index]
        else:
            #Return the (sub)sequence as another Seq object
            return Seq(self._data[index], self.alphabet)

    def __add__(self, other):
        """Add another sequence or string to this sequence.

        If adding a string to a Seq, the alphabet is preserved:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_protein
        >>> Seq("MELKI", generic_protein) + "LV"
        Seq('MELKILV', ProteinAlphabet())

        When adding two Seq (like) objects, the alphabets are important.
        Consider this example:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
        >>> unamb_dna_seq = Seq("ACGT", unambiguous_dna)
        >>> ambig_dna_seq = Seq("ACRGT", ambiguous_dna)
        >>> unamb_dna_seq
        Seq('ACGT', IUPACUnambiguousDNA())
        >>> ambig_dna_seq
        Seq('ACRGT', IUPACAmbiguousDNA())

        If we add the ambiguous and unambiguous IUPAC DNA alphabets, we get
        the more general ambiguous IUPAC DNA alphabet:
        
        >>> unamb_dna_seq + ambig_dna_seq
        Seq('ACGTACRGT', IUPACAmbiguousDNA())

        However, if the default generic alphabet is included, the result is
        a generic alphabet:

        >>> Seq("") + ambig_dna_seq
        Seq('ACRGT', Alphabet())

        You can't add RNA and DNA sequences:
        
        >>> from Bio.Alphabet import generic_dna, generic_rna
        >>> Seq("ACGT", generic_dna) + Seq("ACGU", generic_rna)
        Traceback (most recent call last):
           ...
        TypeError: Incompatible alphabets DNAAlphabet() and RNAAlphabet()

        You can't add nucleotide and protein sequences:

        >>> from Bio.Alphabet import generic_dna, generic_protein
        >>> Seq("ACGT", generic_dna) + Seq("MELKI", generic_protein)
        Traceback (most recent call last):
           ...
        TypeError: Incompatible alphabets DNAAlphabet() and ProteinAlphabet()
        """
        if hasattr(other, "alphabet"):
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]):
                raise TypeError("Incompatible alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            return self.__class__(str(self) + str(other), a)
        elif isinstance(other, basestring):
            #other is a plain string - use the current alphabet
            return self.__class__(str(self) + other, self.alphabet)
        from Bio.SeqRecord import SeqRecord #Lazy to avoid circular imports
        if isinstance(other, SeqRecord):
            #Get the SeqRecord's __radd__ to handle this
            return NotImplemented
        else :
            raise TypeError

    def __radd__(self, other):
        """Adding a sequence on the left.

        If adding a string to a Seq, the alphabet is preserved:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_protein
        >>> "LV" + Seq("MELKI", generic_protein)
        Seq('LVMELKI', ProteinAlphabet())

        Adding two Seq (like) objects is handled via the __add__ method.
        """
        if hasattr(other, "alphabet"):
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]):
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            return self.__class__(str(other) + str(self), a)
        elif isinstance(other, basestring):
            #other is a plain string - use the current alphabet
            return self.__class__(other + str(self), self.alphabet)
        else:
            raise TypeError

    def tostring(self):                            # Seq API requirement
        """Returns the full sequence as a python string (semi-obsolete).

        Although not formally deprecated, you are now encouraged to use
        str(my_seq) instead of my_seq.tostring()."""
        #TODO - Fix all places elsewhere in Biopython using this method,
        #then start deprecation process?
        #import warnings
        #warnings.warn("This method is obsolete; please use str(my_seq) "
        #              "instead of my_seq.tostring().",
        #              PendingDeprecationWarning)
        return str(self)
    
    def tomutable(self):   # Needed?  Or use a function?
        """Returns the full sequence as a MutableSeq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAAL",
        ...              IUPAC.protein)
        >>> my_seq
        Seq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())
        >>> my_seq.tomutable()
        MutableSeq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())

        Note that the alphabet is preserved.
        """
        return MutableSeq(str(self), self.alphabet)

    def _get_seq_str_and_check_alphabet(self, other_sequence):
        """string/Seq/MutableSeq to string, checking alphabet (PRIVATE).

        For a string argument, returns the string.

        For a Seq or MutableSeq, it checks the alphabet is compatible
        (raising an exception if it isn't), and then returns a string.
        """
        try:
            other_alpha = other_sequence.alphabet
        except AttributeError:
            #Assume other_sequence is a string
            return other_sequence

        #Other should be a Seq or a MutableSeq
        if not Alphabet._check_type_compatible([self.alphabet, other_alpha]):
            raise TypeError("Incompatable alphabets %s and %s" \
                            % (repr(self.alphabet), repr(other_alpha)))
        #Return as a string
        return str(other_sequence)
    
    def count(self, sub, start=0, end=sys.maxint):
        """Non-overlapping count method, like that of a python string.

        This behaves like the python string method of the same name,
        which does a non-overlapping count!

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.
    
        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        e.g.

        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("AAAATGA")
        >>> print my_seq.count("A")
        5
        >>> print my_seq.count("ATG")
        1
        >>> print my_seq.count(Seq("AT"))
        1
        >>> print my_seq.count("AT", 2, -1)
        1

        HOWEVER, please note because python strings and Seq objects (and
        MutableSeq objects) do a non-overlapping search, this may not give
        the answer you expect:

        >>> "AAAA".count("AA")
        2
        >>> print Seq("AAAA").count("AA")
        2

        A non-overlapping search would give the answer as three!
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return str(self).count(sub_str, start, end)

    def __contains__(self, char):
        """Implements the 'in' keyword, like a python string.

        e.g.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_dna, generic_rna, generic_protein
        >>> my_dna = Seq("ATATGAAATTTGAAAA", generic_dna)
        >>> "AAA" in my_dna
        True
        >>> Seq("AAA") in my_dna
        True
        >>> Seq("AAA", generic_dna) in my_dna
        True

        Like other Seq methods, this will raise a type error if another Seq
        (or Seq like) object with an incompatible alphabet is used:

        >>> Seq("AAA", generic_rna) in my_dna
        Traceback (most recent call last):
           ...
        TypeError: Incompatable alphabets DNAAlphabet() and RNAAlphabet()
        >>> Seq("AAA", generic_protein) in my_dna
        Traceback (most recent call last):
           ...
        TypeError: Incompatable alphabets DNAAlphabet() and ProteinAlphabet()
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(char)
        return sub_str in str(self)

    def find(self, sub, start=0, end=sys.maxint):
        """Find method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.
        
        e.g. Locating the first typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.find("AUG")
        3
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return str(self).find(sub_str, start, end)

    def rfind(self, sub, start=0, end=sys.maxint):
        """Find from right method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the last (right most) occurrence of
        substring argument sub in the (sub)sequence given by [start:end].

        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.

        e.g. Locating the last typical start codon, AUG, in an RNA sequence:

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.rfind("AUG")
        15
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return str(self).rfind(sub_str, start, end)

    def startswith(self, prefix, start=0, end=sys.maxint):
        """Does the Seq start with the given prefix?  Returns True/False.

        This behaves like the python string method of the same name.

        Return True if the sequence starts with the specified prefix
        (a string or another Seq object), False otherwise.
        With optional start, test sequence beginning at that position.
        With optional end, stop comparing sequence at that position.
        prefix can also be a tuple of strings to try.  e.g.
        
        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.startswith("GUC")
        True
        >>> my_rna.startswith("AUG")
        False
        >>> my_rna.startswith("AUG", 3)
        True
        >>> my_rna.startswith(("UCC","UCA","UCG"),1)
        True
        """
        #If it has one, check the alphabet:
        if isinstance(prefix, tuple):
            #TODO - Once we drop support for Python 2.4, instead of this
            #loop offload to the string method (requires Python 2.5+).
            #Check all the alphabets first...
            prefix_strings = [self._get_seq_str_and_check_alphabet(p) \
                              for p in prefix]
            for prefix_str in prefix_strings:
                if str(self).startswith(prefix_str, start, end):
                    return True
            return False
        else:
            prefix_str = self._get_seq_str_and_check_alphabet(prefix)
            return str(self).startswith(prefix_str, start, end)

    def endswith(self, suffix, start=0, end=sys.maxint):
        """Does the Seq end with the given suffix?  Returns True/False.

        This behaves like the python string method of the same name.

        Return True if the sequence ends with the specified suffix
        (a string or another Seq object), False otherwise.
        With optional start, test sequence beginning at that position.
        With optional end, stop comparing sequence at that position.
        suffix can also be a tuple of strings to try.  e.g.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_rna.endswith("UUG")
        True
        >>> my_rna.endswith("AUG")
        False
        >>> my_rna.endswith("AUG", 0, 18)
        True
        >>> my_rna.endswith(("UCC","UCA","UUG"))
        True
        """        
        #If it has one, check the alphabet:
        if isinstance(suffix, tuple):
            #TODO - Once we drop support for Python 2.4, instead of this
            #loop offload to the string method (requires Python 2.5+).
            #Check all the alphabets first...
            suffix_strings = [self._get_seq_str_and_check_alphabet(p) \
                              for p in suffix]
            for suffix_str in suffix_strings:
                if str(self).endswith(suffix_str, start, end):
                    return True
            return False
        else:
            suffix_str = self._get_seq_str_and_check_alphabet(suffix)
            return str(self).endswith(suffix_str, start, end)


    def split(self, sep=None, maxsplit=-1):
        """Split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done.  If maxsplit is ommited, all
        splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.
        
        e.g.

        >>> from Bio.Seq import Seq
        >>> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
        >>> my_aa = my_rna.translate()
        >>> my_aa
        Seq('VMAIVMGR*KGAR*L', HasStopCodon(ExtendedIUPACProtein(), '*'))
        >>> my_aa.split("*")
        [Seq('VMAIVMGR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('KGAR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('L', HasStopCodon(ExtendedIUPACProtein(), '*'))]
        >>> my_aa.split("*",1)
        [Seq('VMAIVMGR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('KGAR*L', HasStopCodon(ExtendedIUPACProtein(), '*'))]

        See also the rsplit method:

        >>> my_aa.rsplit("*",1)
        [Seq('VMAIVMGR*KGAR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('L', HasStopCodon(ExtendedIUPACProtein(), '*'))]
        """
        #If it has one, check the alphabet:
        sep_str = self._get_seq_str_and_check_alphabet(sep)
        #TODO - If the sep is the defined stop symbol, or gap char,
        #should we adjust the alphabet?
        return [Seq(part, self.alphabet) \
                for part in str(self).split(sep_str, maxsplit)]

    def rsplit(self, sep=None, maxsplit=-1):
        """Right split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done COUNTING FROM THE RIGHT.
        If maxsplit is ommited, all splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.
        
        e.g. print my_seq.rsplit("*",1)

        See also the split method.
        """
        #If it has one, check the alphabet:
        sep_str = self._get_seq_str_and_check_alphabet(sep)
        return [Seq(part, self.alphabet) \
                for part in str(self).rsplit(sep_str, maxsplit)]

    def strip(self, chars=None):
        """Returns a new Seq object with leading and trailing ends stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        ommitted or None (default) then as for the python string method,
        this defaults to removing any white space.
        
        e.g. print my_seq.strip("-")

        See also the lstrip and rstrip methods.
        """
        #If it has one, check the alphabet:
        strip_str = self._get_seq_str_and_check_alphabet(chars)
        return Seq(str(self).strip(strip_str), self.alphabet)

    def lstrip(self, chars=None):
        """Returns a new Seq object with leading (left) end stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        ommitted or None (default) then as for the python string method,
        this defaults to removing any white space.
        
        e.g. print my_seq.lstrip("-")

        See also the strip and rstrip methods.
        """
        #If it has one, check the alphabet:
        strip_str = self._get_seq_str_and_check_alphabet(chars)
        return Seq(str(self).lstrip(strip_str), self.alphabet)

    def rstrip(self, chars=None):
        """Returns a new Seq object with trailing (right) end stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        ommitted or None (default) then as for the python string method,
        this defaults to removing any white space.
        
        e.g. Removing a nucleotide sequence's polyadenylation (poly-A tail):

        >>> from Bio.Alphabet import IUPAC
        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("CGGTACGCTTATGTCACGTAGAAAAAA", IUPAC.unambiguous_dna)
        >>> my_seq
        Seq('CGGTACGCTTATGTCACGTAGAAAAAA', IUPACUnambiguousDNA())
        >>> my_seq.rstrip("A")
        Seq('CGGTACGCTTATGTCACGTAG', IUPACUnambiguousDNA())

        See also the strip and lstrip methods.
        """
        #If it has one, check the alphabet:
        strip_str = self._get_seq_str_and_check_alphabet(chars)
        return Seq(str(self).rstrip(strip_str), self.alphabet)

    def upper(self):
        """Returns an upper case copy of the sequence.

        >>> from Bio.Alphabet import HasStopCodon, generic_protein
        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("VHLTPeeK*", HasStopCodon(generic_protein))
        >>> my_seq
        Seq('VHLTPeeK*', HasStopCodon(ProteinAlphabet(), '*'))
        >>> my_seq.lower()
        Seq('vhltpeek*', HasStopCodon(ProteinAlphabet(), '*'))
        >>> my_seq.upper()
        Seq('VHLTPEEK*', HasStopCodon(ProteinAlphabet(), '*'))

        This will adjust the alphabet if required. See also the lower method.
        """
        return Seq(str(self).upper(), self.alphabet._upper())

    def lower(self):
        """Returns a lower case copy of the sequence.

        This will adjust the alphabet if required. Note that the IUPAC alphabets
        are upper case only, and thus a generic alphabet must be substituted.

        >>> from Bio.Alphabet import Gapped, generic_dna
        >>> from Bio.Alphabet import IUPAC
        >>> from Bio.Seq import Seq
        >>> my_seq = Seq("CGGTACGCTTATGTCACGTAG*AAAAAA", Gapped(IUPAC.unambiguous_dna, "*"))
        >>> my_seq
        Seq('CGGTACGCTTATGTCACGTAG*AAAAAA', Gapped(IUPACUnambiguousDNA(), '*'))
        >>> my_seq.lower()
        Seq('cggtacgcttatgtcacgtag*aaaaaa', Gapped(DNAAlphabet(), '*'))

        See also the upper method.
        """
        return Seq(str(self).lower(), self.alphabet._lower())

    def complement(self):
        """Returns the complement sequence. New Seq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_dna = Seq("CCCCCGATAG", IUPAC.unambiguous_dna)
        >>> my_dna
        Seq('CCCCCGATAG', IUPACUnambiguousDNA())
        >>> my_dna.complement()
        Seq('GGGGGCTATC', IUPACUnambiguousDNA())

        You can of course used mixed case sequences,

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_dna
        >>> my_dna = Seq("CCCCCgatA-GD", generic_dna)
        >>> my_dna
        Seq('CCCCCgatA-GD', DNAAlphabet())
        >>> my_dna.complement()
        Seq('GGGGGctaT-CH', DNAAlphabet())

        Note in the above example, ambiguous character D denotes
        G, A or T so its complement is H (for C, T or A).
        
        Trying to complement a protein sequence raises an exception.

        >>> my_protein = Seq("MAIVMGR", IUPAC.protein)
        >>> my_protein.complement()
        Traceback (most recent call last):
           ...
        ValueError: Proteins do not have complements!
        """
        base = Alphabet._get_base_alphabet(self.alphabet)
        if isinstance(base, Alphabet.ProteinAlphabet):
            raise ValueError("Proteins do not have complements!")
        if isinstance(base, Alphabet.DNAAlphabet):
            ttable = _dna_complement_table
        elif isinstance(base, Alphabet.RNAAlphabet):
            ttable = _rna_complement_table
        elif ('U' in self._data or 'u' in self._data) \
        and ('T' in self._data or 't' in self._data):
            #TODO - Handle this cleanly?
            raise ValueError("Mixed RNA/DNA found")
        elif 'U' in self._data or 'u' in self._data:
            ttable = _rna_complement_table
        else:
            ttable = _dna_complement_table
        #Much faster on really long sequences than the previous loop based one.
        #thx to Michael Palmer, University of Waterloo
        return Seq(str(self).translate(ttable), self.alphabet)

    def reverse_complement(self):
        """Returns the reverse complement sequence. New Seq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_dna = Seq("CCCCCGATAGNR", IUPAC.ambiguous_dna)
        >>> my_dna
        Seq('CCCCCGATAGNR', IUPACAmbiguousDNA())
        >>> my_dna.reverse_complement()
        Seq('YNCTATCGGGGG', IUPACAmbiguousDNA())

        Note in the above example, since R = G or A, its complement
        is Y (which denotes C or T).

        You can of course used mixed case sequences,

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_dna
        >>> my_dna = Seq("CCCCCgatA-G", generic_dna)
        >>> my_dna
        Seq('CCCCCgatA-G', DNAAlphabet())
        >>> my_dna.reverse_complement()
        Seq('C-TatcGGGGG', DNAAlphabet())

        Trying to complement a protein sequence raises an exception:

        >>> my_protein = Seq("MAIVMGR", IUPAC.protein)
        >>> my_protein.reverse_complement()
        Traceback (most recent call last):
           ...
        ValueError: Proteins do not have complements!
        """
        #Use -1 stride/step to reverse the complement
        return self.complement()[::-1]

    def transcribe(self):
        """Returns the RNA sequence from a DNA sequence. New Seq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
        ...                  IUPAC.unambiguous_dna)
        >>> coding_dna
        Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())
        >>> coding_dna.transcribe()
        Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

        Trying to transcribe a protein or RNA sequence raises an exception:

        >>> my_protein = Seq("MAIVMGR", IUPAC.protein)
        >>> my_protein.transcribe()
        Traceback (most recent call last):
           ...
        ValueError: Proteins cannot be transcribed!
        """
        base = Alphabet._get_base_alphabet(self.alphabet)
        if isinstance(base, Alphabet.ProteinAlphabet):
            raise ValueError("Proteins cannot be transcribed!")
        if isinstance(base, Alphabet.RNAAlphabet):
            raise ValueError("RNA cannot be transcribed!")

        if self.alphabet==IUPAC.unambiguous_dna:
            alphabet = IUPAC.unambiguous_rna
        elif self.alphabet==IUPAC.ambiguous_dna:
            alphabet = IUPAC.ambiguous_rna
        else:
            alphabet = Alphabet.generic_rna
        return Seq(str(self).replace('T','U').replace('t','u'), alphabet)
    
    def back_transcribe(self):
        """Returns the DNA sequence from an RNA sequence. New Seq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG",
        ...                     IUPAC.unambiguous_rna)
        >>> messenger_rna
        Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())
        >>> messenger_rna.back_transcribe()
        Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())

        Trying to back-transcribe a protein or DNA sequence raises an
        exception:

        >>> my_protein = Seq("MAIVMGR", IUPAC.protein)
        >>> my_protein.back_transcribe()
        Traceback (most recent call last):
           ...
        ValueError: Proteins cannot be back transcribed!
        """
        base = Alphabet._get_base_alphabet(self.alphabet)
        if isinstance(base, Alphabet.ProteinAlphabet):
            raise ValueError("Proteins cannot be back transcribed!")
        if isinstance(base, Alphabet.DNAAlphabet):
            raise ValueError("DNA cannot be back transcribed!")

        if self.alphabet==IUPAC.unambiguous_rna:
            alphabet = IUPAC.unambiguous_dna
        elif self.alphabet==IUPAC.ambiguous_rna:
            alphabet = IUPAC.ambiguous_dna
        else:
            alphabet = Alphabet.generic_dna
        return Seq(str(self).replace("U", "T").replace("u", "t"), alphabet)

    def translate(self, table="Standard", stop_symbol="*", to_stop=False,
                  cds=False):
        """Turns a nucleotide sequence into a protein sequence. New Seq object.

        This method will translate DNA or RNA sequences, and those with a
        nucleotide or generic alphabet.  Trying to translate a protein
        sequence raises an exception.

        Arguments:
         - table - Which codon table to use?  This can be either a name
                   (string), an NCBI identifier (integer), or a CodonTable
                   object (useful for non-standard genetic codes).  This
                   defaults to the "Standard" table.
         - stop_symbol - Single character string, what to use for terminators.
                         This defaults to the asterisk, "*".
         - to_stop - Boolean, defaults to False meaning do a full translation
                     continuing on past any stop codons (translated as the
                     specified stop_symbol).  If True, translation is
                     terminated at the first in frame stop codon (and the
                     stop_symbol is not appended to the returned protein
                     sequence).
         - cds - Boolean, indicates this is a complete CDS.  If True,
                 this checks the sequence starts with a valid alternative start
                 codon (which will be translated as methionine, M), that the
                 sequence length is a multiple of three, and that there is a
                 single in frame stop codon at the end (this will be excluded
                 from the protein sequence, regardless of the to_stop option).
                 If these tests fail, an exception is raised.
        
        e.g. Using the standard table:

        >>> coding_dna = Seq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
        >>> coding_dna.translate()
        Seq('VAIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
        >>> coding_dna.translate(stop_symbol="@")
        Seq('VAIVMGR@KGAR@', HasStopCodon(ExtendedIUPACProtein(), '@'))
        >>> coding_dna.translate(to_stop=True)
        Seq('VAIVMGR', ExtendedIUPACProtein())

        Now using NCBI table 2, where TGA is not a stop codon:

        >>> coding_dna.translate(table=2)
        Seq('VAIVMGRWKGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
        >>> coding_dna.translate(table=2, to_stop=True)
        Seq('VAIVMGRWKGAR', ExtendedIUPACProtein())

        In fact, GTG is an alternative start codon under NCBI table 2, meaning
        this sequence could be a complete CDS:

        >>> coding_dna.translate(table=2, cds=True)
        Seq('MAIVMGRWKGAR', ExtendedIUPACProtein())

        It isn't a valid CDS under NCBI table 1, due to both the start codon and
        also the in frame stop codons:
        
        >>> coding_dna.translate(table=1, cds=True)
        Traceback (most recent call last):
            ...
        TranslationError: First codon 'GTG' is not a start codon

        If the sequence has no in-frame stop codon, then the to_stop argument
        has no effect:

        >>> coding_dna2 = Seq("TTGGCCATTGTAATGGGCCGC")
        >>> coding_dna2.translate()
        Seq('LAIVMGR', ExtendedIUPACProtein())
        >>> coding_dna2.translate(to_stop=True)
        Seq('LAIVMGR', ExtendedIUPACProtein())

        NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
        or a stop codon.  These are translated as "X".  Any invalid codon
        (e.g. "TA?" or "T-A") will throw a TranslationError.

        NOTE - Does NOT support gapped sequences.

        NOTE - This does NOT behave like the python string's translate
        method.  For that use str(my_seq).translate(...) instead.
        """
        if isinstance(table, str) and len(table)==256:
            raise ValueError("The Seq object translate method DOES NOT take " \
                             + "a 256 character string mapping table like " \
                             + "the python string object's translate method. " \
                             + "Use str(my_seq).translate(...) instead.")
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet):
            raise ValueError("Proteins cannot be translated!")
        try:
            table_id = int(table)
        except ValueError:
            #Assume its a table name
            if self.alphabet==IUPAC.unambiguous_dna:
                #Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_dna_by_name[table]
            elif self.alphabet==IUPAC.unambiguous_rna:
                #Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_rna_by_name[table]
            else:
                #This will use the extended IUPAC protein alphabet with X etc.
                #The same table can be used for RNA or DNA (we use this for
                #translating strings).
                codon_table = CodonTable.ambiguous_generic_by_name[table]
        except (AttributeError, TypeError):
            #Assume its a CodonTable object
            if isinstance(table, CodonTable.CodonTable):
                codon_table = table
            else:
                raise ValueError('Bad table argument')
        else:
            #Assume its a table ID
            if self.alphabet==IUPAC.unambiguous_dna:
                #Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_dna_by_id[table_id]
            elif self.alphabet==IUPAC.unambiguous_rna:
                #Will use standard IUPAC protein alphabet, no need for X
                codon_table = CodonTable.unambiguous_rna_by_id[table_id]
            else:
                #This will use the extended IUPAC protein alphabet with X etc.
                #The same table can be used for RNA or DNA (we use this for
                #translating strings).
                codon_table = CodonTable.ambiguous_generic_by_id[table_id]
        protein = _translate_str(str(self), codon_table, \
                                 stop_symbol, to_stop, cds)
        if stop_symbol in protein:
            alphabet = Alphabet.HasStopCodon(codon_table.protein_alphabet,
                                             stop_symbol = stop_symbol)
        else:
            alphabet = codon_table.protein_alphabet
        return Seq(protein, alphabet)

    def ungap(self, gap=None):
        """Return a copy of the sequence without the gap character(s).

        The gap character can be specified in two ways - either as an explicit
        argument, or via the sequence's alphabet. For example:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_dna
        >>> my_dna = Seq("-ATA--TGAAAT-TTGAAAA", generic_dna)
        >>> my_dna
        Seq('-ATA--TGAAAT-TTGAAAA', DNAAlphabet())
        >>> my_dna.ungap("-")
        Seq('ATATGAAATTTGAAAA', DNAAlphabet())

        If the gap character is not given as an argument, it will be taken from
        the sequence's alphabet (if defined). Notice that the returned sequence's
        alphabet is adjusted since it no longer requires a gapped alphabet:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC, Gapped, HasStopCodon
        >>> my_pro = Seq("MVVLE=AD*", HasStopCodon(Gapped(IUPAC.protein, "=")))
        >>> my_pro
        Seq('MVVLE=AD*', HasStopCodon(Gapped(IUPACProtein(), '='), '*'))
        >>> my_pro.ungap()
        Seq('MVVLEAD*', HasStopCodon(IUPACProtein(), '*'))

        Or, with a simpler gapped DNA example:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC, Gapped
        >>> my_seq = Seq("CGGGTAG=AAAAAA", Gapped(IUPAC.unambiguous_dna, "="))
        >>> my_seq
        Seq('CGGGTAG=AAAAAA', Gapped(IUPACUnambiguousDNA(), '='))
        >>> my_seq.ungap()
        Seq('CGGGTAGAAAAAA', IUPACUnambiguousDNA())

        As long as it is consistent with the alphabet, although it is redundant,
        you can still supply the gap character as an argument to this method:

        >>> my_seq
        Seq('CGGGTAG=AAAAAA', Gapped(IUPACUnambiguousDNA(), '='))
        >>> my_seq.ungap("=")
        Seq('CGGGTAGAAAAAA', IUPACUnambiguousDNA())
        
        However, if the gap character given as the argument disagrees with that
        declared in the alphabet, an exception is raised:

        >>> my_seq
        Seq('CGGGTAG=AAAAAA', Gapped(IUPACUnambiguousDNA(), '='))
        >>> my_seq.ungap("-")
        Traceback (most recent call last):
           ...
        ValueError: Gap '-' does not match '=' from alphabet

        Finally, if a gap character is not supplied, and the alphabet does not
        define one, an exception is raised:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import generic_dna
        >>> my_dna = Seq("ATA--TGAAAT-TTGAAAA", generic_dna)
        >>> my_dna
        Seq('ATA--TGAAAT-TTGAAAA', DNAAlphabet())
        >>> my_dna.ungap()
        Traceback (most recent call last):
           ...
        ValueError: Gap character not given and not defined in alphabet

        """
        if hasattr(self.alphabet, "gap_char"):
            if not gap:
                gap = self.alphabet.gap_char
            elif gap != self.alphabet.gap_char:
                raise ValueError("Gap %s does not match %s from alphabet" \
                                 % (repr(gap), repr(self.alphabet.gap_char)))
            alpha = Alphabet._ungap(self.alphabet)
        elif not gap:
            raise ValueError("Gap character not given and not defined in alphabet")
        else:
            alpha = self.alphabet #modify!
        if len(gap)!=1 or not isinstance(gap, str):
            raise ValueError("Unexpected gap character, %s" % repr(gap))
        return Seq(str(self).replace(gap, ""), alpha)

class UnknownSeq(Seq):
    """A read-only sequence object of known length but unknown contents.

    If you have an unknown sequence, you can represent this with a normal
    Seq object, for example:

    >>> my_seq = Seq("N"*5)
    >>> my_seq
    Seq('NNNNN', Alphabet())
    >>> len(my_seq)
    5
    >>> print my_seq
    NNNNN

    However, this is rather wasteful of memory (especially for large
    sequences), which is where this class is most usefull:

    >>> unk_five = UnknownSeq(5)
    >>> unk_five
    UnknownSeq(5, alphabet = Alphabet(), character = '?')
    >>> len(unk_five)
    5
    >>> print(unk_five)
    ?????

    You can add unknown sequence together, provided their alphabets and
    characters are compatible, and get another memory saving UnknownSeq:

    >>> unk_four = UnknownSeq(4)
    >>> unk_four
    UnknownSeq(4, alphabet = Alphabet(), character = '?')
    >>> unk_four + unk_five
    UnknownSeq(9, alphabet = Alphabet(), character = '?')

    If the alphabet or characters don't match up, the addition gives an
    ordinary Seq object:
    
    >>> unk_nnnn = UnknownSeq(4, character = "N")
    >>> unk_nnnn
    UnknownSeq(4, alphabet = Alphabet(), character = 'N')
    >>> unk_nnnn + unk_four
    Seq('NNNN????', Alphabet())

    Combining with a real Seq gives a new Seq object:

    >>> known_seq = Seq("ACGT")
    >>> unk_four + known_seq
    Seq('????ACGT', Alphabet())
    >>> known_seq + unk_four
    Seq('ACGT????', Alphabet())
    """
    def __init__(self, length, alphabet = Alphabet.generic_alphabet, character = None):
        """Create a new UnknownSeq object.

        If character is ommited, it is determed from the alphabet, "N" for
        nucleotides, "X" for proteins, and "?" otherwise.
        """
        self._length = int(length)
        if self._length < 0:
            #TODO - Block zero length UnknownSeq?  You can just use a Seq!
            raise ValueError("Length must not be negative.")
        self.alphabet = alphabet
        if character:
            if len(character) != 1:
                raise ValueError("character argument should be a single letter string.")
            self._character = character
        else:
            base = Alphabet._get_base_alphabet(alphabet)
            #TODO? Check the case of the letters in the alphabet?
            #We may have to use "n" instead of "N" etc.
            if isinstance(base, Alphabet.NucleotideAlphabet):
                self._character = "N"
            elif isinstance(base, Alphabet.ProteinAlphabet):
                self._character = "X"
            else:
                self._character = "?"

    def __len__(self):
        """Returns the stated length of the unknown sequence."""
        return self._length
    
    def __str__(self):
        """Returns the unknown sequence as full string of the given length."""
        return self._character * self._length

    def __repr__(self):
        return "UnknownSeq(%i, alphabet = %s, character = %s)" \
               % (self._length, repr(self.alphabet), repr(self._character))

    def __add__(self, other):
        """Add another sequence or string to this sequence.

        Adding two UnknownSeq objects returns another UnknownSeq object
        provided the character is the same and the alphabets are compatible.

        >>> from Bio.Seq import UnknownSeq
        >>> from Bio.Alphabet import generic_protein
        >>> UnknownSeq(10, generic_protein) + UnknownSeq(5, generic_protein)
        UnknownSeq(15, alphabet = ProteinAlphabet(), character = 'X')

        If the characters differ, an UnknownSeq object cannot be used, so a
        Seq object is returned:

        >>> from Bio.Seq import UnknownSeq
        >>> from Bio.Alphabet import generic_protein
        >>> UnknownSeq(10, generic_protein) + UnknownSeq(5, generic_protein,
        ...                                              character="x")
        Seq('XXXXXXXXXXxxxxx', ProteinAlphabet())

        If adding a string to an UnknownSeq, a new Seq is returned with the
        same alphabet:
        
        >>> from Bio.Seq import UnknownSeq
        >>> from Bio.Alphabet import generic_protein
        >>> UnknownSeq(5, generic_protein) + "LV"
        Seq('XXXXXLV', ProteinAlphabet())
        """
        if isinstance(other, UnknownSeq) \
        and other._character == self._character:
            #TODO - Check the alphabets match
            return UnknownSeq(len(self)+len(other),
                              self.alphabet, self._character)
        #Offload to the base class...
        return Seq(str(self), self.alphabet) + other

    def __radd__(self, other):
        #If other is an UnknownSeq, then __add__ would be called.
        #Offload to the base class...
        return other + Seq(str(self), self.alphabet)

    def __getitem__(self, index):
        """Get a subsequence from the UnknownSeq object.
        
        >>> unk = UnknownSeq(8, character="N")
        >>> print unk[:]
        NNNNNNNN
        >>> print unk[5:3]
        <BLANKLINE>
        >>> print unk[1:-1]
        NNNNNN
        >>> print unk[1:-1:2]
        NNN
        """
        if isinstance(index, int):
            #TODO - Check the bounds without wasting memory
            return str(self)[index]
        old_length = self._length
        step = index.step
        if step is None or step == 1:
            #This calculates the length you'd get from ("N"*old_length)[index]
            start = index.start
            end = index.stop
            if start is None:
                start = 0
            elif start < 0:
                start = max(0, old_length + start)
            elif start > old_length:
                start = old_length
            if end is None:
                end = old_length
            elif end < 0:
                end = max(0, old_length + end)
            elif end > old_length:
                end = old_length
            new_length = max(0, end-start)
        elif step == 0:
            raise ValueError("slice step cannot be zero")
        else:
            #TODO - handle step efficiently
            new_length = len(("X"*old_length)[index])
        #assert new_length == len(("X"*old_length)[index]), \
        #       (index, start, end, step, old_length,
        #        new_length, len(("X"*old_length)[index]))
        return UnknownSeq(new_length, self.alphabet, self._character)

    def count(self, sub, start=0, end=sys.maxint):
        """Non-overlapping count method, like that of a python string.

        This behaves like the python string (and Seq object) method of the
        same name, which does a non-overlapping count!

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.
    
        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        >>> "NNNN".count("N")
        4
        >>> Seq("NNNN").count("N")
        4
        >>> UnknownSeq(4, character="N").count("N")
        4
        >>> UnknownSeq(4, character="N").count("A")
        0
        >>> UnknownSeq(4, character="N").count("AA")
        0

        HOWEVER, please note because that python strings and Seq objects (and
        MutableSeq objects) do a non-overlapping search, this may not give
        the answer you expect:

        >>> UnknownSeq(4, character="N").count("NN")
        2
        >>> UnknownSeq(4, character="N").count("NNN")
        1
        """
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        if len(sub_str) == 1:
            if str(sub_str) == self._character:
                if start==0 and end >= self._length:
                    return self._length
                else:
                    #This could be done more cleverly...
                    return str(self).count(sub_str, start, end)
            else:
                return 0
        else:
            if set(sub_str) == set(self._character):
                if start==0 and end >= self._length:
                    return self._length // len(sub_str)
                else:
                    #This could be done more cleverly...
                    return str(self).count(sub_str, start, end)
            else:
                return 0

    def complement(self):
        """The complement of an unknown nucleotide equals itself.

        >>> my_nuc = UnknownSeq(8)
        >>> my_nuc
        UnknownSeq(8, alphabet = Alphabet(), character = '?')
        >>> print my_nuc
        ????????
        >>> my_nuc.complement()
        UnknownSeq(8, alphabet = Alphabet(), character = '?')
        >>> print my_nuc.complement()
        ????????
        """
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet):
            raise ValueError("Proteins do not have complements!")
        return self

    def reverse_complement(self):
        """The reverse complement of an unknown nucleotide equals itself.

        >>> my_nuc = UnknownSeq(10)
        >>> my_nuc
        UnknownSeq(10, alphabet = Alphabet(), character = '?')
        >>> print my_nuc
        ??????????
        >>> my_nuc.reverse_complement()
        UnknownSeq(10, alphabet = Alphabet(), character = '?')
        >>> print my_nuc.reverse_complement()
        ??????????
        """
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet):
            raise ValueError("Proteins do not have complements!")
        return self

    def transcribe(self):
        """Returns unknown RNA sequence from an unknown DNA sequence.

        >>> my_dna = UnknownSeq(10, character="N")
        >>> my_dna
        UnknownSeq(10, alphabet = Alphabet(), character = 'N')
        >>> print my_dna
        NNNNNNNNNN
        >>> my_rna = my_dna.transcribe()
        >>> my_rna
        UnknownSeq(10, alphabet = RNAAlphabet(), character = 'N')
        >>> print my_rna
        NNNNNNNNNN
        """
        #Offload the alphabet stuff
        s = Seq(self._character, self.alphabet).transcribe()
        return UnknownSeq(self._length, s.alphabet, self._character)

    def back_transcribe(self):
        """Returns unknown DNA sequence from an unknown RNA sequence.

        >>> my_rna = UnknownSeq(20, character="N")
        >>> my_rna
        UnknownSeq(20, alphabet = Alphabet(), character = 'N')
        >>> print my_rna
        NNNNNNNNNNNNNNNNNNNN
        >>> my_dna = my_rna.back_transcribe()
        >>> my_dna
        UnknownSeq(20, alphabet = DNAAlphabet(), character = 'N')
        >>> print my_dna
        NNNNNNNNNNNNNNNNNNNN
        """
        #Offload the alphabet stuff
        s = Seq(self._character, self.alphabet).back_transcribe()
        return UnknownSeq(self._length, s.alphabet, self._character)

    def upper(self):
        """Returns an upper case copy of the sequence.

        >>> from Bio.Alphabet import generic_dna
        >>> from Bio.Seq import UnknownSeq
        >>> my_seq = UnknownSeq(20, generic_dna, character="n")
        >>> my_seq
        UnknownSeq(20, alphabet = DNAAlphabet(), character = 'n')
        >>> print my_seq
        nnnnnnnnnnnnnnnnnnnn
        >>> my_seq.upper()
        UnknownSeq(20, alphabet = DNAAlphabet(), character = 'N')
        >>> print my_seq.upper()
        NNNNNNNNNNNNNNNNNNNN

        This will adjust the alphabet if required. See also the lower method.
        """
        return UnknownSeq(self._length, self.alphabet._upper(), self._character.upper())

    def lower(self):
        """Returns a lower case copy of the sequence.

        This will adjust the alphabet if required:

        >>> from Bio.Alphabet import IUPAC
        >>> from Bio.Seq import UnknownSeq
        >>> my_seq = UnknownSeq(20, IUPAC.extended_protein)
        >>> my_seq
        UnknownSeq(20, alphabet = ExtendedIUPACProtein(), character = 'X')
        >>> print my_seq
        XXXXXXXXXXXXXXXXXXXX
        >>> my_seq.lower()
        UnknownSeq(20, alphabet = ProteinAlphabet(), character = 'x')
        >>> print my_seq.lower()
        xxxxxxxxxxxxxxxxxxxx

        See also the upper method.
        """
        return UnknownSeq(self._length, self.alphabet._lower(), self._character.lower())

    def translate(self, **kwargs):
        """Translate an unknown nucleotide sequence into an unknown protein.

        e.g.

        >>> my_seq = UnknownSeq(11, character="N")
        >>> print my_seq
        NNNNNNNNNNN
        >>> my_protein = my_seq.translate()
        >>> my_protein
        UnknownSeq(3, alphabet = ProteinAlphabet(), character = 'X')
        >>> print my_protein
        XXX

        In comparison, using a normal Seq object:

        >>> my_seq = Seq("NNNNNNNNNNN")
        >>> print my_seq
        NNNNNNNNNNN
        >>> my_protein = my_seq.translate()
        >>> my_protein
        Seq('XXX', ExtendedIUPACProtein())
        >>> print my_protein
        XXX

        """
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet):
            raise ValueError("Proteins cannot be translated!")
        return UnknownSeq(self._length//3, Alphabet.generic_protein, "X")

    def ungap(self, gap=None):
        """Return a copy of the sequence without the gap character(s).

        The gap character can be specified in two ways - either as an explicit
        argument, or via the sequence's alphabet. For example:

        >>> from Bio.Seq import UnknownSeq
        >>> from Bio.Alphabet import Gapped, generic_dna
        >>> my_dna = UnknownSeq(20, Gapped(generic_dna,"-"))
        >>> my_dna
        UnknownSeq(20, alphabet = Gapped(DNAAlphabet(), '-'), character = 'N')
        >>> my_dna.ungap()
        UnknownSeq(20, alphabet = DNAAlphabet(), character = 'N')
        >>> my_dna.ungap("-")
        UnknownSeq(20, alphabet = DNAAlphabet(), character = 'N')

        If the UnknownSeq is using the gap character, then an empty Seq is
        returned:

        >>> my_gap = UnknownSeq(20, Gapped(generic_dna,"-"), character="-")
        >>> my_gap
        UnknownSeq(20, alphabet = Gapped(DNAAlphabet(), '-'), character = '-')
        >>> my_gap.ungap()
        Seq('', DNAAlphabet())
        >>> my_gap.ungap("-")
        Seq('', DNAAlphabet())

        Notice that the returned sequence's alphabet is adjusted to remove any
        explicit gap character declaration.
        """
        #Offload the alphabet stuff
        s = Seq(self._character, self.alphabet).ungap()
        if s :
            return UnknownSeq(self._length, s.alphabet, self._character)
        else :
            return Seq("", s.alphabet)

class MutableSeq(object):
    """An editable sequence object (with an alphabet).

    Unlike normal python strings and our basic sequence object (the Seq class)
    which are immuatable, the MutableSeq lets you edit the sequence in place.
    However, this means you cannot use a MutableSeq object as a dictionary key.

    >>> from Bio.Seq import MutableSeq
    >>> from Bio.Alphabet import generic_dna
    >>> my_seq = MutableSeq("ACTCGTCGTCG", generic_dna)
    >>> my_seq
    MutableSeq('ACTCGTCGTCG', DNAAlphabet())
    >>> my_seq[5]
    'T'
    >>> my_seq[5] = "A"
    >>> my_seq
    MutableSeq('ACTCGACGTCG', DNAAlphabet())
    >>> my_seq[5]
    'A'
    >>> my_seq[5:8] = "NNN"
    >>> my_seq
    MutableSeq('ACTCGNNNTCG', DNAAlphabet())
    >>> len(my_seq)
    11

    Note that the MutableSeq object does not support as many string-like
    or biological methods as the Seq object.
    """
    def __init__(self, data, alphabet = Alphabet.generic_alphabet):
        if sys.version_info[0] == 3:
            self.array_indicator = "u"
        else:
            self.array_indicator = "c"
        if isinstance(data, str): #TODO - What about unicode?
            self.data = array.array(self.array_indicator, data)
        else:
            self.data = data   # assumes the input is an array
        self.alphabet = alphabet
    
    def __repr__(self):
        """Returns a (truncated) representation of the sequence for debugging."""
        if len(self) > 60:
            #Shows the last three letters as it is often useful to see if there
            #is a stop codon at the end of a sequence.
            #Note total length is 54+3+3=60
            return "%s('%s...%s', %s)" % (self.__class__.__name__,
                                   str(self[:54]), str(self[-3:]),
                                   repr(self.alphabet))
        else:
            return "%s('%s', %s)" % (self.__class__.__name__,
                                   str(self),
                                   repr(self.alphabet))

    def __str__(self):
        """Returns the full sequence as a python string.

        Note that Biopython 1.44 and earlier would give a truncated
        version of repr(my_seq) for str(my_seq).  If you are writing code
        which needs to be backwards compatible with old Biopython, you
        should continue to use my_seq.tostring() rather than str(my_seq).
        """
        #See test_GAQueens.py for an historic usage of a non-string alphabet!
        return "".join(self.data)

    def __cmp__(self, other):
        """Compare the sequence to another sequence or a string (README).

        Currently if compared to another sequence the alphabets must be
        compatible. Comparing DNA to RNA, or Nucleotide to Protein will raise
        an exception. Otherwise only the sequence itself is compared, not the
        precise alphabet.

        A future release of Biopython will change this (and the Seq object etc)
        to use simple string comparison. The plan is that comparing sequences
        with incompatible alphabets (e.g. DNA to RNA) will trigger a warning
        but not an exception.

        During this transition period, please just do explicit comparisons:

        >>> seq1 = MutableSeq("ACGT")
        >>> seq2 = MutableSeq("ACGT")
        >>> id(seq1) == id(seq2)
        False
        >>> str(seq1) == str(seq2)
        True

        This method indirectly supports ==, < , etc.
        """
        if hasattr(other, "alphabet"):
            #other should be a Seq or a MutableSeq
            import warnings
            warnings.warn("In future comparing incompatible alphabets will "
                          "only trigger a warning (not an exception). In " 
                          "the interim please use id(seq1)==id(seq2) or "
                          "str(seq1)==str(seq2) to make your code explicit "
                          "and to avoid this warning.", FutureWarning)
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]):
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            if isinstance(other, MutableSeq):
                #See test_GAQueens.py for an historic usage of a non-string
                #alphabet!  Comparing the arrays supports this.
                return cmp(self.data, other.data)
            else:
                return cmp(str(self), str(other))
        elif isinstance(other, basestring):
            return cmp(str(self), other)
        else:
            raise TypeError

    def __len__(self): return len(self.data)

    def __getitem__(self, index):
        #Note since Python 2.0, __getslice__ is deprecated
        #and __getitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int):
            #Return a single letter as a string
            return self.data[index]
        else:
            #Return the (sub)sequence as another Seq object
            return MutableSeq(self.data[index], self.alphabet)

    def __setitem__(self, index, value):
        #Note since Python 2.0, __setslice__ is deprecated
        #and __setitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int):
            #Replacing a single letter with a new string
            self.data[index] = value
        else:
            #Replacing a sub-sequence
            if isinstance(value, MutableSeq):
                self.data[index] = value.data
            elif isinstance(value, type(self.data)):
                self.data[index] = value
            else:
                self.data[index] = array.array(self.array_indicator,
                                               str(value))

    def __delitem__(self, index):
        #Note since Python 2.0, __delslice__ is deprecated
        #and __delitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        
        #Could be deleting a single letter, or a slice
        del self.data[index]
    
    def __add__(self, other):
        """Add another sequence or string to this sequence.

        Returns a new MutableSeq object."""
        if hasattr(other, "alphabet"):
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]):
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            if isinstance(other, MutableSeq):
                #See test_GAQueens.py for an historic usage of a non-string
                #alphabet!  Adding the arrays should support this.
                return self.__class__(self.data + other.data, a)
            else:
                return self.__class__(str(self) + str(other), a)
        elif isinstance(other, basestring):
            #other is a plain string - use the current alphabet
            return self.__class__(str(self) + str(other), self.alphabet)
        else:
            raise TypeError

    def __radd__(self, other):
        if hasattr(other, "alphabet"):
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]):
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            if isinstance(other, MutableSeq):
                #See test_GAQueens.py for an historic usage of a non-string
                #alphabet!  Adding the arrays should support this.
                return self.__class__(other.data + self.data, a)
            else:
                return self.__class__(str(other) + str(self), a)
        elif isinstance(other, basestring):
            #other is a plain string - use the current alphabet
            return self.__class__(str(other) + str(self), self.alphabet)
        else:
            raise TypeError

    def append(self, c):
        self.data.append(c)

    def insert(self, i, c):
        self.data.insert(i, c)

    def pop(self, i = (-1)):
        c = self.data[i]
        del self.data[i]
        return c

    def remove(self, item):
        for i in range(len(self.data)):
            if self.data[i] == item:
                del self.data[i]
                return
        raise ValueError("MutableSeq.remove(x): x not in list")

    def count(self, sub, start=0, end=sys.maxint):
        """Non-overlapping count method, like that of a python string.

        This behaves like the python string method of the same name,
        which does a non-overlapping count!

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.
    
        Arguments:
         - sub - a string or another Seq object to look for
         - start - optional integer, slice start
         - end - optional integer, slice end

        e.g.
        
        >>> from Bio.Seq import MutableSeq
        >>> my_mseq = MutableSeq("AAAATGA")
        >>> print my_mseq.count("A")
        5
        >>> print my_mseq.count("ATG")
        1
        >>> print my_mseq.count(Seq("AT"))
        1
        >>> print my_mseq.count("AT", 2, -1)
        1
        
        HOWEVER, please note because that python strings, Seq objects and
        MutableSeq objects do a non-overlapping search, this may not give
        the answer you expect:

        >>> "AAAA".count("AA")
        2
        >>> print MutableSeq("AAAA").count("AA")
        2

        A non-overlapping search would give the answer as three!
        """
        try:
            #TODO - Should we check the alphabet?
            search = sub.tostring()
        except AttributeError:
            search = sub

        if not isinstance(search, basestring):
            raise TypeError("expected a string, Seq or MutableSeq")

        if len(search) == 1:
            #Try and be efficient and work directly from the array.
            count = 0
            for c in self.data[start:end]:
                if c == search: count += 1
            return count
        else:
            #TODO - Can we do this more efficiently?
            return self.tostring().count(search, start, end)

    def index(self, item):
        for i in range(len(self.data)):
            if self.data[i] == item:
                return i
        raise ValueError("MutableSeq.index(x): x not in list")

    def reverse(self):
        """Modify the mutable sequence to reverse itself.

        No return value.
        """
        self.data.reverse()

    def complement(self):
        """Modify the mutable sequence to take on its complement.

        Trying to complement a protein sequence raises an exception.

        No return value.
        """
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet):
            raise ValueError("Proteins do not have complements!")
        if self.alphabet in (IUPAC.ambiguous_dna, IUPAC.unambiguous_dna):
            d = ambiguous_dna_complement
        elif self.alphabet in (IUPAC.ambiguous_rna, IUPAC.unambiguous_rna):
            d = ambiguous_rna_complement
        elif 'U' in self.data and 'T' in self.data:
            #TODO - Handle this cleanly?
            raise ValueError("Mixed RNA/DNA found")
        elif 'U' in self.data:
            d = ambiguous_rna_complement
        else:
            d = ambiguous_dna_complement
        c = dict([(x.lower(), y.lower()) for x,y in d.iteritems()])
        d.update(c)
        self.data = map(lambda c: d[c], self.data)
        self.data = array.array(self.array_indicator, self.data)
        
    def reverse_complement(self):
        """Modify the mutable sequence to take on its reverse complement.

        Trying to reverse complement a protein sequence raises an exception.

        No return value.
        """
        self.complement()
        self.data.reverse()

    ## Sorting a sequence makes no sense.
    # def sort(self, *args): self.data.sort(*args)
    
    def extend(self, other):
        if isinstance(other, MutableSeq):
            for c in other.data:
                self.data.append(c)
        else:
            for c in other:
                self.data.append(c)

    def tostring(self):
        """Returns the full sequence as a python string (semi-obsolete).

        Although not formally deprecated, you are now encouraged to use
        str(my_seq) instead of my_seq.tostring().

        Because str(my_seq) will give you the full sequence as a python string,
        there is often no need to make an explicit conversion.  For example,
        
        print "ID={%s}, sequence={%s}" % (my_name, my_seq)

        On Biopython 1.44 or older you would have to have done this:

        print "ID={%s}, sequence={%s}" % (my_name, my_seq.tostring())
        """
        return "".join(self.data)

    def toseq(self):
        """Returns the full sequence as a new immutable Seq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_mseq = MutableSeq("MKQHKAMIVALIVICITAVVAAL", 
        ...                      IUPAC.protein)
        >>> my_mseq
        MutableSeq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())
        >>> my_mseq.toseq()
        Seq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())

        Note that the alphabet is preserved.
        """
        return Seq("".join(self.data), self.alphabet)

# The transcribe, backward_transcribe, and translate functions are
# user-friendly versions of the corresponding functions in Bio.Transcribe
# and Bio.Translate. The functions work both on Seq objects, and on strings.

def transcribe(dna):
    """Transcribes a DNA sequence into RNA.

    If given a string, returns a new string object.

    Given a Seq or MutableSeq, returns a new Seq object with an RNA alphabet.

    Trying to transcribe a protein or RNA sequence raises an exception.

    e.g.
    
    >>> transcribe("ACTGN")
    'ACUGN'
    """
    if isinstance(dna, Seq):
        return dna.transcribe()
    elif isinstance(dna, MutableSeq):
        return dna.toseq().transcribe()
    else:
        return dna.replace('T','U').replace('t','u')

def back_transcribe(rna):
    """Back-transcribes an RNA sequence into DNA.

    If given a string, returns a new string object.
    
    Given a Seq or MutableSeq, returns a new Seq object with an RNA alphabet.

    Trying to transcribe a protein or DNA sequence raises an exception.

    e.g.

    >>> back_transcribe("ACUGN")
    'ACTGN'
    """
    if isinstance(rna, Seq):
        return rna.back_transcribe()
    elif isinstance(rna, MutableSeq):
        return rna.toseq().back_transcribe()
    else:
        return rna.replace('U','T').replace('u','t')
    
def _translate_str(sequence, table, stop_symbol="*", to_stop=False,
                   cds=False, pos_stop="X"):
    """Helper function to translate a nucleotide string (PRIVATE).

    Arguments:
     - sequence    - a string
     - table       - a CodonTable object (NOT a table name or id number)
     - stop_symbol - a single character string, what to use for terminators.
     - to_stop     - boolean, should translation terminate at the first
                     in frame stop codon?  If there is no in-frame stop codon
                     then translation continues to the end.
     - pos_stop    - a single character string for a possible stop codon
                     (e.g. TAN or NNN)
     - cds - Boolean, indicates this is a complete CDS.  If True, this
             checks the sequence starts with a valid alternative start
             codon (which will be translated as methionine, M), that the
             sequence length is a multiple of three, and that there is a
             single in frame stop codon at the end (this will be excluded
             from the protein sequence, regardless of the to_stop option).
             If these tests fail, an exception is raised.

    Returns a string.

    e.g.

    >>> from Bio.Data import CodonTable
    >>> table = CodonTable.ambiguous_dna_by_id[1]
    >>> _translate_str("AAA", table)
    'K'
    >>> _translate_str("TAR", table)
    '*'
    >>> _translate_str("TAN", table)
    'X'
    >>> _translate_str("TAN", table, pos_stop="@")
    '@'
    >>> _translate_str("TA?", table)
    Traceback (most recent call last):
       ...
    TranslationError: Codon 'TA?' is invalid
    >>> _translate_str("ATGCCCTAG", table, cds=True)
    'MP'
    >>> _translate_str("AAACCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    TranslationError: First codon 'AAA' is not a start codon
    >>> _translate_str("ATGCCCTAGCCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    TranslationError: Extra in frame stop codon found.
    """
    sequence = sequence.upper()
    amino_acids = []
    forward_table = table.forward_table
    stop_codons = table.stop_codons
    if table.nucleotide_alphabet.letters is not None:
        valid_letters = set(table.nucleotide_alphabet.letters.upper())
    else:
        #Assume the worst case, ambiguous DNA or RNA:
        valid_letters = set(IUPAC.ambiguous_dna.letters.upper() + \
                            IUPAC.ambiguous_rna.letters.upper())
    if cds:
        if str(sequence[:3]).upper() not in table.start_codons:
            raise CodonTable.TranslationError(\
                "First codon '%s' is not a start codon" % sequence[:3])
        if len(sequence) % 3 != 0:
            raise CodonTable.TranslationError(\
                "Sequence length %i is not a multiple of three" % len(sequence))
        if str(sequence[-3:]).upper() not in stop_codons:
            raise CodonTable.TranslationError(\
                "Final codon '%s' is not a stop codon" % sequence[-3:])
        #Don't translate the stop symbol, and manually translate the M
        sequence = sequence[3:-3]
        amino_acids = ["M"]
    n = len(sequence)
    for i in xrange(0,n-n%3,3):
        codon = sequence[i:i+3]
        try:
            amino_acids.append(forward_table[codon])
        except (KeyError, CodonTable.TranslationError):
            #Todo? Treat "---" as a special case (gapped translation)
            if codon in table.stop_codons:
                if cds:
                    raise CodonTable.TranslationError(\
                        "Extra in frame stop codon found.")
                if to_stop : break
                amino_acids.append(stop_symbol)
            elif valid_letters.issuperset(set(codon)):
                #Possible stop codon (e.g. NNN or TAN)
                amino_acids.append(pos_stop)
            else:
                raise CodonTable.TranslationError(\
                    "Codon '%s' is invalid" % codon)
    return "".join(amino_acids)

def translate(sequence, table="Standard", stop_symbol="*", to_stop=False,
              cds=False):
    """Translate a nucleotide sequence into amino acids.

    If given a string, returns a new string object. Given a Seq or
    MutableSeq, returns a Seq object with a protein alphabet.

    Arguments:
     - table - Which codon table to use?  This can be either a name (string),
               an NCBI identifier (integer), or a CodonTable object (useful
               for non-standard genetic codes).  Defaults to the "Standard"
               table.
     - stop_symbol - Single character string, what to use for any
                     terminators, defaults to the asterisk, "*".
     - to_stop - Boolean, defaults to False meaning do a full
                 translation continuing on past any stop codons
                 (translated as the specified stop_symbol).  If
                 True, translation is terminated at the first in
                 frame stop codon (and the stop_symbol is not
                 appended to the returned protein sequence).
     - cds - Boolean, indicates this is a complete CDS.  If True, this
                 checks the sequence starts with a valid alternative start
                 codon (which will be translated as methionine, M), that the
                 sequence length is a multiple of three, and that there is a
                 single in frame stop codon at the end (this will be excluded
                 from the protein sequence, regardless of the to_stop option).
                 If these tests fail, an exception is raised.
    
    A simple string example using the default (standard) genetic code:
    
    >>> coding_dna = "GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    >>> translate(coding_dna)
    'VAIVMGR*KGAR*'
    >>> translate(coding_dna, stop_symbol="@")
    'VAIVMGR@KGAR@'
    >>> translate(coding_dna, to_stop=True)
    'VAIVMGR'
     
    Now using NCBI table 2, where TGA is not a stop codon:

    >>> translate(coding_dna, table=2)
    'VAIVMGRWKGAR*'
    >>> translate(coding_dna, table=2, to_stop=True)
    'VAIVMGRWKGAR'

    In fact this example uses an alternative start codon valid under NCBI table 2,
    GTG, which means this example is a complete valid CDS which when translated
    should really start with methionine (not valine):
    
    >>> translate(coding_dna, table=2, cds=True)
    'MAIVMGRWKGAR'

    Note that if the sequence has no in-frame stop codon, then the to_stop
    argument has no effect:

    >>> coding_dna2 = "GTGGCCATTGTAATGGGCCGC"
    >>> translate(coding_dna2)
    'VAIVMGR'
    >>> translate(coding_dna2, to_stop=True)
    'VAIVMGR'
    
    NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
    or a stop codon.  These are translated as "X".  Any invalid codon
    (e.g. "TA?" or "T-A") will throw a TranslationError.

    NOTE - Does NOT support gapped sequences.
    
    It will however translate either DNA or RNA.
    """
    if isinstance(sequence, Seq):
        return sequence.translate(table, stop_symbol, to_stop, cds)
    elif isinstance(sequence, MutableSeq):
        #Return a Seq object
        return sequence.toseq().translate(table, stop_symbol, to_stop, cds)
    else:
        #Assume its a string, return a string
        try:
            codon_table = CodonTable.ambiguous_generic_by_id[int(table)]
        except ValueError:
            codon_table = CodonTable.ambiguous_generic_by_name[table]
        except (AttributeError, TypeError):
            if isinstance(table, CodonTable.CodonTable):
                codon_table = table
            else:
                raise ValueError('Bad table argument')
        return _translate_str(sequence, codon_table, stop_symbol, to_stop, cds)
      
def reverse_complement(sequence):
    """Returns the reverse complement sequence of a nucleotide string.

    If given a string, returns a new string object.
    Given a Seq or a MutableSeq, returns a new Seq object with the same alphabet.

    Supports unambiguous and ambiguous nucleotide sequences.

    e.g.

    >>> reverse_complement("ACTG-NH")
    'DN-CAGT'
    """
    if isinstance(sequence, Seq):
        #Return a Seq
        return sequence.reverse_complement()
    elif isinstance(sequence, MutableSeq):
        #Return a Seq
        #Don't use the MutableSeq reverse_complement method as it is 'in place'.
        return sequence.toseq().reverse_complement()

    #Assume its a string.
    #In order to avoid some code duplication, the old code would turn the string
    #into a Seq, use the reverse_complement method, and convert back to a string.
    #This worked, but is over five times slower on short sequences!
    if ('U' in sequence or 'u' in sequence) \
    and ('T' in sequence or 't' in sequence):
        raise ValueError("Mixed RNA/DNA found")
    elif 'U' in sequence or 'u' in sequence:
        ttable = _rna_complement_table
    else:
        ttable = _dna_complement_table
    return sequence.translate(ttable)[::-1]

def _test():
    """Run the Bio.Seq module's doctests (PRIVATE)."""
    if sys.version_info[0:2] == (3,1):
        print "Not running Bio.Seq doctest on Python 3.1"
        print "See http://bugs.python.org/issue7490"
    else:
        print "Runing doctests..."
        import doctest
        doctest.testmod(optionflags=doctest.IGNORE_EXCEPTION_DETAIL)
        print "Done"

if __name__ == "__main__":
    _test()
