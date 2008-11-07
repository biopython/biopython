# Copyright 2000-2002 Brad Chapman.
# Copyright 2004-2005 by M de Hoon.
# Copyright 2007-2008 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Represent a sequence or mutable sequence, with an alphabet."""

import string, array
import sys

#TODO - Remove this work around once we drop python 2.3 support
try:
    set = set
except NameError:
    from sets import Set as set

import Alphabet
from Alphabet import IUPAC
from Data.IUPACData import ambiguous_dna_complement, ambiguous_rna_complement
from Bio.Data import CodonTable

class Seq(object):
    """A read-only sequence object (essentially a string with an alphabet).

    Like normal python strings, our basic sequence object is immutable.
    This prevents you from doing my_seq[5] = "A" for example, but does allow
    Seq objects to be used as dictionary keys.

    The Seq object provides a number of string like methods (such as count,
    find, split and strip), which are alphabet aware where appropriate.

    The Seq object also provides some biological methods, such as complement,
    reverse_complement, transcribe, back_transcribe and translate (which are
    not applicable to sequences with a protein alphabet).
    """
    def __init__(self, data, alphabet = Alphabet.generic_alphabet):
        """Create a Seq object.

        Arguments:
        seq      - Sequence, required (string)
        alphabet - Optional argument, an Alphabet object from Bio.Alphabet
        
        You will typically use Bio.SeqIO to read in sequences from files as
        SeqRecord objects, whose sequence will be exposed as a Seq object via
        the seq property.

        However, will often want to create your own Seq objects directly:

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF", \
                          IUPAC.protein)
        >>> my_seq
        Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())
        >>> print my_seq
        MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
        """
        # Enforce string storage
        assert (type(data) == type("") or # must use a string
                type(data) == type(u""))  # but can be a unicode string

        self._data = data
        self.alphabet = alphabet  # Seq API requirement
 
    # A data property is/was a Seq API requirement
    def _set_data(self, value) :
        #TODO - In the next release, actually raise an exception?
        #The Seq object is like a python string, it should be read only!
        import warnings
        warnings.warn("Writing to the Seq object's .data propery is deprecated.",
                      DeprecationWarning)
        self._data = value
    data = property(fget= lambda self : self._data,
                    fset=_set_data,
                    doc="Sequence as a string (DEPRECATED)")

    def __repr__(self):
        """Returns a (truncated) representation of the sequence for debugging."""
        if len(self) > 60 :
            #Shows the last three letters as it is often useful to see if there
            #is a stop codon at the end of a sequence.
            #Note total length is 54+3+3=60
            return "%s('%s...%s', %s)" % (self.__class__.__name__,
                                   str(self)[:54], str(self)[-3:],
                                   repr(self.alphabet))
        else :
            return "%s(%s, %s)" % (self.__class__.__name__,
                                  repr(self.data),
                                   repr(self.alphabet))
    def __str__(self):
        """Returns the full sequence as a python string.

        Note that Biopython 1.44 and earlier would give a truncated
        version of repr(my_seq) for str(my_seq).  If you are writing code
        which need to be backwards compatible with old Biopython, you
        should continue to use my_seq.tostring() rather than str(my_seq)
        """
        return self._data

    """
    TODO - Work out why this breaks test_Restriction.py
    (Comparing Seq objects would be nice to have.  May need to think about
    hashes and the in operator for when have list/dictionary of Seq objects...)
    def __cmp__(self, other):
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            return cmp(str(self), str(other))
        elif isinstance(other, basestring) :
            return cmp(str(self), other)
        else :
            raise TypeError
    """

    def __len__(self): return len(self._data)       # Seq API requirement

    def __getitem__(self, index) :                 # Seq API requirement
        #Note since Python 2.0, __getslice__ is deprecated
        #and __getitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int) :
            #Return a single letter as a string
            return self._data[index]
        else :
            #Return the (sub)sequence as another Seq object
            return Seq(self._data[index], self.alphabet)

    def __add__(self, other):
        """Add another sequence or string to this sequence."""
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            return self.__class__(str(self) + str(other), a)
        elif isinstance(other, basestring) :
            #other is a plain string - use the current alphabet
            return self.__class__(str(self) + other, self.alphabet)
        else :
            raise TypeError

    def __radd__(self, other):
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            return self.__class__(str(other) + str(self), a)
        elif isinstance(other, basestring) :
            #other is a plain string - use the current alphabet
            return self.__class__(other + str(self), self.alphabet)
        else :
            raise TypeError

    def tostring(self):                            # Seq API requirement
        """Returns the full sequence as a python string.

        Although not formally deprecated, you are now encouraged to use
        str(my_seq) instead of my_seq.tostring()."""
        return self._data

    def tomutable(self):   # Needed?  Or use a function?
        """Returns the full sequence as a MutableSeq object.

        >>> from Bio.Seq import Seq
        >>> from Bio.Alphabet import IUPAC
        >>> my_seq = Seq("MKQHKAMIVALIVICITAVVAAL", \
                          IUPAC.protein)
        >>> my_seq
        Seq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())
        >>> my_seq.tomutable()
        MutableSeq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())

        Note that the alphabet is preserved.
        """
        return MutableSeq(str(self), self.alphabet)

    def _get_seq_str_and_check_alphabet(self, other_sequence) :
        """string/Seq/MutableSeq to string, checking alphabet (PRIVATE).

        For a string argument, returns the string.

        For a Seq or MutableSeq, it checks the alphabet is compatible
        (raising an exception if it isn't), and then returns a string.
        """
        try :
            other_alpha = other_sequence.alphabet
        except AttributeError :
            #Assume other_sequence is a string
            return other_sequence

        #Other should be a Seq or a MutableSeq
        if not Alphabet._check_type_compatible([self.alphabet, other_alpha]) :
            raise TypeError("Incompatable alphabets %s and %s" \
                            % (repr(self.alphabet), repr(other_alpha)))
        #Return as a string
        return str(other_sequence)
    
    def count(self, sub, start=0, end=sys.maxint):
        """Count method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.
    
        sub - a string or another Seq object to look for
        start - optional integer, slice start
        end - optional integer, slice end

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
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return str(self).count(sub_str, start, end)

    def find(self, sub, start=0, end=sys.maxint):
        """Find method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        sub - a string or another Seq object to look for
        start - optional integer, slice start
        end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return str(self).find(sub_str, start, end)

    def rfind(self, sub, start=0, end=sys.maxint):
        """Find from right method, like that of a python string.

        This behaves like the python string method of the same name.

        Returns an integer, the index of the last (right most) occurrence of
        substring argument sub in the (sub)sequence given by [start:end].

        sub - a string or another Seq object to look for
        start - optional integer, slice start
        end - optional integer, slice end

        Returns -1 if the subsequence is NOT found.
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return str(self).rfind(sub_str, start, end)
    
    def split(self, sep=None, maxsplit=-1) :
        """Split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done.  If maxsplit is ommited, all
        splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.
        
        e.g. print my_seq.split("*")

        See also the rsplit method.
        """
        #If it has one, check the alphabet:
        sep_str = self._get_seq_str_and_check_alphabet(sep)
        return [Seq(part, self.alphabet) \
                for part in str(self).split(sep_str, maxsplit)]

    def rsplit(self, sep=None, maxsplit=-1) :
        """Right split method, like that of a python string.

        This behaves like the python string method of the same name.

        Return a list of the 'words' in the string (as Seq objects),
        using sep as the delimiter string.  If maxsplit is given, at
        most maxsplit splits are done COUNTING FROM THE RIGHT.
        If maxsplit is ommited, all splits are made.

        Following the python string method, sep will by default be any
        white space (tabs, spaces, newlines) but this is unlikely to
        apply to biological sequences.
        
        e.g. print my_seq.split("*")

        See also the split method.
        """
        #If it has one, check the alphabet:
        sep_str = self._get_seq_str_and_check_alphabet(sep)
        try :
            return [Seq(part, self.alphabet) \
                    for part in str(self).rsplit(sep_str, maxsplit)]
        except AttributeError :
            #Python 2.3 doesn't have a string rsplit method, which we can
            #word around by reversing the sequence, using (left) split,
            #and then reversing the answer. Not very efficient!
            words = [Seq(word[::-1], self.alphabet) for word \
                     in str(self)[::-1].split(sep_str[::-1], maxsplit)]
            words.reverse()
            return words

    def strip(self, chars=None) :
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

    def lstrip(self, chars=None) :
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

    def rstrip(self, chars=None) :
        """Returns a new Seq object with trailing (right) end stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        ommitted or None (default) then as for the python string method,
        this defaults to removing any white space.
        
        e.g. print my_seq.lstrip("-")

        See also the strip and lstrip methods.
        """
        #If it has one, check the alphabet:
        strip_str = self._get_seq_str_and_check_alphabet(chars)
        return Seq(str(self).rstrip(strip_str), self.alphabet)

    def __maketrans(self, alphabet) :
        """Seq.__maketrans(alphabet) -> translation table (PRIVATE).

        Return a translation table for use with complement()
        and reverse_complement().

        Compatible with lower case and upper case sequences.

        alphabet is a dictionary as implement in Data.IUPACData

        For internal use only.
        """
        before = ''.join(alphabet.keys())
        after  = ''.join(alphabet.values())
        before = before + before.lower()
        after  = after + after.lower()
        return string.maketrans(before, after)

    def complement(self):
        """Returns the complement sequence. New Seq object.

        Trying to complement a protein sequence raises an exception.
        """
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet) :
            raise ValueError("Proteins do not have complements!")
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.DNAAlphabet) :
            d = ambiguous_dna_complement
        elif isinstance(Alphabet._get_base_alphabet(self.alphabet),
                        Alphabet.RNAAlphabet) :
            d = ambiguous_rna_complement
        elif 'U' in self._data and 'T' in self._data:
            #TODO - Handle this cleanly?
            raise ValueError("Mixed RNA/DNA found")
        elif 'U' in self._data:
            d = ambiguous_rna_complement
        else:
            d = ambiguous_dna_complement
        ttable = self.__maketrans(d)
        #Much faster on really long sequences than the previous loop based one.
        #thx to Michael Palmer, University of Waterloo
        s = str(self).translate(ttable)
        return Seq(s, self.alphabet)

    def reverse_complement(self):
        """Returns the reverse complement sequence. New Seq object.

        Trying to complement a protein sequence raises an exception.
        """
        #Use -1 stride/step to reverse the complement
        return self.complement()[::-1]

    def transcribe(self):
        """Returns the RNA sequence from a DNA sequence. New Seq object.

        Trying to transcribe a protein or RNA sequence raises an exception.
        """
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet) :
            raise ValueError("Proteins cannot be transcribed!")
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.RNAAlphabet) :
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

        Trying to back-transcribe a protein or DNA sequence raises an
        exception.
        """
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet) :
            raise ValueError("Proteins cannot be back transcribed!")
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.DNAAlphabet) :
            raise ValueError("DNA cannot be back transcribed!")

        if self.alphabet==IUPAC.unambiguous_rna:
            alphabet = IUPAC.unambiguous_dna
        elif self.alphabet==IUPAC.ambiguous_rna:
            alphabet = IUPAC.ambiguous_dna
        else:
            alphabet = Alphabet.generic_dna
        return Seq(str(self).replace("U", "T").replace("u", "t"), alphabet)

    def translate(self, table="Standard", stop_symbol="*", to_stop=False):
        """Turns a nucleotide sequence into a protein sequence. New Seq object.

        Trying to back-transcribe a protein sequence raises an exception.
        This method will translate DNA or RNA sequences.

        Trying to translate a protein sequence raises an exception.

        table - Which codon table to use?  This can be either a name
                (string) or an NCBI identifier (integer).  This defaults
                to the "Standard" table.
        stop_symbol - Single character string, what to use for terminators.
                This defaults to the asterisk, "*".
        to_stop - Boolean, defaults to False meaning do a full translation
                continuing on past any stop codons (translated as the
                specified stop_symbol).  If True, translation is terminated
                at the first in frame stop codon (and the stop_symbol is
                not appended to the returned protein sequence).

        e.g. Using the standard table,
        
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
        try:
            table_id = int(table)
        except ValueError:
            table_id = None
        if isinstance(table, str) and len(table)==256 :
            raise ValueError("The Seq object translate method DOES NOT take " \
                             + "a 256 character string mapping table like " \
                             + "the python string object's translate method. " \
                             + "Use str(my_seq).translate(...) instead.")
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet) :
            raise ValueError("Proteins cannot be translated!")
        if self.alphabet==IUPAC.unambiguous_dna:
            if table_id is None:
                codon_table = CodonTable.unambiguous_dna_by_name[table]
            else:
                codon_table = CodonTable.unambiguous_dna_by_id[table_id]
        elif self.alphabet==IUPAC.ambiguous_dna:
            if table_id is None:
                codon_table = CodonTable.ambiguous_dna_by_name[table]
            else:
                codon_table = CodonTable.ambiguous_dna_by_id[table_id]
        elif self.alphabet==IUPAC.unambiguous_rna:
            if table_id is None:
                codon_table = CodonTable.unambiguous_rna_by_name[table]
            else:
                codon_table = CodonTable.unambiguous_rna_by_id[table_id]
        elif self.alphabet==IUPAC.ambiguous_rna:
            if table_id is None:
                codon_table = CodonTable.ambiguous_rna_by_name[table]
            else:
                codon_table = CodonTable.ambiguous_rna_by_id[table_id]
        else:
            if table_id is None:
                codon_table = CodonTable.ambiguous_generic_by_name[table]
            else:
                codon_table = CodonTable.ambiguous_generic_by_id[table_id]
        protein = _translate_str(str(self), codon_table, stop_symbol, to_stop)
        if stop_symbol in protein :
            alphabet = Alphabet.HasStopCodon(codon_table.protein_alphabet,
                                             stop_symbol = stop_symbol)
        else :
            alphabet = codon_table.protein_alphabet
        return Seq(protein, alphabet)

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
        if type(data) == type(""):
            self.data = array.array("c", data)
        else:
            self.data = data   # assumes the input is an array
        self.alphabet = alphabet
    
    def __repr__(self):
        """Returns a (truncated) representation of the sequence for debugging."""
        if len(self) > 60 :
            #Shows the last three letters as it is often useful to see if there
            #is a stop codon at the end of a sequence.
            #Note total length is 54+3+3=60
            return "%s('%s...%s', %s)" % (self.__class__.__name__,
                                   str(self[:54]), str(self[-3:]),
                                   repr(self.alphabet))
        else :
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
        """Compare the sequence for to another sequence or a string.

        If compared to another sequence the alphabets must be compatible.
        Comparing DNA to RNA, or Nucleotide to Protein will raise an
        exception.

        Otherwise only the sequence itself is compared, not the precise
        alphabet.

        This method indirectly supports ==, < , etc."""
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            if isinstance(other, MutableSeq):
                #See test_GAQueens.py for an historic usage of a non-string
                #alphabet!  Comparing the arrays supports this.
                return cmp(self.data, other.data)
            else :
                return cmp(str(self), str(other))
        elif isinstance(other, basestring) :
            return cmp(str(self), other)
        else :
            raise TypeError

    def __len__(self): return len(self.data)

    def __getitem__(self, index) :
        #Note since Python 2.0, __getslice__ is deprecated
        #and __getitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int) :
            #Return a single letter as a string
            return self.data[index]
        else :
            #Return the (sub)sequence as another Seq object
            return MutableSeq(self.data[index], self.alphabet)

    def __setitem__(self, index, value):
        #Note since Python 2.0, __setslice__ is deprecated
        #and __setitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int) :
            #Replacing a single letter with a new string
            self.data[index] = value
        else :
            #Replacing a sub-sequence
            if isinstance(value, MutableSeq):
                self.data[index] = value.data
            elif isinstance(value, type(self.data)):
                self.data[index] = value
            else:
                self.data[index] = array.array("c", str(value))

    def __delitem__(self, index):
        #Note since Python 2.0, __delslice__ is deprecated
        #and __delitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        
        #Could be deleting a single letter, or a slice
        del self.data[index]
    
    def __add__(self, other):
        """Add another sequence or string to this sequence.

        Returns a new MutableSeq object."""
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            if isinstance(other, MutableSeq):
                #See test_GAQueens.py for an historic usage of a non-string
                #alphabet!  Adding the arrays should support this.
                return self.__class__(self.data + other.data, a)
            else :
                return self.__class__(str(self) + str(other), a)
        elif isinstance(other, basestring) :
            #other is a plain string - use the current alphabet
            return self.__class__(str(self) + str(other), self.alphabet)
        else :
            raise TypeError

    def __radd__(self, other):
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError("Incompatable alphabets %s and %s" \
                                % (repr(self.alphabet), repr(other.alphabet)))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            if isinstance(other, MutableSeq):
                #See test_GAQueens.py for an historic usage of a non-string
                #alphabet!  Adding the arrays should support this.
                return self.__class__(other.data + self.data, a)
            else :
                return self.__class__(str(other) + str(self), a)
        elif isinstance(other, basestring) :
            #other is a plain string - use the current alphabet
            return self.__class__(str(other) + str(self), self.alphabet)
        else :
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
        """Count method, like that of a python string.

        Return an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.
    
        sub - a string or another Seq object to look for
        start - optional integer, slice start
        end - optional integer, slice end

        e.g.
        >>> from Bio.Seq import Seq, MutableSeq
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
        """
        try :
            #TODO - Should we check the alphabet?
            search = sub.tostring()
        except AttributeError :
            search = sub

        if len(search) == 1 :
            #Try and be efficient and work directly from the array.
            count = 0
            for c in self.data[start:end]:
                if c == search: count += 1
            return count
        else :
            #TODO - Can we do this more efficiently?
            return self.tostring().count(search, start, end)

    def index(self, item):
        for i in range(len(self.data)):
            if self.data[i] == item:
                return i
        raise ValueError("MutableSeq.index(x): x not in list")

    def reverse(self):
        """Modify the mutable sequence to reverse itself.

        No return value."""
        self.data.reverse()

    def complement(self):
        """Modify the mutable sequence to take on its complement.

        Trying to complement a protein sequence raises an exception.

        No return value"""
        if isinstance(Alphabet._get_base_alphabet(self.alphabet),
                      Alphabet.ProteinAlphabet) :
            raise ValueError("Proteins do not have complements!")
        if self.alphabet in (IUPAC.ambiguous_dna, IUPAC.unambiguous_dna):
            d = ambiguous_dna_complement
        elif self.alphabet in (IUPAC.ambiguous_rna, IUPAC.unambiguous_rna):
            d = ambiguous_rna_complement
        elif 'U' in self.data and 'T' in self.data :
            #TODO - Handle this cleanly?
            raise ValueError("Mixed RNA/DNA found")
        elif 'U' in self.data:
            d = ambiguous_rna_complement
        else:
            d = ambiguous_dna_complement
        c = dict([(x.lower(), y.lower()) for x,y in d.iteritems()])
        d.update(c)
        self.data = map(lambda c: d[c], self.data)
        self.data = array.array('c', self.data)
        
    def reverse_complement(self):
        """Modify the mutable sequence to take on its reverse complement.

        Trying to reverse complement a protein sequence raises an exception.

        No return value."""
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
        """Returns the full sequence as a python string.

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
        >>> my_mseq = MutableSeq("MKQHKAMIVALIVICITAVVAAL", \
                                 IUPAC.protein)
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
    if isinstance(dna, Seq) :
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
    if isinstance(rna, Seq) :
        return rna.back_transcribe()
    elif isinstance(rna, MutableSeq):
        return rna.toseq().back_transcribe()
    else:
        return rna.replace('U','T').replace('u','t')
    
def _translate_str(sequence, table, stop_symbol="*",
                   to_stop=False, pos_stop="X") :
    """Helper function to translate a nucleotide string (PRIVATE).

    sequence    - a string
    table       - a CodonTable object (NOT a table name or id number)
    stop_symbol - a single character string, what to use for terminators.
    to_stop     - boolean, should translation terminate at the first
                  in frame stop codon?  If there is no in-frame stop codon
                  then translation continues to the end.
    pos_stop    - a single character string for a possible stop codon
                  (e.g. TAN or NNN)

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
    """
    sequence = sequence.upper()
    amino_acids = []
    forward_table = table.forward_table
    stop_codons = table.stop_codons
    if table.nucleotide_alphabet.letters is not None :
        valid_letters = set(table.nucleotide_alphabet.letters.upper())
    else :
        #Assume the worst case, ambiguous DNA or RNA:
        valid_letters = set(IUPAC.ambiguous_dna.letters.upper() + \
                            IUPAC.ambiguous_rna.letters.upper())

    n = len(sequence)
    for i in xrange(0,n-n%3,3) :
        codon = sequence[i:i+3]
        try :
            amino_acids.append(forward_table[codon])
        except (KeyError, CodonTable.TranslationError) :
            #Todo? Treat "---" as a special case (gapped translation)
            if codon in table.stop_codons :
                if to_stop : break
                amino_acids.append(stop_symbol)
            elif valid_letters.issuperset(set(codon)) :
                #Possible stop codon (e.g. NNN or TAN)
                amino_acids.append(pos_stop)
            else :
                raise CodonTable.TranslationError(\
                    "Codon '%s' is invalid" % codon)
    return "".join(amino_acids)

def translate(sequence, table="Standard", stop_symbol="*", to_stop=False):
    """Translate a nucleotide sequence into amino acids.

    If given a string, returns a new string object.
    Given a Seq or MutableSeq, returns a Seq object with a protein
    alphabet.

    table - Which codon table to use?  This can be either a name (string) or
            an NCBI identifier (integer).  Defaults to the "Standard" table.
    stop_symbol - Single character string, what to use for terminators.
            This defaults to the asterisk, "*".
    to_stop - Boolean, defaults to False meaning do a full translation
            continuing on past any stop codons (translated as the
            specified stop_symbol).  If True, translation is terminated
            at the first in frame stop codon (and the stop_symbol is
            not appended to the returned protein sequence).

    A simple string example using the default (standard) genetic code,
    
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
    if isinstance(sequence, Seq) :
        return sequence.translate(table, stop_symbol, to_stop)
    elif isinstance(sequence, MutableSeq):
        #Return a Seq object
        return sequence.toseq().translate(table, stop_symbol, to_stop)
    else:
        #Assume its a string, return a string
        try :
            codon_table = CodonTable.ambiguous_generic_by_id[int(table)]
        except ValueError :
            codon_table = CodonTable.ambiguous_generic_by_name[table]
        return _translate_str(sequence, codon_table, stop_symbol, to_stop)
      
def reverse_complement(sequence):
    """Returns the reverse complement sequence of a nucleotide string.

    If given a string, returns a new string object.
    Given a Seq or a MutableSeq, returns a new Seq object with the same alphabet.

    Supports unambiguous and ambiguous nucleotide sequences.

    e.g.
    >>> reverse_complement("ACTGN")
    'NCAGT'
    """
    if isinstance(sequence, Seq) :
        #Return a Seq
        return sequence.reverse_complement()
    elif isinstance(sequence, MutableSeq) :
        #Return a Seq
        #Don't use the MutableSeq reverse_complement method as it is 'in place'.
        return sequence.toseq().reverse_complement()
    else :
        #Assume its a string, turn it into a Seq,
        #do the reverse complement, and turn this back to a string
        #TODO - Find a more efficient way to do this without code duplication?
        return Seq(sequence).reverse_complement().tostring()

def _test():
    """Run the Bio.Seq module's doctests."""
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
