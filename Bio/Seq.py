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

import Alphabet
from Alphabet import IUPAC
from Data.IUPACData import ambiguous_dna_complement, ambiguous_rna_complement
from Bio.Data import CodonTable

class Seq:
    """A read-only sequence object (essentially a string with an alphabet).

    Like normal python strings, our basic sequence object is immuatable.
    This prevents you from doing my_seq[5] = "A" for example, but does allow
    Seq objects to be used as dictionary keys."""
    def __init__(self, data, alphabet = Alphabet.generic_alphabet):
        # Enforce string storage
        assert (type(data) == type("") or # must use a string
                type(data) == type(u""))  # but can be a unicode string

        self.data = data                           # Seq API requirement
        self.alphabet = alphabet                   # Seq API requirement

    def __repr__(self):
        """Returns a (truncated) representation of the sequence for debugging."""
        if len(self) > 60 :
            #Shows the last three letters as it is often useful to see if there
            #is a stop codon at the end of a sequence.
            #Note total length is 54+3+3=60
            return "%s('%s...%s', %s)" % (self.__class__.__name__,
                                   self.data[:54], self.data[-3:],
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
        return self.data

    """
    TODO - Work out why this breaks test_Restriction.py
    def __cmp__(self, other):
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError, ("incompatable alphabets", str(self.alphabet),
                                  str(other.alphabet))
            #They should be the same sequence type (or one of them is generic)
            return cmp(str(self), str(other))
        elif isinstance(other, basestring) :
            return cmp(str(self), str(other))
        else :
            raise TypeError
    """

    def __len__(self): return len(self.data)       # Seq API requirement

    def __getitem__(self, index) :                 # Seq API requirement
        #Note since Python 2.0, __getslice__ is deprecated
        #and __getitem__ is used instead.
        #See http://docs.python.org/ref/sequence-methods.html
        if isinstance(index, int) :
            #Return a single letter as a string
            return self.data[index]
        else :
            #Return the (sub)sequence as another Seq object
            return Seq(self.data[index], self.alphabet)

    def __add__(self, other):
        """Add another sequence or string to this sequence."""
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError, ("incompatable alphabets", str(self.alphabet),
                                  str(other.alphabet))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
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
                raise TypeError, ("incompatable alphabets", str(self.alphabet),
                                  str(other.alphabet))
            #They should be the same sequence type (or one of them is generic)
            a = Alphabet._consensus_alphabet([self.alphabet, other.alphabet])
            return self.__class__(str(other) + str(self), a)
        elif isinstance(other, basestring) :
            #other is a plain string - use the current alphabet
            return self.__class__(str(other) + str(self), self.alphabet)
        else :
            raise TypeError

    def tostring(self):                            # Seq API requirement
        """Returns the full sequence as a python string.

        Although not formally deprecated, you are now encouraged to use
        str(my_seq) instead of my_seq.tostring()."""
        return self.data

    def tomutable(self):   # Needed?  Or use a function?
        return MutableSeq(self.data, self.alphabet)

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
            raise TypeError, ("incompatable alphabets", str(self.alphabet),
                              str(other_alpha))
        #Return as a string
        return other_sequence.tostring()
    
    def count(self, sub, start=0, end=sys.maxint):
        """Count method, like that of a python string.

        Returns an integer, the number of occurrences of substring
        argument sub in the (sub)sequence given by [start:end].
        Optional arguments start and end are interpreted as in slice
        notation.
    
        sub - a string or another Seq object to look for
        start - optional integer, slice start
        end - optional integer, slice end

        e.g.
        from Bio.Seq import Seq
        my_seq = Seq("AAAATGA")
        print my_seq.count("A")
        print my_seq.count("ATG")
        print my_seq.count(Seq("AT"))
        print my_seq.count("AT", 2, -1)
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return self.tostring().count(sub_str, start, end)

    def find(self, sub, start=0, end=sys.maxint):
        """Find method, like that of a python string.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        Returns -1 if the subsequence is NOT found.
        """
        #If it has one, check the alphabet:
        sub_str = self._get_seq_str_and_check_alphabet(sub)
        return self.tostring().find(sub_str, start, end)
    
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
        """
        #If it has one, check the alphabet:
        sep_str = self._get_seq_str_and_check_alphabet(sep)
        return [Seq(chunk, self.alphabet) \
                for chunk in str(self).split(sep_str, maxsplit)]

    def strip(self, chars=None) :
        """Returns a new Seq object with leading and trailing ends stripped.

        This behaves like the python string method of the same name.

        Optional argument chars defines which characters to remove.  If
        ommitted or None (default) then as for the python string method,
        this defaults to removing any white space.
        
        e.g. print my_seq.strip("-")
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
        """
        if isinstance(self.alphabet, Alphabet.ProteinAlphabet) :
            raise ValueError, "Proteins do not have complements!"
        if isinstance(self.alphabet, Alphabet.DNAAlphabet) :
            d = ambiguous_dna_complement
        elif isinstance(self.alphabet, Alphabet.RNAAlphabet) :
            d = ambiguous_rna_complement
        elif 'U' in self.data:
            d = ambiguous_rna_complement
        else:
            d = ambiguous_dna_complement
        ttable = self.__maketrans(d)
        #Much faster on really long sequences than the previous loop based one.
        #thx to Michael Palmer, University of Waterloo
        s = self.data.translate(ttable)
        return Seq(s, self.alphabet)

    def reverse_complement(self):
        """Returns the reverse complement sequence. New Seq object.
        """
        #Use -1 stride to reverse the complement
        return self.complement()[::-1]

class MutableSeq:
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
                raise TypeError, ("incompatable alphabets", str(self.alphabet),
                                  str(other.alphabet))
            #They should be the same sequence type (or one of them is generic)
            if isinstance(other, MutableSeq):
                #See test_GAQueens.py for an historic usage of a non-string
                #alphabet!  Comparing the arrays supports this.
                return cmp(self.data, other.data)
            else :
                return cmp(str(self), str(other))
        elif isinstance(other, basestring) :
            return cmp(str(self), str(other))
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
        """Add another sequence or string to this sequence."""
        if hasattr(other, "alphabet") :
            #other should be a Seq or a MutableSeq
            if not Alphabet._check_type_compatible([self.alphabet,
                                                    other.alphabet]) :
                raise TypeError, ("incompatable alphabets", str(self.alphabet),
                                  str(other.alphabet))
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
                raise TypeError, ("incompatable alphabets", str(self.alphabet),
                                  str(other.alphabet))
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
        raise ValueError, "MutableSeq.remove(x): x not in list"

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
        from Bio.Seq import MutableSeq
        my_mseq = MutableSeq("AAAATGA")
        print my_mseq.count("A")
        print my_mseq.count("ATG")
        print my_mseq.count(Seq("AT"))
        print my_mseq.count("AT", 2, -1)
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
        raise ValueError, "MutableSeq.index(x): x not in list"

    def reverse(self):
        """Modify the MutableSequence to reverse itself.

        No return value."""
        self.data.reverse()

    def complement(self):
        """Modify the MutableSequence to take on its complement.

        No return value"""
        if isinstance(self.alphabet, Alphabet.ProteinAlphabet) :
            raise ValueError, "Proteins do not have complements!"
        if self.alphabet in (IUPAC.ambiguous_dna, IUPAC.unambiguous_dna):
            d = ambiguous_dna_complement
        elif self.alphabet in (IUPAC.ambiguous_rna, IUPAC.unambiguous_rna):
            d = ambiguous_rna_complement
        elif 'U' in self.data:
            d = ambiguous_rna_complement
        else:
            d = ambiguous_dna_complement
        c = dict([(x.lower(), y.lower()) for x,y in d.iteritems()])
        d.update(c)
        self.data = map(lambda c: d[c], self.data)
        self.data = array.array('c', self.data)
        
    def reverse_complement(self):
        """Modify the MutableSequence to take on its reverse complement.

        No return value."""
        if isinstance(self.alphabet, Alphabet.ProteinAlphabet) :
            raise ValueError, "Proteins do not have complements!"
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
        """Returns the full sequence as a new immutable Seq object"""
        return Seq("".join(self.data), self.alphabet)


# The transcribe, backward_transcribe, and translate functions are
# user-friendly versions of the corresponding functions in Bio.Transcribe
# and Bio.Translate. The functions work both on Seq objects, and on strings.

def transcribe(dna):
    """Transcribes a DNA sequence into RNA.

    If given a string, returns a new string object.
    Given a Seq or MutableSeq, returns a new Seq object with the same alphabet.
    """
    if isinstance(dna, Seq) or isinstance(dna, MutableSeq):
        if isinstance(dna.alphabet, Alphabet.ProteinAlphabet) :
            raise ValueError, "Proteins cannot be transcribed!"
        if isinstance(dna.alphabet, Alphabet.RNAAlphabet) :
            raise ValueError, "RNA cannot be transcribed!"
        rna = dna.tostring().replace('T','U').replace('t','u')
        if dna.alphabet==IUPAC.unambiguous_dna:
            alphabet = IUPAC.unambiguous_rna
        elif dna.alphabet==IUPAC.ambiguous_dna:
            alphabet = IUPAC.ambiguous_rna
        else:
            alphabet = Alphabet.generic_rna
        return Seq(rna, alphabet)
    else:
        rna = dna.replace('T','U').replace('t','u')
        return rna


def back_transcribe(rna):
    """Back-transcribes an RNA sequence into DNA.

    If given a string, returns a new string object.
    Given a Seq or MutableSeq, returns a new Seq object with the same alphabet.
    """
    if isinstance(rna, Seq) or isinstance(rna, MutableSeq):
        if isinstance(rna.alphabet, Alphabet.ProteinAlphabet) :
            raise ValueError, "Proteins cannot be (back)transcribed!"
        if isinstance(rna.alphabet, Alphabet.DNAAlphabet) :
            raise ValueError, "DNA cannot be back transcribed!"
        dna = rna.tostring().replace('U','T').replace('u','t')
        if rna.alphabet==IUPAC.unambiguous_rna:
            alphabet = IUPAC.unambiguous_dna
        elif rna.alphabet==IUPAC.ambiguous_rna:
            alphabet = IUPAC.ambiguous_dna
        else:
            alphabet = Alphabet.generic_dna
        return Seq(dna, alphabet)
    else:
        dna = rna.replace('U','T').replace('u','t')
        return dna


def translate(sequence, table = "Standard", stop_symbol = "*"):
    """Translate a nucleotide sequence into amino acids.

    If given a string, returns a new string object.
    Given a Seq or MutableSeq, returns a Seq object.

    table - Which codon table to use?  This can be either a name
           (string) or an identifier (integer)

    NOTE - Does NOT support ambiguous nucleotide sequences which
    could code for a stop codon.  This will throw a TranslationError.

    NOTE - Does NOT support gapped sequences.
    
    It will however translate either DNA or RNA."""
    try:
        id = int(table)
    except:
        id = None
    if isinstance(sequence, Seq) or isinstance(sequence, MutableSeq):
        if isinstance(sequence.alphabet, Alphabet.ProteinAlphabet) :
            raise ValueError, "Proteins cannot be translated!"
        if sequence.alphabet==IUPAC.unambiguous_dna:
            if id==None:
                table = CodonTable.unambiguous_dna_by_name[table]
            else:
                table = CodonTable.unambiguous_dna_by_id[id]
        elif sequence.alphabet==IUPAC.ambiguous_dna:
            if id==None:
                table = CodonTable.ambiguous_dna_by_name[table]
            else:
                table = CodonTable.ambiguous_dna_by_id[id]
        elif sequence.alphabet==IUPAC.unambiguous_rna:
            if id==None:
                table = CodonTable.unambiguous_rna_by_name[table]
            else:
                table = CodonTable.unambiguous_rna_by_id[id]
        elif sequence.alphabet==IUPAC.ambiguous_rna:
            if id==None:
                table = CodonTable.ambiguous_rna_by_name[table]
            else:
                table = CodonTable.ambiguous_rna_by_id[id]
        else:
            if id==None:
                table = CodonTable.ambiguous_generic_by_name[table]
            else:
                table = CodonTable.ambiguous_generic_by_id[id]
        sequence = sequence.tostring().upper()
        n = len(sequence)
        get = table.forward_table.get
        protein = [get(sequence[i:i+3], stop_symbol) for i in xrange(0,n-n%3,3)]
        protein = "".join(protein)
        if stop_symbol in protein :
            alphabet = Alphabet.HasStopCodon(table.protein_alphabet,
                                             stop_symbol = stop_symbol)
        else :
            alphabet = table.protein_alphabet
        return Seq(protein, alphabet)
    else:
        if id==None:
            table = CodonTable.ambiguous_generic_by_name[table]
        else:
            table = CodonTable.ambiguous_generic_by_id[id]
        get = table.forward_table.get
        sequence = sequence.upper()
        n = len(sequence)
        protein = [get(sequence[i:i+3], stop_symbol) for i in xrange(0,n-n%3,3)]
        protein = "".join(protein)
        return protein


def reverse_complement(sequence):
    """Returns the reverse complement sequence of a nucleotide string.

    If given a string, returns a new string object.
    Given a Seq or a MutableSeq, returns a new Seq object with the same alphabet.

    Supports unambiguous nucleotide sequences
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

if __name__ == "__main__" :
    print "Quick self test"
    from Bio.Data.IUPACData import ambiguous_dna_values, ambiguous_rna_values#
    from Bio.Alphabet import generic_dna, generic_rna

    #TODO - Remove this work around once we drop python 2.3 support
    try:
       set = set
    except NameError:
       from sets import Set as set

    print ambiguous_dna_complement
    for ambig_char, values in ambiguous_dna_values.iteritems() :
        compl_values = reverse_complement(values)[::-1]
        print "%s={%s} --> {%s}=%s" % \
            (ambig_char, values, compl_values, ambiguous_dna_complement[ambig_char])
        assert set(compl_values) == set(ambiguous_dna_values[ambiguous_dna_complement[ambig_char]])

    for s in ["".join(ambiguous_dna_values),
              Seq("".join(ambiguous_dna_values)),
              Seq("".join(ambiguous_dna_values), generic_dna),
              "".join(ambiguous_rna_values),
              Seq("".join(ambiguous_rna_values)),
              Seq("".join(ambiguous_dna_values), generic_rna)]:
        print "%s -> %s [RC]" % (repr(s), repr(reverse_complement(s)))
        print "%s -> %s [RNA]" % (repr(s), repr(transcribe(s)))
        print "%s -> %s [DNA]" % (repr(s), repr(back_transcribe(s)))

    #Quick check of the count method
    for letter in "ABCDEFGHIjklmnopqrstuvwyz" :
        assert 1 == Seq(letter).count(letter)
        assert 1 == MutableSeq(letter).count(letter)
    my_str = "AAAATGACGGCGGCGGCT"
    for my_obj in [my_str, Seq(my_str), MutableSeq(my_str)] :
        assert 5 == my_obj.count("A")
        assert 1 == my_obj.count("ATG")
        assert 3 == my_obj.count("CG")
        assert 2 == my_obj.count("A", 3, -5)
    for my_obj in [Seq(my_str), MutableSeq(my_str)] :
        assert 1 == my_obj.count(Seq("AT"))
        assert 5 == my_obj.count(Seq("A"))
        assert 3 == my_obj.count(Seq("CG"))
        assert 2 == my_obj.count(Seq("A"), 3, -5)
        assert 1 == my_obj.count(MutableSeq("AT"))
        assert 5 == my_obj.count(MutableSeq("A"))
        assert 3 == my_obj.count(MutableSeq("CG"))
        assert 2 == my_obj.count(MutableSeq("A"), 3, -5)
        for start in range(-len(my_str), len(my_str)) :
            for end in range(-len(my_str), len(my_str)) :
                c = my_str.count("A",start,end)
                assert c == my_str[start:end].count("A")
                assert c == my_obj.count("A",start,end)
                assert c == my_obj[start:end].count("A")
                #This one is a bit silly:
                assert my_str[start:end:-1].count("A") == my_obj[start:end:-1].count("A")

    print repr(translate(Seq("GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGT",
                             IUPAC.unambiguous_dna)))
    print repr(translate(Seq("GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG",
                             IUPAC.unambiguous_dna)))
    print repr(translate(Seq("GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG",
                             IUPAC.unambiguous_dna), stop_symbol="@"))
    
    ambig = set(IUPAC.IUPACAmbiguousDNA.letters)
    for c1 in ambig :
        for c2 in ambig :
            for c3 in ambig :
                values = set([translate(a+b+c, table=1) \
                              for a in ambiguous_dna_values[c1] \
                              for b in ambiguous_dna_values[c2] \
                              for c in ambiguous_dna_values[c3]])
                try :
                    t = translate(c1+c2+c3)
                except CodonTable.TranslationError :
                    assert "*" in values
                    continue
                if t=="*" :
                    assert values == set("*")
                elif t=="X" :
                    assert len(values) > 1, \
                        "translate('%s') = '%s' not '%s'" \
                        % (c1+c2+c3, t, ",".join(values))
                elif t=="Z" :
                    assert values == set("EQ")
                elif t=="J" :
                    assert values == set("IL")
                elif t=="B" :
                    assert values == set("DN")
                else :
                    assert values == set(t)

    print "Checking addition"
    p = Seq("PKLPAK", Alphabet.generic_protein)
    q = Seq("PKL-PAK", Alphabet.Gapped(Alphabet.generic_protein,"-"))
    r = Seq("PKL-PAK*", Alphabet.Gapped(Alphabet.HasStopCodon(Alphabet.generic_protein,"*"),"-"))
    s = Seq("PKL-PAK*", Alphabet.HasStopCodon(Alphabet.Gapped(Alphabet.generic_protein,"-"),"*"))
    t = Seq("PKLPAK*", Alphabet.HasStopCodon(Alphabet.generic_protein,"*"))
    u = "PAPKXALOA"
    v = Seq("PKLPAK!", Alphabet.HasStopCodon(Alphabet.generic_protein,"!"))
    w = Seq("PKL.PAK", Alphabet.Gapped(Alphabet.generic_protein,"."))
    for a in [p,q,r,s,t,u,v,w] :
        for b in [p,q,r,s,t,u,v,w] :
            try :
                c = a+b
            except (TypeError,ValueError), e :
                print "%s and %s -> %s" % (a.alphabet, b.alphabet, str(e))
        try :
            c = a+Seq("ACTG", Alphabet.generic_dna)
            assert isinstance(a,str), "Should have failed"
        except TypeError, e :
            pass

    print "Checking split, strip, count and find"
    for s in [p,q,r,s,t,v,w] :
        assert [x.tostring() for x in s.split()] == s.tostring().split()
        assert [x.tostring() for x in s.split("-")] == s.tostring().split("-")
        for sep, max_split in [(None,-1),(None,1),("-",-1),("L-",-1),
                               (Seq("L",Alphabet.generic_protein),2)] :
            assert [x.tostring() for x in s.split(sep,max_split)] \
                   == s.tostring().split(str(sep),max_split)
            
        assert s.strip().tostring() == s.tostring().strip()
        assert s.strip("*").tostring() == s.tostring().strip("*")
        assert s.rstrip("*").tostring() == s.tostring().rstrip("*")
        assert s.lstrip("PK").tostring() == s.tostring().lstrip("PK")
        assert s.lstrip(Seq("PK",Alphabet.generic_protein)).tostring() \
               == s.tostring().lstrip("PK")

        for text in ["-","-P", "P"] :
            assert s.find(text) == s.tostring().find(text)
            assert s.find(text,1,-2) == s.tostring().find(text,1,-2)
            assert s.find(text,2,-2) == s.tostring().find(text,2,-2)
            assert s.count(text) == s.tostring().count(text)
            assert s.count(text,1,-2) == s.tostring().count(text,1,-2)
            assert s.count(text,2,-2) == s.tostring().count(text,2,-2)
    print

    print "Checking comparisons"
    for a in [Alphabet.generic_protein, Alphabet.HasStopCodon(Alphabet.Gapped(Alphabet.generic_protein,"-"),"*")] :
        assert MutableSeq("A",a) == MutableSeq("A",a)
        assert MutableSeq("A",a) == Seq("A",a)
        assert MutableSeq("A",a) == Seq("A")
        assert MutableSeq("A",a) == "A"
        assert MutableSeq("ABC",a) == MutableSeq("A") \
               + MutableSeq("B", Alphabet.generic_protein) + MutableSeq("C",a)
        assert MutableSeq("A",a) <> MutableSeq("B")
        assert MutableSeq("ABC",a) <= MutableSeq("ABD")
        assert MutableSeq("ABC",a) < MutableSeq("ABD")
        assert MutableSeq("ABC",a) < Seq("ABD")
        assert Seq("ABC",a) < MutableSeq("ABD")
        try :
            assert MutableSeq("A",a) == MutableSeq("A",Alphabet.generic_dna)
            assert False
        except TypeError :
            pass
        """
        #TODO - Support __cmp__ for Seq object
        assert Seq("A",a) == Seq("A",a)
        assert Seq("A",a) == Seq("A")
        assert Seq("ABC",a) == Seq("A") + Seq("B", Alphabet.generic_protein) + Seq("C",a)
        assert Seq("A",a) <> Seq("B")
        assert Seq("ABC",a) <= Seq("ABD")
        assert Seq("ABC",a) < Seq("ABD")
        try :
            assert Seq("A",a) == Seq("A",Alphabet.generic_dna)
            assert False
        except TypeError :
            pass
        """
    print "Done"
