# Copyright 2000-2002 Brad Chapman.
# Copyright 2004-2005 by M de Hoon.
# Copyright 2007 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import string, array

import Alphabet
from Alphabet import IUPAC
from Data.IUPACData import ambiguous_dna_complement, ambiguous_rna_complement
from Bio.Data import CodonTable

class Seq:
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

    # I don't think I like this method...
##    def __cmp__(self, other):
##        if isinstance(other, Seq):
##            return cmp(self.data, other.data)
##        else:
##            return cmp(self.data, other)

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
        if type(other) == type(' '):
            return self.__class__(self.data + other, self.alphabet)
        elif self.alphabet.contains(other.alphabet):
            return self.__class__(self.data + other.data, self.alphabet)
        elif other.alphabet.contains(self.alphabet):
            return self.__class__(self.data + other.data, other.alphabet)
        else:
            raise TypeError, ("incompatable alphabets", str(self.alphabet),
                              str(other.alphabet))
        
    def __radd__(self, other):
        if self.alphabet.contains(other.alphabet):
            return self.__class__(other.data + self.data, self.alphabet)
        elif other.alphabet.contains(self.alphabet):
            return self.__class__(other.data + self.data, other.alphabet)
        else:
            raise TypeError, ("incompatable alphabets", str(self.alphabet),
                              str(other.alphabet))


    def tostring(self):                            # Seq API requirement
        """Returns the full sequence as a python string.

        Although not formally deprecated, you are now encouraged to use
        str(my_seq) instead of my_seq.tostring()."""
        return self.data

    def tomutable(self):   # Needed?  Or use a function?
        return MutableSeq(self.data, self.alphabet)
    
    def count(self, sub, start=None, end=None):
        """Count method, like that of a python string.

        Return an integer, the number of occurrences of substring
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
        try :
            #Assume sub is a Seq or MutableSeq object, so pass
            #it to the string's count method as a string.
            #TODO - Should we check the alphabet?
            search = sub.tostring()
        except AttributeError:
            #Assume sub is a string.
            search = sub

        #TODO - More elegant way of doing this splice notation
        if start is None and end is None :
            return self.data.count(search)
        elif end is None :
            return self.data.count(search, start)
        else :
            return self.data.count(search, start, end)

    def __maketrans(self, alphabet) :
        """Seq.__maketrans(alphabet) -> translation table.

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
        which need to be backwards compatible with old Biopython, you
        should continue to use my_seq.tostring() rather than str(my_seq).
        """
        return "".join(self.data)

    def __cmp__(self, other):
        if isinstance(other, MutableSeq):
            x = cmp(self.alphabet, other.alphabet)
            if x == 0:
                return cmp(self.data, other.data)
            return x
        elif type(other) == type(""):
            return cmp(self.data.tostring(), other)
        elif isinstance(other, Seq):
            x = cmp(self.alphabet, other.alphabet)
            if x == 0:
                return cmp(self.data.tostring(), other.data)
            return x
        else:
            return cmp(self.data, other)

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
        if self.alphabet.contains(other.alphabet):
            return self.__class__(self.data + other.data, self.alphabet)
        elif other.alphabet.contains(self.alphabet):
            return self.__class__(self.data + other.data, other.alphabet)
        else:
            raise TypeError, ("incompatable alphabets", str(self.alphabet),
                              str(other.alphabet))
    def __radd__(self, other):
        if self.alphabet.contains(other.alphabet):
            return self.__class__(other.data + self.data, self.alphabet)
        elif other.alphabet.contains(self.alphabet):
            return self.__class__(other.data + self.data, other.alphabet)
        else:
            raise TypeError, ("incompatable alphabets", str(self.alphabet),
                              str(other.alphabet))

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

    def count(self, sub, start=None, end=None):
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
        if len(sub) == 1 :
            #Try and be efficient and work directly from the array.
            try :
                #Assume item is a single letter Seq/MutableSeq
                #TODO - Should we check the alphabet?
                letter = sub.tostring()
            except AttributeError :
                letter = sub

            count = 0
            #TODO - More elegant way of doing this splice notation
            if start is None and end is None :
                for c in self.data:
                    if c == letter: count += 1
            elif end is None :
                for c in self.data[start:]:
                    if c == letter: count += 1
            else :
                for c in self.data[start:end]:
                    if c == letter: count += 1
            return count
        else :
            #TODO - Can we do this more efficiently?
            #We could use the self.tostring().count(...) method
            return self.toseq().count(sub, start, end)

    def index(self, item):
        for i in range(len(self.data)):
            if self.data[i] == item:
                return i
        raise ValueError, "MutableSeq.index(x): x not in list"

    def reverse(self):
        """Modify the MutableSequence to reverse itself.

        No return value"""
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

        No return value"""
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
        str(my_seq) instead of my_seq.tostring()."""
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
        #TODO - Raise an error if already is RNA alphabet?
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
        #TODO - Raise an error if already is DNA alphabet?
        dna = rna.data.replace('U','T').replace('u','t')
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
    from sets import Set
    print ambiguous_dna_complement
    for ambig_char, values in ambiguous_dna_values.iteritems() :
        compl_values = reverse_complement(values)[::-1]
        print "%s={%s} --> {%s}=%s" % \
            (ambig_char, values, compl_values, ambiguous_dna_complement[ambig_char])
        assert Set(compl_values) == Set(ambiguous_dna_values[ambiguous_dna_complement[ambig_char]])

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
    
    ambig = Set(IUPAC.IUPACAmbiguousDNA.letters)
    for c1 in ambig :
        for c2 in ambig :
            for c3 in ambig :
                values = Set([translate(a+b+c, table=1) \
                              for a in ambiguous_dna_values[c1] \
                              for b in ambiguous_dna_values[c2] \
                              for c in ambiguous_dna_values[c3]])
                try :
                    t = translate(c1+c2+c3)
                except CodonTable.TranslationError :
                    assert "*" in values
                    continue
                if t=="*" :
                    assert values == Set(["*"])
                elif t=="X" :
                    assert len(values) > 1, \
                        "translate('%s') = '%s' not '%s'" \
                        % (c1+c2+c3, t, ",".join(values))
                elif t=="Z" :
                    assert values == Set(("E", "Q"))
                elif t=="B" :
                    assert values == Set(["D", "N"])
                else :
                    assert values == Set(t)
 
