import string, array

import Alphabet

class Seq:
    def __init__(self, data, alphabet = Alphabet.generic_alphabet):
        # Enforce string storage
        assert (type(data) == type("") or # must use a string
                type(data) == type(u""))  # but can be a unicode string

        self.data = data                           # Seq API requirement
        self.alphabet = alphabet                   # Seq API requirement

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__,
                               repr(self.data),
                               repr(self.alphabet))
    def __str__(self):
        if len(self.data) > 60:
            s = repr(self.data[:60] + " ...")
        else:
            s = repr(self.data)
        return "%s(%s, %s)" % (self.__class__.__name__, s,
                               repr(self.alphabet))
    # I don't think I like this method...
##    def __cmp__(self, other):
##        if isinstance(other, Seq):
##            return cmp(self.data, other.data)
##        else:
##            return cmp(self.data, other)
    def __len__(self): return len(self.data)       # Seq API requirement
    def __getitem__(self, i): return self.data[i]  # Seq API requirement
    def __getslice__(self, i, j):                  # Seq API requirement
        i = max(i, 0); j = max(j, 0)
        return Seq(self.data[i:j], self.alphabet)
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
        return self.data

    def tomutable(self):   # Needed?  Or use a function?
        return MutableSeq(self.data, self.alphabet)
    
    def count(self, item):
        return len([x for x in self.data if x == item])

class MutableSeq:
    def __init__(self, data, alphabet = Alphabet.generic_alphabet):
        if type(data) == type(""):
            self.data = array.array("c", data)
        else:
            self.data = data   # assumes the input is an array
        self.alphabet = alphabet
    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__,
                               repr(self.data),
                               repr(self.alphabet))

    def __str__(self):
        if len(self.data) > 60:
            s = repr(string.join(self.data[:60], "") + " ...")
        else:
            s = repr(string.join(self.data, ""))
        return "%s(%s, %s)" % (self.__class__.__name__, s,
                               repr(self.alphabet))
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
    def __getitem__(self, i): return self.data[i]
    def __setitem__(self, i, item): self.data[i] = item
    def __delitem__(self, i): del self.data[i]
    def __getslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        return self.__class__(self.data[i:j], self.alphabet)
    def __setslice__(self, i, j, other):
        i = max(i, 0); j = max(j, 0)
        if isinstance(other, MutableSeq):
            self.data[i:j] = other.data
        elif isinstance(other, type(self.data)):
            self.data[i:j] = other
        else:
            self.data[i:j] = array.array("c", str(other))
    def __delslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        del self.data[i:j]
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
    def count(self, item):
        count = 0
        for c in self.data:
            if c == item:
                count = count + 1
        return count
    def index(self, item):
        for i in range(len(self.data)):
            if self.data[i] == item:
                return i
        raise ValueError, "MutableSeq.index(x): x not in list"
    def reverse(self):
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
        return string.join(self.data, "")

    def toseq(self):
        return Seq(string.join(self.data, ""), self.alphabet)
