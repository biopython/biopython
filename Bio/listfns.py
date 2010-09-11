# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides useful general functions for working with lists (DEPRECATED).

This module is considered to be deprecated, and is likely to be removed in a
future release of Biopython.  Its C code implementation has already been
removed. Please get in touch via the mailing list if this will affect you.

Many of these functions can be avoided using the python set object.

Functions:
asdict        Make the list into a dictionary (for fast testing of membership).
items         Get one of each item in a list.
count         Count the number of times each item appears.
contents      Calculate percentage each item appears in a list.
itemindex     Make an index of the items in the list.
intersection  Get the items in common between 2 lists.
difference    Get the items in 1 list, but not the other.
indexesof     Get a list of the indexes of some items in a list.
take          Take some items from a list.

"""
import warnings
import Bio
warnings.warn("Bio.listfns and its C code equivalent Bio.clistfns are" \
              +" deprecated, and will be removed in a future release of"\
              +" Biopython.  If you want to continue to use this code,"\
              +" please get in contact with the Biopython developers via"\
              +" the mailing lists to avoid its permanent removal from"\
              +" Biopython. See also the Python built in set datatype.", \
              Bio.BiopythonDeprecationWarning)

def asdict(l):
    """asdict(l) -> dictionary

    Return a dictionary where the keys are the items in the list, with
    arbitrary values.  This is useful for quick testing of membership.

    """
    return count(l)

def items(l):
    """items(l) -> list of items

    Generate a list of one of each item in l.  The items are returned
    in arbitrary order.

    """
    try:
        return asdict(l).keys()
    except TypeError, x:
        if str(x).find("unhashable") == -1:
            raise
    # asdict failed because l is unhashable.  Back up to a naive
    # implementation.
    l = l[:]
    l.sort()
    i = 0
    while i < len(l)-1:
        if l[i] == l[i+1]:
            del l[i]
        else:
            i += 1
    return l

def count(items):
    """count(items) -> dict of counts of each item

    Count the number of times each item appears in a list of data.

    """
    c = {}
    for i in items:
        c[i] = c.get(i, 0) + 1
    return c

def contents(items):
    """contents(items) -> dict of item:percentage

    Summarize the contents of the list in terms of the percentages of each
    item.  For example, if an item appears 3 times in a list with 10 items,
    it is in 0.3 of the list.

    """
    counts = count(items)
    l = float(len(items))
    contents = {}
    for i, c in counts.iteritems():
        contents[i] = c / l
    return contents

def intersection(l1, l2):
    """intersection(l1, l2) -> list of common items

    Return a list of the items in both l1 and l2.  The list is in
    arbitrary order.

    """
    inter = []
    words1 = count(l1)
    for w in l2:
        if words1.has_key(w):
            inter.append(w)
            del words1[w]  # don't add the same word twice
    return inter

def difference(l1, l2):
    """difference(l1, l2) -> list of items in l1, but not l2
    
    Return a list of the items in l1, but not l2.  The list is in
    arbitrary order.

    """
    diff = []
    words2 = count(l2)
    for w in l1:
        if not words2.has_key(w):
            diff.append(w)
            words2[w] = 1   # don't add the same word twice
    return diff

def itemindex(l):
    """itemindex(l) -> dict of item : index of item

    Make an index of the items in the list.  The dictionary contains
    the items in the list as the keys, and the index of the first
    occurrence of the item as the value.

    """
    dict = {}
    for i in range(len(l)):
        if not dict.has_key(l[i]):
            dict[l[i]] = i
    return dict

def indexesof(l, fn, opposite=0):
    """indexesof(l, fn) -> list of indexes

    Return a list of indexes i where fn(l[i]) is true.

    """
    indexes = []
    for i in range(len(l)):
        f = fn(l[i])
        if (not opposite and f) or (opposite and not f):
            indexes.append(i)
    return indexes

def take(l, indexes):
    """take(l, indexes) -> list of just the indexes from l"""
    items = []
    for i in indexes:
        items.append(l[i])
    return items

def take_byfn(l, fn, opposite=0):
    indexes = indexesof(l, fn, opposite=opposite)
    return take(l, indexes)
