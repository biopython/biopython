# Copyright 2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This provides useful general functions for working with lists.

Functions:
asdict        Make the list into a dictionary (for fast testing of membership).
items         Get one of each item in a list.
count         Count the number of times each item appears.
contents      Calculate percentage each item appears in a list.
itemindex     Make an index of the items in the list.
intersection  Get the items in common between 2 lists.
difference    Get the items in 1 list, but not the other.

"""

from clistfns import *

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
    return asdict(l).keys()

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
