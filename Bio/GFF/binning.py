#!/usr/bin/env python
#
# Copyright 2002 by Michael Hoffman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# note:
# I used Lincoln Stein's perl Bio::DB::GFF::Util::Binning as a model for this

"""Binning support for Bio.GFF (DEPRECATED)

This is part of the "old" Bio.GFF module by Michael Hoffman, which offered
access to a MySQL database holding GFF data loaded by BioPerl. This code has
now been deprecated, and will probably be removed in order to free the Bio.GFF
namespace for a new GFF parser in Biopython (including GFF3 support).
"""

class Meta(object):
    MIN_BIN = 1000
    MAX_BIN = 100000000
    
def query(start=0, stop=Meta.MAX_BIN, minbin=Meta.MIN_BIN, maxbin=Meta.MAX_BIN):
    args = []
    bins = []
    tier = maxbin

    if start is None:
        start=0
    if stop is None:
        stop = Meta.MAX_BIN

    while tier >= minbin:
        tier_start, tier_stop = bot(tier, start), top(tier, stop)
        if (tier_start == tier_stop):
            bins.append('fbin=%s')
            args.append(tier_start)
        else:
            bins.append('fbin BETWEEN %s AND %s')
            args.extend([tier_start, tier_stop])
        tier /= 10
    query = "\n\t OR ".join(bins)
    
    return query % tuple(args)

def bin(start, stop, min, wantarray=0):
    tier = min
    while 1:
        bin_start = int(float(start)/tier)
        bin_end = int(float(stop)/tier)
        if bin_start == bin_end:
            break
        tier *= 10
    if wantarray:
        return tier, bin_start
    else:
        return bin_name(tier, bin_start)

def bot(tier, pos):
    return name(tier, int(pos/tier))

top = bot

def name(x, y):
    return "%d.%06d" % (x, y)

def _test(*args, **keywds):
    import doctest, sys
    doctest.testmod(sys.modules[__name__], *args, **keywds)

if __name__ == "__main__":
    if __debug__:
        _test()
