#!/usr/bin/env python

__version__ = "$Revision: 1.6 $"

from __future__ import division

import commands
import itertools
import os
import re

from Bio import Wise

_SCORE_MATCH = "4"
_SCORE_MISMATCH = "-1"
_PENALTY_GAP_START = "5"
_PENALTY_GAP_EXTENSION = "1"

_CMDLINE_DNAL = ["dnal", "-alb", "-nopretty",
                "-match", _SCORE_MATCH, "-mis", _SCORE_MISMATCH,
                "-gap", _PENALTY_GAP_START, "-ext", _PENALTY_GAP_EXTENSION]

_CMDLINE_FGREP_COUNT = "fgrep -c '%s' %s"
def _fgrep_count(pattern, file):
    return int(commands.getoutput(_CMDLINE_FGREP_COUNT % (pattern, file)))

_re_alb_line2coords = re.compile(r"^\[([^:]+):[^\[]+\[([^:]+):")
def _alb_line2coords(line):
    return tuple([int(coord)+1 # one-based -> zero-based
                  for coord
                  in _re_alb_line2coords.match(line).groups()])

def _get_coords(filename):
    alb = file(filename)

    start_line = None
    end_line = None

    for line in alb:
        if line.startswith("["):
            if not start_line:
                start_line = line # rstrip not needed
            else:
                end_line = line

    return zip(*map(_alb_line2coords, [start_line, end_line])) # returns [(start0, end0), (start1, end1)]

def _any(seq, pred=bool):
    "Returns True if pred(x) is True at least one element in the iterable"
    return True in itertools.imap(pred, seq)

class Statistics(object):
    """
    Calculate statistics from an ALB report
    """
    def __init__(self, filename):
        self.matches = _fgrep_count('"SEQUENCE" %s' % _SCORE_MATCH, filename)
        self.mismatches = _fgrep_count('"SEQUENCE" %s' % _SCORE_MISMATCH, filename)
        self.gaps = _fgrep_count('"INSERT" -%s' % _PENALTY_GAP_START, filename)
        self.extensions = _fgrep_count('"INSERT" -%s' % _PENALTY_GAP_EXTENSION, filename)

        if _any([self.matches, self.mismatches, self.gaps, self.extensions]):
            self.coords = _get_coords(filename)
        else:
            self.coords = [(0, 0), (0,0)]

    def identity_fraction(self):
        return self.matches/(self.matches+self.mismatches)

    header = "identity_fraction\tmatches\tmismatches\tgaps\textensions"

    def __str__(self):
        return "\t".join([str(x) for x in (self.identity_fraction(), self.matches, self.mismatches, self.gaps, self.extensions)])

def align(pair, *args, **keywds):
    temp_file = Wise.align(_CMDLINE_DNAL, pair, *args, **keywds)
    try:
        return Statistics(temp_file.name)
    except AttributeError:
        try:
            keywds['dry_run']
            return None
        except KeyError:
            raise

def main():
    import sys
    stats = align(sys.argv[1:3])
    print "\n".join(["%s: %s" % (attr, getattr(stats, attr))
                     for attr in
                     ("matches", "mismatches", "gaps", "extensions")])
    print "identity_fraction: %s" % stats.identity_fraction()
    print "coords: %s" % stats.coords

def _test(*args, **keywds):
    import doctest, sys
    doctest.testmod(sys.modules[__name__], *args, **keywds)

if __name__ == "__main__":
    if __debug__:
        _test()
    main()
