#!/usr/bin/env python

__version__ = "$Revision: 1.2 $"

from __future__ import division

import commands
import os

from Bio import Wise

_SCORE_MATCH = "4"
_SCORE_MISMATCH = "-1"
_PENALTY_GAP_START = "5"
_PENALTY_GAP_EXTENSION = "1"

_CMDLINE_DNAL = ["dnal", "-alb", "-nopretty",
                "-match", _SCORE_MATCH, "-mis", _SCORE_MISMATCH,
                "-gap", _PENALTY_GAP_START, "-ext", _PENALTY_GAP_EXTENSION]

try:
    import lsf
    if lsf.jobid:
        _CMDLINE_DNAL.extend(["-silent", "-quiet"])
except ImportError:
    pass

_CMDLINE_FGREP_COUNT = "fgrep -c '%s' %s"
def _fgrep_count(pattern, file):
    return int(commands.getoutput(_CMDLINE_FGREP_COUNT % (pattern, file)))

class Statistics(object):
    """
    Calculate statistics from an ALB report
    """
    def __init__(self, filename):
        self.matches = _fgrep_count('"SEQUENCE" %s' % _SCORE_MATCH, filename)
        self.mismatches = _fgrep_count('"SEQUENCE" %s' % _SCORE_MISMATCH, filename)
        self.gaps = _fgrep_count('"INSERT" -%s' % _PENALTY_GAP_START, filename)
        self.extensions = _fgrep_count('"INSERT" -%s' % _PENALTY_GAP_EXTENSION, filename)

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

def _test(*args, **keywds):
    import doctest, sys
    doctest.testmod(sys.modules[__name__], *args, **keywds)

if __name__ == "__main__":
    if __debug__:
        _test()
    main()
