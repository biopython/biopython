# Copyright 2013 Antony Lee
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Utilities for downloading and parsing REBASE data files in EMBOSS format.

Importing this module will have the side-effect of trying to fetch the REBASE
data files if they are not present!
"""

from collections import namedtuple
import os
import urllib2


_dir = os.path.dirname(__file__)


def suppliers():
    """Return a dict of supplier keys to supplier names.
    """
    with open(os.path.join(_dir, "emboss_s.txt")) as f:
        return dict(line.strip().split(None, 1)
                    for line in f if not line.startswith("#"))


def patterns():
    """Return a dict of enzyme names to enzyme patterns.
    """
    Pattern = namedtuple("Pattern", ("name", "site", "cuts"))
    def to_restriction(name, site, length, ncuts, blunt, c1, c2, c3, c4):
        c1, c2, c3, c4 = map(int, (c1, c2, c3, c4))
        if ncuts == "0":
            cuts = []
        elif ncuts == "2":
            cuts = [(c1, c2)]
        elif ncuts == "4":
            cuts = [(c1, c2), (c3, c4)]
        return name, Pattern(name, site.upper(), cuts)
    with open(os.path.join(_dir, "emboss_e.txt")) as f:
        return dict(to_restriction(*line.strip().split())
                    for line in f if not line.startswith("#"))


def information():
    """Return a dict of enzyme names to enzyme informations.
    """
    Information = namedtuple(
        "Information",
        ("name", "organism", "isoschizomers", "methylation", "source",
         "suppliers", "references"))
    informations = {}
    with open(os.path.join(_dir, "emboss_r.txt")) as f:
        for enz in "".join(line for line in f
                           if not line.startswith("#")).split("//\n")[:-1]:
            (name, organism, isoschizomers, _methylation, source, suppliers,
             nref, references) = enz.split("\n", 7)
            references = references.splitlines()
            methylation = []
            if _methylation:
                for spec in _methylation.split(","):
                    base = spec[:spec.index("(")]
                    base = None if base == "?" else int(base)
                    type = int(spec[-2])
                    methylation.append((base, type))
            informations[name] = Information(
                name, organism, isoschizomers, methylation, source, suppliers,
                references)
    return informations


def update_db():
    """Download the EMBOSS files from the NEB server and check their validity.
    """
    for fname in ["emboss_e.txt", "emboss_r.txt", "emboss_s.txt"]:
        file = urllib2.urlopen("ftp://ftp.neb.com/pub/rebase/" + fname)
        with open(os.path.join(_dir, fname), "w") as out:
            out.write(file.read())
        file.close()
    e_names = []
    r_names = []
    with open(os.path.join(_dir, "emboss_e.txt")) as f:
        for line in f:
            if line.startswith("#"):
                continue
            name, site, length, ncuts, blunt, c1, c2, c3, c4 = line.strip().split()
            e_names.append(name)
            assert len(site) == int(length)
            c1, c2, c3, c4 = map(int, (c1, c2, c3, c4))
            if ncuts == "0":
                assert c1 == c2 == c3 == c4 == 0
            elif ncuts == "2":
                assert c3 == c4 == 0
            elif ncuts == "4":
                pass
            else:
                raise AssertionError
            if blunt == "1":
                assert ncuts == "2" and c1 == c2
    with open(os.path.join(_dir, "emboss_r.txt")) as f:
        for enz in "".join(line for line in f
                           if not line.startswith("#")).split("//\n")[:-1]:
            (name, organism, isoschizomers, methylation, source, suppliers,
             nref, references) = enz.split("\n", 7)
            r_names.append(name)
            references = references.splitlines()
            assert len(references) == int(nref)
    assert sorted(e_names) == sorted(r_names)


def ensure_db():
    """Download the missing EMBOSS files, if needed.
    """
    required = ["emboss_{}.txt".format(c) for c in "ers"]
    if any(fname not in os.listdir(_dir) for fname in required):
        update_db()


ensure_db()
