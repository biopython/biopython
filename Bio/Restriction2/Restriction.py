# Copyright 2013 Antony Lee
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Restriction enzyme and sequence objects, and analysis utilities.
"""

from __future__ import division, print_function

from collections import namedtuple
import re

from . import RestrictionDB


_codes = {"A": "A", "C": "C", "G": "G", "T": "T",
          "R": "GA", "Y": "TC", "K": "GT", "M": "AC", "S": "GC", "W": "AT",
          "B": "GTC", "D": "GAT", "H": "ACT", "V": "GCA", "N": "ACGT"}


_complement = {"A": "T", "C": "G", "G": "C", "T": "A",
               "R": "Y", "Y": "R", "M": "K", "K": "M", "S": "S", "W": "W",
               "B": "D", "D": "B", "H": "V", "V": "H", "N": "N"}


class Seq(str):
    """A DNA sequence object, including linear/circular information.
    """

    def __new__(cls, s, circular=False):
        self = super(Seq, cls).__new__(Seq, s)
        self.circular = circular
        return self

    @property
    def reverse_complement(self):
        return Seq("".join(_complement[base] for base in reversed(self)),
                   self.circular)


class Restriction(namedtuple("Restriction",
                             ("name", "site", "cuts", "methylation", "suppliers"))):
    """Type of all enzymes.

    The following attributes are defined:
        - name: enzyme name
        - site: recongnized site
        - cuts: list of (5', 3') cutting sites
        - methylation: list of methylation site and type pairs
        - suppliers: string of REBASE codes

    The following class attributes are defined:
        - all: all defined enzymes
        - commercial: all commercially available enzymes

    Enzymes with unknown cutting sites are ignored.
    """

    all = {}
    commercial = {}

    @classmethod
    def load_all(cls):
        """Populate cls.all and cls.commercial from the database.
        """
        patterns = RestrictionDB.patterns()
        information = RestrictionDB.information()
        for name, pattern in patterns.items():
            if not pattern.cuts:
                continue # ignore enzymes with unknown cutting site
            if len(pattern.cuts) == 2:
                (c1, c2), (c3, c4) = pattern.cuts
                if c2 - c1 != c4 - c3:
                    continue # ignore HaeIV (EMBOSS v301)
            info = information[name]
            cls.all[name] = cls(name, pattern.site, pattern.cuts,
                                info.methylation, info.suppliers)
            if info.suppliers:
                cls.commercial[name] = cls.all[name]

    def __len__(self):
        """Return the length of the recognition site.
        """
        return len(self.site)

    def isoschizomer(self, other):
        """Return whether two enzymes have the same site.
        """
        return isinstance(other, Restriction) and self.site == other.site

    def neoschizomer(self, other):
        """Return whether two enzymes have the same site but cut differently.
        """
        return self.isoschizomer(other) and self.cuts != other.cuts

    def isoschizomers(self, batch=None):
        """Lists the isoschizomers of an enzyme from a batch, or all enzymes.
        """
        if batch is None:
            batch = self.all
        return sorted(x for x in batch if x.isoschizomer(self) and x != self)

    def neoschizomers(self, batch=None):
        """Lists the neoschizomers of an enzyme from a batch, or all enzymes.
        """
        if batch is None:
            batch = self.all
        return sorted(x for x in batch if x.neoschizomer(self))

    def __str__(self):
        """Return the enzyme name followed by its pattern.
        """
        cuts_pos, cuts_neg = zip(*self.cuts)
        markers = sorted([(cut, "^") for cut in cuts_pos] +
                         [(cut, "_") for cut in cuts_neg])
        offset5 = max(-markers[0][0], 0)
        offset3 = max(markers[-1][0] - len(self.site), 0)
        site = "N" * offset5 + self.site + "N" * offset3
        str = ""
        previdx = 0
        for idx, char in markers:
            idx += offset5
            str += site[previdx:idx] + char
            previdx = idx
        str += site[previdx:]
        return self.name + ":" + str

    @property
    def overhang(self):
        """Return an overhang type and sequence pair.

        The overhang type is one of 5, 3 or None and the overhang sequence is
        given on the + strand, 5'-to-3'.
        """
        str_self = str(self)
        cut_pos = str_self.find("^")
        cut_neg = str_self.find("_")
        if cut_neg == cut_pos + 1:
            return None, ""
        if cut_pos > cut_neg:
            return 3, str_self[cut_neg+1:cut_pos]
        else:
            return 5, str_self[cut_pos+1:cut_neg]

    def always_compatible(self, other):
        """Test whether the product of two digests are always compatible.
        """
        type1, ov1 = self.overhang
        type2, ov2 = other.overhang
        return (type1 == type2 and ov1 == ov2 and
                "N" not in ov1 and "N" not in ov2)

    def perhaps_compatible(self, other):
        """Test whether the product of two digests may be compatible.
        """
        type1, ov1 = self.overhang
        type2, ov2 = other.overhang
        return (type1 == type2 and len(ov1) == len(ov2) and
                all(set(_codes[b1]).intersection(_codes[b2])
                    for b1, b2 in zip(ov1, ov2)))

    @property
    def frequency(self):
        """Return the inverse of the frequency of a site.
        """
        f = 1
        for b in self.site:
            f *= 4 / len(_codes[b])
        return f

    def search(self, seq):
        """List in order the indices after which an enzyme cuts on a + strand.
        """
        if seq.circular:
            seq += seq[:len(self.site) - 1]
        regexp_pos = "".join("[" + _codes[b] + "]" for b in self.site)
        regexp_neg = "".join("[" + _codes[b] + "]"
                             for b in Seq(self.site).reverse_complement)
        matches = [match.start() + cut_pos
                   for match in re.finditer(regexp_pos, seq)
                   for cut_pos, cut_neg in self.cuts]
        if regexp_pos != regexp_neg:
            matches += [match.start() + len(self.site) - cut_neg
                        for match in re.finditer(regexp_neg, seq)
                        for cut_pos, cut_neg in self.cuts]
        return sorted(matches)

    def catalyze(self, seq):
        """List the + strands obtained by a digest.
        """
        fragments = []
        previdx = 0
        for idx in self.search(seq):
            fragments.append(seq[previdx:idx])
            previdx = idx
        fragments.append(seq[previdx:])
        if seq.circular:
            fragments[0] = fragments[-1] + fragments[0]
            fragments.pop()
        return fragments


Restriction.load_all()


def analyze(seq, enzymes=Restriction.all.values()):
    """Maps enzymes (by default, all) to the points they cut in a sequence.
    """
    return {enzyme.name: enzyme.search(seq) for enzyme in enzymes}


def pretty_str(seq, analysis):
    """Creates a prettified string of a restriction analysis.
    """
    width = 80
    start = 0
    pretty = ""
    while seq:
        stop = start + width
        cps = {cp: [name for name, cps in analysis.items() if cp in cps]
               for cp in {cp for cps in analysis.values() for cp in cps
                          if start <= cp < stop}}
        ss = [[" "] * width]
        markers = [set()]
        for cp in sorted(cps):
            for name in cps[cp]:
                for s, marker in zip(ss, markers):
                    if all(c == " " for c in s[cp:cp+len(name)]):
                        s[cp:cp+len(name)] = name
                        marker.add(cp)
                        break
                else:
                    ss.append([" "] * len(seq))
                    ss[-1][cp:cp+len(name)] = name
                    markers.append({cp})
        ss.insert(0, [" "] * len(seq))
        for i, s in enumerate(ss):
            for j, c in enumerate(s):
                if c == " " and any(j in marker for marker in markers[i:]):
                    s[j] = u"\u258f"
        pretty += ("\n".join("".join(s) for s in reversed(ss)) + "\n"
                   + seq + "\n" + Seq(seq).reverse_complement[::-1] + "\n")
        start = stop
        seq = seq[width:]
    return pretty
