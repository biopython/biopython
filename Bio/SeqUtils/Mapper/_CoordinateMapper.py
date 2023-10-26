# AUTHOR: Reece Hart <reecehart@gmail.com>
# Modifications copyright 2012 Lenna X. Peterson <arklenna@gmail.com>

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

r"""CoordinateMapper -- map between genomic, cds, and protein coordinates.

Examples::
    AB026906.1:g.7872G>T
    AB026906.1:c.274G>T
    BA...:p.Asp92Tyr

All refer to congruent variants. A coordinate mapper is needed to
translate between at least these three coordinate frames.  The mapper
should deal with specialized syntax for splicing and UTR (e.g., `88+1`,
`89-2`, `-14`, `*46`). In addition, care should be taken to ensure consistent 0-
or 1-based numbering (0 internally, as with Python/BioPython and
Perl/BioPerl).

.. code-block::

    g   -----------00000000000-----1111111----------22222222222*222-----
                   s0        e0    s1    e1         s2            e2
                   \         \     |     |          /             /
                    +--+      +--+ |     | +-------+     +-------+
                        \         \|     |/             /
    c                   00000000000111111122222222222*222
                                  c0     c1             c2
                        aaabbbcccdddeeefffggghhhiiijj*kkk
    p                     A  B  C  D  E  F  G  H  I  J  K
                          p0 p1 p2 ...                  pn

"""

from functools import wraps
from math import floor
import typing
import warnings

from Bio import BiopythonParserWarning
from Bio.SeqFeature import FeatureLocation
from Bio.SeqUtils.Mapper import GenomePosition, CDSPosition, ProteinPosition
from Bio.SeqUtils.Mapper.MapPositions import CDSPositionError, ProteinPositionError
from Bio.SeqUtils.Mapper import MapPositions  # for decorator


def pos_factory(pos_type: str) -> typing.Callable:
    """
    Convert string or int pos to appropriate Position object.

    Parameters
    ----------
    pos_type : str
        Position type (Genome, CDS, Protein)
    """

    def wrapper(fn):
        @wraps(fn)
        def make_pos(self, pos, dialect=None):
            # retrieve Position object
            _obj = getattr(MapPositions, pos_type + "Position")
            # if pos is not already Position object, make it one
            if not isinstance(pos, _obj):
                # no dialect: use default constructor
                if dialect is None:
                    pos = _obj(pos)
                # use dialect alternate constructor
                else:
                    pos = _obj.from_dialect(dialect, pos)
            # call function with new pos
            return fn(self, pos, dialect)

        return make_pos

    return wrapper


class CoordinateMapper(object):
    """Convert positions between genomic, CDS, and protein coordinates."""

    def __init__(self, selist):
        """Set exons to be used for mapping.

        Parameters
        ----------
        selist : SeqRecord, SeqFeature, list
            Object containing exon information
        """
        self._exons = self._get_exons(selist)

    @staticmethod
    def _get_exons(seq):
        """Extract exons from SeqRecord, SeqFeature, or list of pairs.

        Parameters
        ----------
        seq : SeqRecord, SeqFeature, list
            Object containing exon information.

        Returns
        -------
        SeqFeature.FeatureLocation
        """
        # Try as SeqRecord
        if hasattr(seq, "features"):
            # generator
            cdsf = next(f for f in seq.features if f.type == "CDS")
            return cdsf.location
        # Try as SeqFeature
        elif hasattr(seq, "location"):
            if seq.type != "CDS":
                # FIXME should this be a fatal error?
                warnings.warn(
                    "Provided SeqFeature should be CDS", BiopythonParserWarning
                )
            return seq.location
        # Try as list of pairs
        return sum([FeatureLocation(s, e) for s, e in seq])

    @property  # read-only
    def exons(self):
        """Get exons."""
        return self._exons

    @property  # read-only
    def exon_list(self):
        """Get list of exons."""
        return list(self.exons)

    @pos_factory("Genome")
    def g2c(self, gpos, dialect=None):
        """Convert integer from genomic to CDS coordinates.

        Parameters
        ----------
        gpos : int
            Genomic position
        dialect : str
            Coordinate dialect (GenBank or HGVS, default None)

        Returns
        -------
        CDSPosition
        """
        gpos = int(gpos)
        fmts = CDSPosition.fmt_dict

        def _simple_g2c(g):
            return self.exon_list.index(g)

        # within exon
        if gpos in self.exons:
            return CDSPosition(_simple_g2c(gpos))
        # before CDS
        if gpos < self.exons.start:
            return CDSPosition(fmts["post-CDS"].format(offset=gpos - self.exons.start))
        # after CDS
        if gpos >= self.exons.end:
            return CDSPosition(fmts["pre-CDS"].format(offset=gpos - self.exons.end + 1))
        # intron
        # set start of first intron
        prev_end = self.exons.parts[0].end
        for part in self.exons.parts[1:]:
            # not in this intron
            if gpos > part.start:
                prev_end = part.end
                continue
            len_intron = part.start - prev_end
            # second half (exclusive) of intron
            #       offset     >     middle of intron
            if gpos - prev_end > floor(len_intron / 2.0):
                anchor = _simple_g2c(part.start)
                offset = gpos - part.start
                assert offset < 0
            # first half (inclusive) of intron
            else:
                anchor = _simple_g2c(prev_end - 1)
                offset = gpos - prev_end + 1
                assert offset > 0
            assert self.check_intron(anchor, offset)
            return CDSPosition(fmts["intron"].format(pos=anchor, offset=offset))

        raise AssertionError  # function should return for every integer

    # TODO verify that values of offset are sane
    def check_intron(self, anchor, offset):
        """
        Verify that CDS-relative intron position is valid with given exons.

        Parameters
        ----------
        anchor : int
            Intron anchor (closest CDS position)
        offset : int
            Intron offset (distance to anchor)

        Returns
        -------
        bool
        """
        for exon in self.exons.parts:
            start = int(self.g2c(exon.start))
            if anchor == start:
                if offset > 0:
                    raise CDSPositionError(
                        "Invalid intron: offset from exon start must be negative."
                    )
                return True
            end = int(self.g2c(exon.end - 1))
            if anchor == end:
                if offset < 0:
                    raise CDSPositionError(
                        "Invalid intron: offset from exon end must be positive."
                    )
                return True
        raise CDSPositionError(
            "Invalid intron: anchor must be start or end of an exon."
        )

    def get_strand(self, gpos):
        """Get strand."""
        for exon in self.exons.parts:
            if gpos in exon:
                return exon.strand
        raise ValueError("Provided gpos must be exon")

    @pos_factory("CDS")
    def c2g(self, cpos, dialect=None):
        """Convert from CDS to genomic coordinates.

        Parameters
        ----------
        cpos : int
            CDS position
        dialect : str
            Coordinate dialect (GenBank or HGVS, default None)

        Returns
        -------
        GenomePosition
        """
        if cpos.pos_type == "pre-CDS":
            return GenomePosition(self.exons.start + cpos.offset)
        elif cpos.pos_type == "post-CDS":
            return GenomePosition(self.exons.end - 1 + cpos.offset)

        g_anchor = self.exon_list[cpos.anchor]
        if cpos.pos_type == "exon":
            strand = self.get_strand(g_anchor)
            return GenomePosition(g_anchor, strand=strand)
        elif cpos.pos_type == "intron":
            offset = cpos.offset
            if self.check_intron(cpos.anchor, offset):
                return GenomePosition(g_anchor + offset)

        raise AssertionError  # all positions should be one of the 4 types

    @pos_factory("CDS")
    def c2p(self, cpos, dialect=None):
        """Convert from CDS to protein coordinates.

        Parameters
        ----------
        cpos : int
            CDS position
        dialect : str
            Coordinate dialect (GenBank or HGVS, default None)

        Returns
        -------
        ProteinPosition
        """
        try:
            cpos = int(cpos)
        except TypeError:
            raise CDSPositionError("'%s' does not correspond to a protein" % repr(cpos))
        return ProteinPosition(int(cpos / 3.0))

    @pos_factory("Genome")
    def g2p(self, gpos, dialect=None):
        """Convert integer from genomic to protein coordinates.

        Parameters
        ----------
        gpos : int
            Genomic position
        dialect : str
            Coordinate dialect (GenBank or HGVS, default None)

        Returns
        -------
        ProteinPosition
        """
        return self.c2p(self.g2c(gpos))

    @pos_factory("Protein")
    def p2c(self, ppos, dialect=None):
        """Convert integer from protein coordinate to CDS closed range.

        Parameters
        ----------
        ppos : int
            Protein position
        dialect : str
            Coordinate dialect (GenBank or HGVS, default None)

        Returns
        -------
        CDSPosition
        """
        try:
            ppos = int(ppos)
        except TypeError:
            return None
        if ppos < 0:
            raise ProteinPositionError("'%s' should not be negative")
        # FIXME is CDS guaranteed to have len % 3 == 0?
        first_base = (ppos) * 3
        last_base = first_base + 2
        if last_base > len(self.exons):
            raise ProteinPositionError("'%s' is too large")
        return (CDSPosition(first_base), CDSPosition(last_base))

    @pos_factory("Protein")
    def p2g(self, ppos, dialect=None):
        """Convert integer from protein to genomic coordinates.

        Parameters
        ----------
        ppos : int
            Protein position
        dialect : str
            Coordinate dialect (GenBank or HGVS, default None)

        Returns
        -------
        GenomePosition
        """
        return tuple(self.c2g(x) for x in self.p2c(ppos))


if __name__ == "__main__":
    # The following exons are from AB026906.1.
    # test case: g.7872 -> c.274 -> p.92
    # N.B. These are python counting coordinates (0-based)
    exons = [(5808, 5860), (6757, 6874), (7767, 7912), (13709, 13785)]

    def test_list(g_range):
        """Test a list."""
        cm = CoordinateMapper(exons)
        for g1 in g_range:
            print(g1, end=" ")
            c1 = cm.g2c(g1)
            print(c1, end=" ")
            p1 = cm.c2p(c1)
            print(p1, end=" ")
            if p1:
                c2 = cm.p2c(p1)[0]
            else:
                c2 = c1
            print(" | ", c2, end=" ")
            g2 = cm.c2g(c2)
            print(g2)

        print(cm.g2p(7872))
        print(cm.p2g(92))

    def test_simple():
        """Simple test."""
        from Bio.SeqFeature import SeqFeature

        location = sum(
            [
                FeatureLocation(2, 4, +1),
                FeatureLocation(8, 11, +1),
                FeatureLocation(16, 18, +1),
            ]
        )
        simple_exons = SeqFeature(location, type="CDS")
        cm = CoordinateMapper(simple_exons)
        print(cm.exons)
        print(list(cm.exons))
        print(range(len(cm.exons)))
        for i in range(len(cm.exons)):
            print("%3s" % cm.c2g(i), end=" ")
        print()
        for i in range(20):
            print("%3d" % i, end=" ")
        print()
        for i in range(20):
            print("%3s" % cm.g2c(i), end=" ")
        print()

    r1 = (7870, 7871, 7872, 7873, 7874)
    r2 = (
        5807,
        5808,
        5809,
        6871,
        6872,
        6873,
        6874,
        6875,
        7766,
        7767,
        7768,
        7769,
        13784,
        13785,
        13786,
    )
    test_list(r1)
    # test_list(r2)
    print()
    test_simple()
