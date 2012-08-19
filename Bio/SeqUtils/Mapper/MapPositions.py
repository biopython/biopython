# Copyright 2012 Lenna X. Peterson <arklenna@gmail.com>

# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
from copy import copy
import re


class InvalidPositionError(ValueError):
    "Exception for bad coordinates"


class GenomePositionError(InvalidPositionError):
    "Exception for bad genome coordinates"


class CDSPositionError(InvalidPositionError):
    "Exception for bad CDS coordinates"


class ProteinPositionError(InvalidPositionError):
    "Exception for bad protein coordinates"


class MapPosition(object):
    """Generic position for coordinate mapping."""
    def __init__(self, pos, index=None, **kwargs):
        """Init class from HGVS position.

        Parameters
        ----------
        pos : int, str
            Position to convert
        index : int
            Start of counting (e.g. 1 for GenBank or HGVS)
        """
        self.index = index
        if pos and self.index:
            pos -= index
        self.pos = pos

    @classmethod
    def from_dialect(cls, dialect, pos, strand=None):
        """Init class from given dialect.

        Parameters
        ----------
        dialect : str
            input dialect (HGVS or GenBank)
        pos : int, str
            dialect position to convert
        strand : int
            Strand of position (-1 or +1, default None)

        Returns
        -------
        cls
        """
        if dialect is None:
            return cls(pos, strand=strand)
        try:
            from_func = getattr(cls, "from_" + dialect.lower())
        except AttributeError:
            raise ValueError("Dialect '%s' not valid" % dialect)
        return from_func(pos, strand)

    @classmethod
    def from_hgvs(cls, hgvs_pos, strand=None):
        """Init class from HGVS position.

        Parameters
        ----------
        hgvs_pos : int, str
            HGVS position to convert
        strand : int
            Strand of position (-1 or +1, default None)

        Returns
        -------
        cls
        """
        return cls(hgvs_pos, index=1, strand=strand, post_fmt="*%d")

    @classmethod
    def from_genbank(cls, gbk_pos, strand=None):
        """Init class from GenBank position.

        Parameters
        ----------
        gbk_pos : int, str
            GenBank position to convert
        strand : int
            Strand of position (-1 or +1, default None)

        Returns
        -------
        cls
        """
        return cls(gbk_pos, index=1, strand=strand)

    def to(self, dialect):
        """Convert position to specified dialect.

        Parameters
        ----------
        dialect : str
            Output dialect (HGVS or GenBank)
        """
        if dialect is None:
            return self.to_str()
        try:
            to_func = getattr(self, "to_" + dialect.lower())
        except AttributeError:
            raise ValueError("Dialect '%s' not valid" % dialect)
        return to_func()

    def to_hgvs(self):
        """Convert position to HGVS"""
        if self.pos or self.pos == 0:
            return self.pos + 1
        return None

    def to_genbank(self):
        """Convert position to GenBank"""
        if self.pos or self.pos == 0:
            return self.pos + 1
        return None

    def to_str(self):
        """Make string representation without conversion"""
        return self.pos

    def __str__(self):
        """String representation"""
        return str(self.to_str())

    def __int__(self):
        """Integer representation"""
        return self.pos

    def __repr__(self):
        """Detailed representation for debugging"""
        return "%s(%s)" % (self.__class__.__name__, self.to_str())


class GenomePosition(MapPosition):
    """Genome position for coordinate mapping"""
    def __init__(self, gpos, index=None, strand=None, **kwargs):
        # FIXME if index is string, error may be raised
        if gpos < (index or 0):
            raise GenomePositionError("Genome position cannot be negative.")
        # call superclass constructor
        MapPosition.__init__(self, gpos, index)
        self.strand = strand

    def __eq__(self, other):
        """Compare equal to other GenomePosition with same pos
        or integer equal to pos"""
        if isinstance(other, int):
            return self.pos == other
        return isinstance(other, GenomePosition) and self.pos == other.pos


class ProteinPosition(MapPosition):
    """Protein position for coordinate mapping"""
    def __init__(self, ppos, index=None, **kwargs):
        # call superclass constructor
        MapPosition.__init__(self, ppos, index)

    def __eq__(self, other):
        """Compare equal to other ProteinPosition with same pos
        or integer equal to pos"""
        if isinstance(other, int):
            return self.pos == other
        return isinstance(other, ProteinPosition) and self.pos == other.pos


class CDSPosition(MapPosition):
    """CDS position for coordinate mapping"""
    def __init__(self, cpos, index=None,
                 pre_fmt=None, post_fmt=None, **kwargs):
        # Dispatch types and return anchor, offset
        if isinstance(cpos, int):
            anchor, offset = self.parse_int(cpos)
        elif isinstance(cpos, str):
            anchor, offset = self.parse_str(cpos, pre_fmt, post_fmt)
        else:
            raise CDSPositionError("'%s' is of unknown type" % cpos)
        # Set instance anchor and offset
        # call superclass constructor
        MapPosition.__init__(self, anchor, index)
        self._offset = offset
        self.validate()

    @classmethod
    def from_anchor(cls, anchor=None, offset=None):
        """Init CDSPosition with anchor, offset pair.

        Parameters
        ----------
        anchor : int
            CDS anchor (coordinate of nearest exon position)
        offset : int
            Offset from nearest exon position

        Returns
        -------
        CDSPosition
        """
        if not anchor:
            pos = cls("%+d" % offset)
        elif anchor < 0:
            raise CDSPositionError("Anchor cannot be negative.")
        else:
            pos = cls(anchor)
            pos.offset = offset
        return pos

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, val):
        "Validate new offset, then update"
        self.validate(offset=val)
        self._offset = val

    @property
    def anchor(self):
        return self.pos

    @anchor.setter
    def anchor(self, val):
        "Validate new anchor, then update pos"
        self.validate(anchor=val)
        self.pos = val

    def validate(self, anchor="sentinel", offset="sentinel"):
        """Check whether anchor and offset yield a valid position.

        Parameters
        ----------
        anchor : int
            CDS anchor (coordinate of nearest exon position)
        offset : int
            Offset from nearest exon position

        Returns
        -------
        bool
        """

        if anchor == "sentinel":
            anchor = self.anchor
        if offset == "sentinel":
            offset = self.offset
        if offset == 0:
            raise CDSPositionError(
                "Offset may not be 0. For no offset, use None.")
        if not anchor and anchor != 0 and not offset:
            raise CDSPositionError(
                "At least one of pos or offset must be defined")
        if anchor and anchor < 0:
            raise CDSPositionError("CDS anchor may not be negative.")
        return True

    @property
    def pos_type(self):
        "Type of CDS position, dynamically determined from values"
        # inside CDS
        if self.pos or self.pos == 0:
            if not self.offset:
                return "exon"
            return "intron"
        # outside CDS
        elif self.offset > 0:
            return "post-CDS"
        return "pre-CDS"
        assert False  # all integers should return

    @property
    def sub_dict(self):
        if self.pos_type == 'intron':
            return {'pos': self.pos, 'offset': self.offset}
        if self.pos_type == 'exon':
            return {'pos': self.pos}
        if self.pos_type == 'post-CDS' or self.pos_type == 'pre-CDS':
            return {'offset': self.offset}

    fmt_dict = {
        'exon': "%d",
        'intron': "%d%+d",
        'post-CDS': "%+d",
        'pre-CDS': "%+d",
    }

    @staticmethod
    def _shift_index(pos_dict, idx):
        """Increment value of dict key 'pos' by given index.

        Parameters
        ----------
        pos_dict : dict
            Dictionary to search for 'pos'
        idx : int
            Index to add to pos_dict['pos']

        Returns
        -------
        dict
        """
        if 'pos' in pos_dict and idx:
            pos_dict['pos'] += idx
        return pos_dict

    def _make_str(self, val_dict=None, fmt_dict=None):
        """Retrieve format string and substitute values.

        Parameters
        ----------
        val_dict : dict
            Dictionary of values to substitute into string
        fmt_dict : dict
            Dictionary of format strings for each pos type

        Returns
        -------
        str
        """
        # set default dicts if parameter dicts are false
        fmt_dict = fmt_dict or self.fmt_dict
        val_dict = val_dict or self.sub_dict
        return fmt_dict[self.pos_type] % tuple(val_dict.values())

    @staticmethod
    def parse_int(cpos):
        """Parse int to anchor, offset pair

        Parameters
        ----------
        cpos : int
            Integer position to convert

        Returns
        -------
        tuple
        """
        # treat negative int as offset
        if cpos < 0:
            return (None, cpos)
        # treat positive int as anchor
        return (cpos, None)

    @staticmethod
    def parse_str(cpos, pre_fmt, post_fmt):
        """Parse string to anchor, offset pair

        Parameters
        ----------
        cpos : str
            String position to convert
        pre_fmt : str
            Format string for pre-CDS positions
        post_fmt : str
            Format string for post-CDS positions

        Returns
        -------
        tuple
        """
        delimiters = "\+\-"
        if post_fmt and "*" in post_fmt:
            delimiters += "\*"
        # parenth causes split pattern to be kept
        delim_rx = re.compile("([%s])" % delimiters)
        parsed = delim_rx.split(cpos, 1)
        if len(parsed) == 1:  # may be int
            return CDSPosition.parse_int(int(cpos))
        # 1 split is normally length 2 but delimiter is also kept
        elif len(parsed) != 3:
            raise CDSPositionError(
                "String '%s' not parseable for this position." % cpos)
        if parsed[0] == "":
            anchor = None
        else:
            anchor = int(parsed[0])
        if post_fmt and parsed[1] == "*" == post_fmt[0]:
            parsed[1] = "+"
        offset = int("".join(parsed[1:]))
        return (anchor, offset)

    def to_hgvs(self):
        """Convert CDS position to HGVS"""
        fmt_dict = copy(self.fmt_dict)
        fmt_dict['post-CDS'] = "*%d"
        sub_dict = self._shift_index(self.sub_dict, 1)
        return self._make_str(sub_dict, fmt_dict)

    def to_genbank(self):
        """Convert CDS position to GenBank"""
        sub_dict = self._shift_index(self.sub_dict, 1)
        return self._make_str(sub_dict)

    def to_str(self):
        """Make string representation of CDS position"""
        return self._make_str()

    def __int__(self):
        """Integer representation of CDS exon, otherwise NotImplemented"""
        if self.pos_type == "exon":
            return MapPosition.__int__(self)
        return NotImplemented

    def __eq__(self, other):
        """Compare equal to other MapPosition with same pos and offset
        or int if exon"""
        if isinstance(other, int) and self.pos_type == "exon":
            return self.pos == other
        return isinstance(other, CDSPosition) and \
                self.pos == other.pos and \
                self.offset == other.offset


if __name__ == "__main__":

    def print_pos(pos_obj):
        print "object:", pos_obj
        print "repr:", repr(pos_obj)
        print "HGVS:", pos_obj.to_hgvs()
        print

    g = GenomePosition.from_hgvs(6)
    print_pos(g)

    test_g = GenomePosition(5)
    #test_c = CDSPosition("6+1")
    #test_c = CDSPosition.from_hgvs("6+1")
    test_c = CDSPosition.from_hgvs("*1")
    #test_c = CDSPosition(6)
    #test_c = CDSPosition(-1)
    print_pos(test_g)
    print test_c.pos_type
    print_pos(test_c)
