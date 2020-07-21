# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Abstract base classes for the SearchIO object model."""


from Bio.SearchIO._utils import getattr_str


class _BaseSearchObject:
    """Abstract class for SearchIO objects."""

    _NON_STICKY_ATTRS = ()

    def _transfer_attrs(self, obj):
        """Transfer instance attributes to the given object (PRIVATE).

        This method is used to transfer attributes set externally (for example
        using ``setattr``) to a new object created from this one (for example
        from slicing).

        The reason this method is necessary is because different parsers will
        set different attributes for each QueryResult, Hit, HSP, or HSPFragment
        objects, depending on the attributes they found in the search output
        file. Ideally, we want these attributes to 'stick' with any new instance
        object created from the original one.

        """
        # list of attribute names we don't want to transfer
        for attr in self.__dict__:
            if attr not in self._NON_STICKY_ATTRS:
                setattr(obj, attr, self.__dict__[attr])


class _BaseHSP(_BaseSearchObject):
    """Abstract base class for HSP objects."""

    def _str_hsp_header(self):
        """Print the alignment header info (PRIVATE)."""
        lines = []
        # set query id line
        qid_line = "      Query: %s %s" % (self.query_id, self.query_description)
        qid_line = qid_line[:77] + "..." if len(qid_line) > 80 else qid_line
        # set hit id line
        hid_line = "        Hit: %s %s" % (self.hit_id, self.hit_description)
        hid_line = hid_line[:77] + "..." if len(hid_line) > 80 else hid_line
        lines.append(qid_line)
        lines.append(hid_line)

        # coordinates
        query_start = getattr_str(self, "query_start")
        query_end = getattr_str(self, "query_end")
        hit_start = getattr_str(self, "hit_start")
        hit_end = getattr_str(self, "hit_end")

        # strands
        try:
            qstrand = self.query_strand
            hstrand = self.hit_strand
        except ValueError:
            qstrand = self.query_strand_all[0]
            hstrand = self.hit_strand_all[0]
        lines.append("Query range: [%s:%s] (%r)" % (query_start, query_end, qstrand))
        lines.append("  Hit range: [%s:%s] (%r)" % (hit_start, hit_end, hstrand))

        return "\n".join(lines)
