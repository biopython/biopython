# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Abstract base classes for the SearchIO object model."""


class _BaseSearchObject(object):

    """Abstract class for SearchIO objects."""

    _NON_STICKY_ATTRS = ()

    def _transfer_attrs(self, obj):
        """Transfer instance attributes to the given object.

        This method is used to transfer attributes set externally (for example
        using `setattr`) to a new object created from this one (for example
        from slicing).

        The reason this method is necessary is because different parsers will
        set different attributes for each QueryResult, Hit, HSP, or HSPFragment
        objects, depending on the attributes they found in the search output
        file. Ideally, we want these attributes to 'stick' with any new instance
        object created from the original one.

        """
        # list of attribute names we don't want to transfer
        for attr in self.__dict__.keys():
            if attr not in self._NON_STICKY_ATTRS:
                setattr(obj, attr, self.__dict__[attr])

    def _trunc_display(string, max_len, concat_char):
        """Truncates the given string for display."""
        if len(string) > max_len:
            return string[:max_len - len(concat_char)] + concat_char
        return string

    _trunc_display = staticmethod(_trunc_display)

    def _attr_display(obj, attr, fmt=None, fallback='?'):
        """Returns a string of the given object's attribute."""
        if hasattr(obj, attr):
            if fmt is not None:
                return fmt % getattr(obj, attr)
            return str(getattr(obj, attr))
        return fallback

    _attr_display = staticmethod(_attr_display)


class _BaseHSP(_BaseSearchObject):

    """Abstract base class for HSP objects."""

    def _str_hsp_header(self):
        """Prints the alignment header info."""
        lines = []
        # set query id line
        qid_line = self._trunc_display('      Query: %s %s' %
                (self.query_id, self.query_description), 80, '...')
        # set hit id line
        hid_line = self._trunc_display('        Hit: %s %s' %
                (self.hit_id, self.hit_description), 80, '...')
        lines.append(qid_line)
        lines.append(hid_line)

        # coordinates
        query_start = _BaseHSP._attr_display(self, 'query_start')
        query_end = _BaseHSP._attr_display(self, 'query_end')
        hit_start = _BaseHSP._attr_display(self, 'hit_start')
        hit_end = _BaseHSP._attr_display(self, 'hit_end')

        # strands
        try:
            qstrand = self.query_strand
            hstrand = self.hit_strand
        except ValueError:
            qstrand = self.query_strands[0]
            hstrand = self.hit_strands[0]
        lines.append('Query range: [%s:%s] (%r)' % (query_start, query_end,
                qstrand))
        lines.append('  Hit range: [%s:%s] (%r)' % (hit_start, hit_end,
                hstrand))

        return '\n'.join(lines)
