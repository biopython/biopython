# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Common SearchIO utility functions."""

from Bio._py3k import basestring


def get_processor(format, mapping):
    """Returns the object to process the given format according to the mapping.

    Arguments:
    format -- Lower case string denoting one of the supported formats.
    mapping -- Dictionary of format and object name mapping.

    """
    # map file format to iterator name
    try:
        obj_info = mapping[format]
    except KeyError:
        # handle the errors with helpful messages
        if format is None:
            raise ValueError("Format required (lower case string)")
        elif not isinstance(format, basestring):
            raise TypeError("Need a string for the file format (lower case)")
        elif format != format.lower():
            raise ValueError("Format string %r should be lower case" %
                    format)
        else:
            raise ValueError("Unknown format %r. Supported formats are "
                    "%r" % (format, "', '".join(mapping)))

    mod_name, obj_name = obj_info
    mod = __import__('Bio.SearchIO.%s' % mod_name, fromlist=[''])

    return getattr(mod, obj_name)


def singleitem(attr=None, doc=''):
    """Returns a property that fetches the given attribute from
    the first item in a SearchIO container object."""
    def getter(self):
        if len(self._items) > 1:
            raise ValueError("More than one HSPFragment objects "
                    "found in HSP")
        if attr is None:
            return self._items[0]
        return getattr(self._items[0], attr)
    return property(fget=getter, doc=doc)


def allitems(attr=None, doc=''):
    """Returns a property that fetches the given attributes from
    all items in a SearchIO container object."""
    def getter(self):
        if attr is None:
            return self._items
        return [getattr(frag, attr) for frag in self._items]
    return property(fget=getter, doc=doc)


def fullcascade(attr, doc=''):
    """Returns a getter property with a cascading setter.

    This is similar to `optionalcascade`, but for SearchIO containers that have
    at least one item (HSP). The getter always retrieves the attribute
    value from the first item. If the items have more than one attribute values,
    an error will be raised. The setter behaves like `partialcascade`, except
    that it only sets attributes to items in the object, not the object itself.

    """
    def getter(self):
        return getattr(self._items[0], attr)

    def setter(self, value):
        for item in self:
            setattr(item, attr, value)

    return property(fget=getter, fset=setter, doc=doc)


def optionalcascade(cont_attr, item_attr, doc=''):
    """Returns a getter property with a cascading setter.

    This is used for the `id` and `description` properties of the container
    objects with zero or more items. These items have their own private
    attributes that stores query and/or hit ID and description. When the
    container has zero items, attribute values are always retrieved from the
    container's attribute. Otherwise, the first item's attribute is used.

    To keep the container items' query and/or hit ID and description in-sync,
    the setter cascades any new value given to the items' values.

    """
    def getter(self):
        if self._items:
            # don't use self._items here, so QueryResult can use this property
            # as well (the underlying OrderedDict is not integer-indexable)
            return getattr(self[0], item_attr)
        else:
            return getattr(self, cont_attr)

    def setter(self, value):
        setattr(self, cont_attr, value)
        for item in self:
            setattr(item, item_attr, value)

    return property(fget=getter, fset=setter, doc=doc)


def fragcascade(attr, seq_type, doc=''):
    """Returns a getter property with cascading setter, for HSPFragment objects.

    Similar to `partialcascade`, but for HSPFragment objects and acts on `query`
    or `hit` properties of the object if they are not None.

    """
    assert seq_type in ('hit', 'query')
    attr_name = '_%s_%s' % (seq_type, attr)

    def getter(self):
        return getattr(self, attr_name)

    def setter(self, value):
        setattr(self, attr_name, value)
        seq = getattr(self, seq_type)
        if seq is not None:
            setattr(seq, attr, value)

    return property(fget=getter, fset=setter, doc=doc)
