# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Common SearchIO utility functions."""


def getattr_str(obj, attr, fmt=None, fallback="?"):
    """Return string of the given object's attribute.

    Defaults to the given fallback value if attribute is not present.
    """
    try:
        value = getattr(obj, attr)
    except AttributeError:
        return fallback
    if fmt is None:
        return str(value)
    return fmt % value


def read_forward(handle):
    """Read through whitespaces, return the first non-whitespace line."""
    while True:
        line = handle.readline()
        # if line is empty or line has characters and stripping does not remove
        # them, return the line
        if (not line) or (line and line.strip()):
            return line


def get_processor(format, mapping):
    """Return the object to process the given format according to the mapping.

    :param format: format name
    :type format: string, lower case
    :param mapping: mapping of format name and its processor object
    :type mapping: dictionary {string: object}

    """
    # map file format to iterator name
    try:
        obj_info = mapping[format]
    except KeyError:
        # handle the errors with helpful messages
        if format is None:
            raise ValueError("Format required (lower case string)") from None
        elif not isinstance(format, str):
            raise TypeError("Need a string for the file format (lower case)") from None
        elif format != format.lower():
            raise ValueError("Format string %r should be lower case" % format) from None
        else:
            raise ValueError(
                "Unknown format %r. Supported formats are %r"
                % (format, "', '".join(mapping))
            ) from None

    mod_name, obj_name = obj_info
    mod = __import__("Bio.SearchIO.%s" % mod_name, fromlist=[""])

    return getattr(mod, obj_name)


def singleitem(attr=None, doc=""):
    """Property for fetching attribute from first entry of container.

    Returns a property that fetches the given attribute from
    the first item in a SearchIO container object.
    """

    def getter(self):
        if len(self._items) > 1:
            raise ValueError("More than one HSPFragment objects found in HSP")
        if attr is None:
            return self._items[0]
        return getattr(self._items[0], attr)

    return property(fget=getter, doc=doc)


def allitems(attr=None, doc=""):
    """Property for fetching attribute from all entries of container.

    Returns a property that fetches the given attributes from
    all items in a SearchIO container object.
    """

    def getter(self):
        if attr is None:
            return self._items
        return [getattr(frag, attr) for frag in self._items]

    return property(fget=getter, doc=doc)


def fullcascade(attr, doc=""):
    """Return a getter property with a cascading setter.

    This is similar to ``optionalcascade``, but for SearchIO containers that have
    at least one item (HSP). The getter always retrieves the attribute
    value from the first item. If the items have more than one attribute values,
    an error will be raised. The setter behaves like ``partialcascade``, except
    that it only sets attributes to items in the object, not the object itself.

    """

    def getter(self):
        return getattr(self._items[0], attr)

    def setter(self, value):
        for item in self:
            setattr(item, attr, value)

    return property(fget=getter, fset=setter, doc=doc)


def optionalcascade(cont_attr, item_attr, doc=""):
    """Return a getter property with a cascading setter.

    This is used for the ``id`` and ``description`` properties of the container
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
            # as well (the underlying dict is not integer-indexable)
            return getattr(self[0], item_attr)
        else:
            return getattr(self, cont_attr)

    def setter(self, value):
        setattr(self, cont_attr, value)
        for item in self:
            setattr(item, item_attr, value)

    return property(fget=getter, fset=setter, doc=doc)


def fragcascade(attr, seq_type, doc=""):
    """Return a getter property with cascading setter, for HSPFragment objects.

    Similar to ``partialcascade``, but for HSPFragment objects and acts on ``query``
    or ``hit`` properties of the object if they are not None.

    """
    assert seq_type in ("hit", "query")
    attr_name = f"_{seq_type}_{attr}"

    def getter(self):
        return getattr(self, attr_name)

    def setter(self, value):
        setattr(self, attr_name, value)
        seq = getattr(self, seq_type)
        if seq is not None:
            setattr(seq, attr, value)

    return property(fget=getter, fset=setter, doc=doc)
