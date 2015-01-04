# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Abstract base classes for the SearchIO object model."""


from Bio._utils import getattr_str, trim_str


__docformat__ = "restructuredtext en"


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
        for attr in self.__dict__:
            if attr not in self._NON_STICKY_ATTRS:
                setattr(obj, attr, self.__dict__[attr])
