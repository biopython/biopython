# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# Contact:       Leighton Pritchard, The James Hutton Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                Leighton.Pritchard@hutton.ac.uk
################################################################################
#
# Thanks to Peter Cock for the impetus to write the get_features() code to
# subselect Features.
#
################################################################################

"""FeatureSet module.

Provides:
 - FeatureSet - container for Feature objects

For drawing capabilities, this module uses reportlab to draw and write
the diagram: http://www.reportlab.com
"""


# GenomeDiagram
from ._Feature import Feature

# Builtins
import re


class FeatureSet:
    """FeatureSet object."""

    def __init__(self, set_id=None, name=None, parent=None):
        """Create the object.

        Arguments:
         - set_id: Unique id for the set
         - name: String identifying the feature set

        """
        self.parent = parent
        self.id = id  # Unique id for the set
        self.next_id = 0  # counter for unique feature ids
        self.features = {}  # Holds features, keyed by ID
        self.name = name  # String describing the set

    def add_feature(self, feature, **kwargs):
        """Add a new feature.

        Arguments:
         - feature: Bio.SeqFeature object
         - kwargs: Keyword arguments for Feature.  Named attributes
           of the Feature

        Add a Bio.SeqFeature object to the diagram (will be stored
        internally in a Feature wrapper).
        """
        id = self.next_id  # get id number
        f = Feature(self, id, feature)
        self.features[id] = f  # add feature
        for key in kwargs:
            if key == "colour" or key == "color":
                # Deal with "colour" as a special case by also mapping to color.
                # If Feature.py used a python property we wouldn't need to call
                # set_color explicitly.  However, this is important to make sure
                # every color gets mapped to a colors object - for example color
                # numbers, or strings (may not matter for PDF, but does for PNG).
                self.features[id].set_color(kwargs[key])
                continue
            setattr(self.features[id], key, kwargs[key])
        self.next_id += 1  # increment next id
        return f

    def del_feature(self, feature_id):
        """Delete a feature.

        Arguments:
         - feature_id: Unique id of the feature to delete

        Remove a feature from the set, indicated by its id.
        """
        del self.features[feature_id]

    def set_all_features(self, attr, value):
        """Set an attribute of all the features.

        Arguments:
         - attr: An attribute of the Feature class
         - value: The value to set that attribute to

        Set the passed attribute of all features in the set to the
        passed value.
        """
        for feature in self.features.values():
            if hasattr(feature, attr):
                # If the feature has the attribute, set it to the passed value
                setattr(feature, attr, value)

        # For backwards compatibility, we support both colour and color.
        # As a quick hack, make "colour" set both "colour" and "color".
        # if attr=="colour":
        #    self.set_all_feature("color",value)

    def get_features(self, attribute=None, value=None, comparator=None):
        """Retrieve features.

        Arguments:
         - attribute: String, attribute of a Feature object
         - value: The value desired of the attribute
         - comparator: String, how to compare the Feature attribute to the
           passed value

        If no attribute or value is given, return a list of all features in the
        feature set.  If both an attribute and value are given, then depending
        on the comparator, then a list of all features in the FeatureSet
        matching (or not) the passed value will be returned.  Allowed comparators
        are: 'startswith', 'not', 'like'.

        The user is expected to make a responsible decision about which feature
        attributes to use with which passed values and comparator settings.
        """
        # If no attribute or value specified, return all features
        if attribute is None or value is None:
            return list(self.features.values())
        # If no comparator is specified, return all features where the attribute
        # value matches that passed
        if comparator is None:
            return [
                feature
                for feature in self.features.values()
                if getattr(feature, attribute) == value
            ]
        # If the comparator is 'not', return all features where the attribute
        # value does not match that passed
        elif comparator == "not":
            return [
                feature
                for feature in self.features.values()
                if getattr(feature, attribute) != value
            ]
        # If the comparator is 'startswith', return all features where the attribute
        # value does not match that passed
        elif comparator == "startswith":
            return [
                feature
                for feature in self.features.values()
                if getattr(feature, attribute).startswith(value)
            ]
        # If the comparator is 'like', use a regular expression search to identify
        # features
        elif comparator == "like":
            return [
                feature
                for feature in self.features.values()
                if re.search(value, getattr(feature, attribute))
            ]
        # As a final option, just return an empty list
        return []

    def get_ids(self):
        """Return a list of all ids for the feature set."""
        return list(self.features.keys())

    def range(self):
        """Return the lowest and highest base (or mark) numbers as a tuple."""
        lows, highs = [], []
        for feature in self.features.values():
            for start, end in feature.locations:
                lows.append(start)
                highs.append(end)
        if len(lows) != 0 and len(highs) != 0:  # Default in case there is
            return (min(lows), max(highs))  # nothing in the set
        return 0, 0

    def to_string(self, verbose=0):
        """Return a formatted string with information about the set.

        Arguments:
         - verbose: Boolean indicating whether a short (default) or
           complete account of the set is required

        """
        if not verbose:  # Short account only required
            return f"{self}"
        else:  # Long account desired
            outstr = [f"\n<{self.__class__}: {self.name}>"]
            outstr.append("%d features" % len(self.features))
            for key in self.features:
                outstr.append(f"feature: {self.features[key]}")
            return "\n".join(outstr)

    def __len__(self):
        """Return the number of features in the set."""
        return len(self.features)

    def __getitem__(self, key):
        """Return a feature, keyed by id."""
        return self.features[key]

    def __str__(self):
        """Return a formatted string with information about the feature set."""
        outstr = [
            "\n<%s: %s %d features>" % (self.__class__, self.name, len(self.features))
        ]
        return "\n".join(outstr)
