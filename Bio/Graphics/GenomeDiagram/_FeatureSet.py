# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Contact:       Leighton Pritchard, Scottish Crop Research Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                L.Pritchard@scri.ac.uk
################################################################################
#
# Thanks to Peter Cock for the impetus to write the get_features() code to
# subselect Features.
#
################################################################################

""" FeatureSet module

    Provides:

    o FeatureSet - container for Feature objects

    For drawing capabilities, this module uses reportlab to draw and write
    the diagram:

    http://www.reportlab.com

    For dealing with biological information, the package expects BioPython
    objects:

    http://www.biopython.org
"""

#------------------------------------------------------------------------------
# IMPORTS

# ReportLab
from reportlab.pdfbase import _fontdata
from reportlab.lib import colors

# GenomeDiagram
from _Feature import Feature

# Builtins
import re

#------------------------------------------------------------------------------
# CLASSES

#------------------------------------------------------------
# FeatureSet

class FeatureSet(object):
    """ FeatureSet

        Provides:

        Methods:

        o __init__(self, set_id=None, name=None) Called on instantiation

        o add_feature(self, feature, color=colors.lightgreen)  Add a Feature
                        object to the set

        o del_feature(self, feature_id) Remove a feature from the set, by id

        o set_all_features(self, attr, value)   Set the passed attribute to the
                        passed value in all features in the set

        o get_features(self)    Returns a list of Features from the set
    
        o get_ids(self)     Returns a list of unique ids for features in the set

        o range(self)       Returns the range of bases covered by features in 
                            the set

        o to_string(self, verbose=0)    Returns a string describing the set

        o __len__(self)     Returns the length of sequence covered by the set

        o __getitem__(self, key)    Returns a feature from the set, keyed by id

        o __str__(self)     Returns a string describing the set

        Attributes:

        o id    Unique id for the set

        o name  String describing the set
    """
    def __init__(self, set_id=None, name=None, parent=None):
        """ __init__(self, set_id=None, name=None)

            o set_id    Unique id for the set

            o name      String identifying the feature set
        """
        self.parent = parent
        self.id = id            # Unique id for the set
        self.next_id = 0       # counter for unique feature ids
        self.features = {}     # Holds features, keyed by ID
        self.name = name        # String describing the set


    def add_feature(self, feature, **kwargs):
        """ add_feature(self, feature, **args)

            o feature       Bio.SeqFeature object

            o **kwargs      Keyword arguments for Feature.  Named attributes
                            of the Feature
                                                        

            Add a Bio.SeqFeature object to the diagram (will be stored
            internally in a Feature wrapper
        """
        id = self.next_id                                  # get id number
        f = Feature(self, id, feature)
        self.features[id] = f # add feature
        for key in kwargs:
            if key == "colour" or key == "color":
                #Deal with "colour" as a special case by also mapping to color.
                #If Feature.py used a python property we wouldn't need to call
                #set_color explicitly.  However, this is important to make sure
                #every color gets mapped to a colors object - for example color
                #numbers, or strings (may not matter for PDF, but does for PNG).
                self.features[id].set_color(kwargs[key])
                continue
            setattr(self.features[id], key, kwargs[key])
        self.next_id += 1                                  # increment next id
        return f

    def del_feature(self, feature_id):
        """ del_feature(self, feature_id)

            o feature_id        Unique id of the feature to delete

            Remove a feature from the set, indicated by its id
        """
        del self.features[feature_id]


    def set_all_features(self, attr, value):
        """ set_all_features(self, attr, value)

            o attr      An attribute of the Feature class

            o value     The value to set that attribute

            Set the passed attribute of all features in the set to the
            passed value
        """
        changed = 0
        for feature in self.features.values():
            # If the feature has the attribute, and the value should change
            if hasattr(feature, attr):    
                if getattr(feature, attr) != value:
                    setattr(feature, attr, value) # set it to the passed value

        #For backwards compatibility, we support both colour and color.
        #As a quick hack, make "colour" set both "colour" and "color".
        #if attr=="colour":
        #    self.set_all_feature("color",value)


    def get_features(self, attribute=None, value=None, comparator=None):
        """ get_features(self, attribute=None, value=None, comparator=None) ->
                                            [Feature, Feature, ...]

            o attribute        String, attribute of a Feature object

            o value            The value desired of the attribute

            o comparator       String, how to compare the Feature attribute to the
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
            return self.features.values()
        # If no comparator is specified, return all features where the attribute
        # value matches that passed
        if comparator is None:
            return [feature for feature in self.features.values() if\
                    getattr(feature, attribute) == value]
        # If the comparator is 'not', return all features where the attribute
        # value does not match that passed
        elif comparator == 'not':
            return [feature for feature in self.features.values() if\
                    getattr(feature, attribute) != value]
        # If the comparator is 'startswith', return all features where the attribute
        # value does not match that passed
        elif comparator == 'startswith':
            return [feature for feature in self.features.values() if\
                    getattr(feature, attribute).startswith(value)]
        # If the comparator is 'like', use a regular expression search to identify
        # features
        elif comparator == 'like':
            return [feature for feature in self.features.values() if\
                    re.search(value, getattr(feature, attribute))]
        # As a final option, just return an empty list
        return []



    def get_ids(self):
        """ get_ids(self) -> [int, int, ...]

            Return a list of all ids for the feature set
        """
        return self.features.keys()


    def range(self):
        """ range(self)

            Returns the lowest and highest base (or mark) numbers as a tuple
        """
        lows, highs = [], []
        for feature in self.features.values():
            for start, end in feature.locations:
                lows.append(start)
                highs.append(end)
        if len(lows) != 0 and len(highs) != 0:      # Default in case there is 
            return (min(lows), max(highs))          # nothing in the set
        return 0, 0


    def to_string(self, verbose=0):
        """ to_string(self, verbose=0) -> ""

            o verbose       Boolean indicating whether a short or complete 
                            account of the set is required

            Returns a formatted string with information about the set
        """
        if not verbose:         # Short account only required
            return "%s" % self
        else:                   # Long account desired
            outstr = ["\n<%s: %s>" % (self.__class__, self.name)]
            outstr.append("%d features" % len(self.features))
            for key in self.features:
                outstr.append("feature: %s" % self.features[key])
            return "\n".join(outstr)

    def __len__(self):
        """ __len__(self) -> int

            Return the number of features in the set
        """
        return len(self.features)


    def __getitem__(self, key):
        """ __getitem__(self, key) -> Feature

            Return a feature, keyed by id
        """
        return self.features[key]


    def __str__(self):
        """ __str__(self) -> ""

            Returns a formatted string with information about the feature set
        """
        outstr = ["\n<%s: %s %d features>" % (self.__class__, self.name, 
                                              len(self.features))]
        return "\n".join(outstr)

################################################################################
# RUN AS SCRIPT
################################################################################

if __name__ == '__main__':

    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature
    
    genbank_entry = SeqIO.read('/data/Genomes/Bacteria/Nanoarchaeum_equitans/NC_005213.gbk', 'gb')

    # Test code
    gdfs = FeatureSet(0, 'Nanoarchaeum equitans CDS')
    for feature in genbank_entry.features:
        if feature.type == 'CDS':
            gdfs.add_feature(feature)

    #print len(gdfs)
    #print gdfs.get_ids()
    #gdfs.del_feature(560)
    #print gdfs.get_ids()
    #print gdfs.get_features()
    #for feature in gdfs.get_features():
    #    print feature.id, feature.start, feature.end
    #print gdfs[500]   

    
