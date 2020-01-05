# Copyright 2011-2017 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Genome Diagram Feature cross-link module."""

from reportlab.lib import colors


class CrossLink:
    """Hold information for drawing a cross link between features."""

    def __init__(
        self, featureA, featureB, color=colors.lightgreen, border=None, flip=False
    ):
        """Create a new cross link.

        Arguments featureA and featureB should GenomeDiagram feature objects,
        or 3-tuples (track object, start, end), and currently must be on
        different tracks.

        The color and border arguments should be ReportLab colour objects, or
        for border use a boolean False for no border, otherwise it defaults to
        the same as the main colour.

        The flip argument draws an inverted cross link, useful for showing a
        mapping where one sequence has been reversed. It is conventional to
        also use a different colour (e.g. red for simple links, blue for any
        flipped links).
        """
        # Initialize attributes
        self.featureA = featureA
        self.featureB = featureB
        self.color = color  # default color to draw the feature
        self.border = border
        self.flip = flip

    @property
    def startA(self):
        """Start position of Feature A."""
        try:
            return self.featureA.start
        except AttributeError:
            track, start, end = self.featureA
            return start

    @property
    def endA(self):
        """End position of Feature A."""
        try:
            return self.featureA.end
        except AttributeError:
            track, start, end = self.featureA
            return end

    def _trackA(self, tracks):
        try:
            track, start, end = self.featureA
            assert track in tracks
            return track
        except TypeError:
            for track in tracks:
                for feature_set in track.get_sets():
                    if hasattr(feature_set, "features"):
                        if self.featureA in feature_set.features.values():
                            return track
            return None

    @property
    def startB(self):
        """Start position of Feature B."""
        try:
            return self.featureB.start
        except AttributeError:
            track, start, end = self.featureB
            return start

    @property
    def endB(self):
        """End position of Feature B."""
        try:
            return self.featureB.end
        except AttributeError:
            track, start, end = self.featureB
            return end

    def _trackB(self, tracks):
        try:
            track, start, end = self.featureB
            assert track in tracks
            return track
        except TypeError:
            for track in tracks:
                for feature_set in track.get_sets():
                    if hasattr(feature_set, "features"):
                        if self.featureB in feature_set.features.values():
                            return track
            return None
