# Copyright 2001, 2003 by Brad Chapman.  All rights reserved.
# Revisions copyright 2011 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Draw representations of organism chromosomes with added information.

These classes are meant to model the drawing of pictures of chromosomes.
This can be useful for lots of things, including displaying markers on
a chromosome (ie. for genetic mapping) and showing syteny between two
chromosomes.

The structure of these classes is intended to be a Composite, so that
it will be easy to plug in and switch different parts without
breaking the general drawing capabilities of the system. The
relationship between classes is that everything derives from
_ChromosomeComponent, which specifies the overall interface. The parts
then are related so that an Organism contains Chromosomes, and these
Chromosomes contain ChromosomeSegments. This representation differents
from the canonical composite structure in that we don't really have
'leaf' nodes here -- all components can potentially hold sub-components.

Most of the time the ChromosomeSegment class is what you'll want to
customize for specific drawing tasks.

For providing drawing capabilities, these classes use reportlab:

http://www.reportlab.com

This provides nice output in PDF, SVG and postscript.  If you have
reportlab's renderPM module installed you can also use PNG etc.
"""

# reportlab
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.pdfbase.pdfmetrics import stringWidth

from reportlab.graphics.shapes import Drawing, String, Line, Rect, Wedge, ArcPath
from reportlab.graphics.widgetbase import Widget

from Bio.Graphics import _write
from Bio.Graphics.GenomeDiagram import _Colors


_color_trans = _Colors.ColorTranslator()


class _ChromosomeComponent(Widget):
    """Base class specifying the interface for a component of the system.

    This class should not be instantiated directly, but should be used
    from derived classes.
    """

    def __init__(self):
        """Initialize a chromosome component.

        Attributes:
        - _sub_components -- Any components which are contained under
        this parent component. This attribute should be accessed through
        the add() and remove() functions.

        """
        self._sub_components = []

    def add(self, component):
        """Add a sub_component to the list of components under this item."""
        if not isinstance(component, _ChromosomeComponent):
            raise TypeError(
                "Expected a _ChromosomeComponent object, got %s" % component
            )

        self._sub_components.append(component)

    def remove(self, component):
        """Remove the specified component from the subcomponents.

        Raises a ValueError if the component is not registered as a
        sub_component.
        """
        try:
            self._sub_components.remove(component)
        except ValueError:
            raise ValueError(
                "Component %s not found in sub_components." % component
            ) from None

    def draw(self):
        """Draw the specified component."""
        raise AssertionError("Subclasses must implement.")


class Organism(_ChromosomeComponent):
    """Top level class for drawing chromosomes.

    This class holds information about an organism and all of it's
    chromosomes, and provides the top level object which could be used
    for drawing a chromosome representation of an organism.

    Chromosomes should be added and removed from the Organism via the
    add and remove functions.
    """

    def __init__(self, output_format="pdf"):
        """Initialize."""
        _ChromosomeComponent.__init__(self)

        # customizable attributes
        self.page_size = letter
        self.title_size = 20

        # Do we need this given we don't draw a legend?
        # If so, should be a public API...
        self._legend_height = 0  # 2 * inch

        self.output_format = output_format

    def draw(self, output_file, title):
        """Draw out the information for the Organism.

        Arguments:
         - output_file -- The name of a file specifying where the
           document should be saved, or a handle to be written to.
           The output format is set when creating the Organism object.
           Alternatively, output_file=None will return the drawing using
           the low-level ReportLab objects (for further processing, such
           as adding additional graphics, before writing).
         - title -- The output title of the produced document.

        """
        width, height = self.page_size
        cur_drawing = Drawing(width, height)

        self._draw_title(cur_drawing, title, width, height)

        cur_x_pos = inch * 0.5
        if len(self._sub_components) > 0:
            x_pos_change = (width - inch) / len(self._sub_components)
        # no sub_components
        else:
            pass

        for sub_component in self._sub_components:
            # set the drawing location of the chromosome
            sub_component.start_x_position = cur_x_pos + 0.05 * x_pos_change
            sub_component.end_x_position = cur_x_pos + 0.95 * x_pos_change
            sub_component.start_y_position = height - 1.5 * inch
            sub_component.end_y_position = self._legend_height + 1 * inch

            # do the drawing
            sub_component.draw(cur_drawing)

            # update the locations for the next chromosome
            cur_x_pos += x_pos_change

        self._draw_legend(cur_drawing, self._legend_height + 0.5 * inch, width)

        if output_file is None:
            # Let the user take care of writing to the file...
            return cur_drawing

        return _write(cur_drawing, output_file, self.output_format)

    def _draw_title(self, cur_drawing, title, width, height):
        """Write out the title of the organism figure (PRIVATE)."""
        title_string = String(width / 2, height - inch, title)
        title_string.fontName = "Helvetica-Bold"
        title_string.fontSize = self.title_size
        title_string.textAnchor = "middle"

        cur_drawing.add(title_string)

    def _draw_legend(self, cur_drawing, start_y, width):
        """Draw a legend for the figure (PRIVATE).

        Subclasses should implement this (see also self._legend_height) to
        provide specialized legends.
        """
        pass


class Chromosome(_ChromosomeComponent):
    """Class for drawing a chromosome of an organism.

    This organizes the drawing of a single organisms chromosome. This
    class can be instantiated directly, but the draw method makes the
    most sense to be called in the context of an organism.
    """

    def __init__(self, chromosome_name):
        """Initialize a Chromosome for drawing.

        Arguments:
         - chromosome_name - The label for the chromosome.

        Attributes:
         - start_x_position, end_x_position - The x positions on the page
           where the chromosome should be drawn. This allows multiple
           chromosomes to be drawn on a single page.
         - start_y_position, end_y_position - The y positions on the page
           where the chromosome should be contained.

        Configuration Attributes:
         - title_size - The size of the chromosome title.
         - scale_num - A number of scale the drawing by. This is useful if
           you want to draw multiple chromosomes of different sizes at the
           same scale. If this is not set, then the chromosome drawing will
           be scaled by the number of segements in the chromosome (so each
           chromosome will be the exact same final size).

        """
        _ChromosomeComponent.__init__(self)

        self._name = chromosome_name

        self.start_x_position = -1
        self.end_x_position = -1
        self.start_y_position = -1
        self.end_y_position = -1

        self.title_size = 20
        self.scale_num = None

        self.label_size = 6
        self.chr_percent = 0.25
        self.label_sep_percent = self.chr_percent * 0.5
        self._color_labels = False

    def subcomponent_size(self):
        """Return the scaled size of all subcomponents of this component."""
        total_sub = 0
        for sub_component in self._sub_components:
            total_sub += sub_component.scale

        return total_sub

    def draw(self, cur_drawing):
        """Draw a chromosome on the specified template.

        Ideally, the x_position and y_*_position attributes should be
        set prior to drawing -- otherwise we're going to have some problems.
        """
        for position in (
            self.start_x_position,
            self.end_x_position,
            self.start_y_position,
            self.end_y_position,
        ):
            assert position != -1, "Need to set drawing coordinates."

        # first draw all of the sub-sections of the chromosome -- this
        # will actually be the picture of the chromosome
        cur_y_pos = self.start_y_position
        if self.scale_num:
            y_pos_change = (
                self.start_y_position * 0.95 - self.end_y_position
            ) / self.scale_num
        elif len(self._sub_components) > 0:
            y_pos_change = (
                self.start_y_position * 0.95 - self.end_y_position
            ) / self.subcomponent_size()
        # no sub_components to draw
        else:
            pass

        left_labels = []
        right_labels = []
        for sub_component in self._sub_components:
            this_y_pos_change = sub_component.scale * y_pos_change

            # set the location of the component to draw
            sub_component.start_x_position = self.start_x_position
            sub_component.end_x_position = self.end_x_position
            sub_component.start_y_position = cur_y_pos
            sub_component.end_y_position = cur_y_pos - this_y_pos_change

            # draw the sub component
            sub_component._left_labels = []
            sub_component._right_labels = []
            sub_component.draw(cur_drawing)
            left_labels += sub_component._left_labels
            right_labels += sub_component._right_labels

            # update the position for the next component
            cur_y_pos -= this_y_pos_change

        self._draw_labels(cur_drawing, left_labels, right_labels)
        self._draw_label(cur_drawing, self._name)

    def _draw_label(self, cur_drawing, label_name):
        """Draw a label for the chromosome (PRIVATE)."""
        x_position = 0.5 * (self.start_x_position + self.end_x_position)
        y_position = self.end_y_position

        label_string = String(x_position, y_position, label_name)
        label_string.fontName = "Times-BoldItalic"
        label_string.fontSize = self.title_size
        label_string.textAnchor = "middle"

        cur_drawing.add(label_string)

    def _draw_labels(self, cur_drawing, left_labels, right_labels):
        """Layout and draw sub-feature labels for the chromosome (PRIVATE).

        Tries to place each label at the same vertical position as the
        feature it applies to, but will adjust the positions to avoid or
        at least reduce label overlap.

        Draws the label text and a coloured line linking it to the
        location (i.e. feature) it applies to.
        """
        if not self._sub_components:
            return
        color_label = self._color_labels

        segment_width = (self.end_x_position - self.start_x_position) * self.chr_percent
        label_sep = (
            self.end_x_position - self.start_x_position
        ) * self.label_sep_percent
        segment_x = self.start_x_position + 0.5 * (
            self.end_x_position - self.start_x_position - segment_width
        )

        y_limits = []
        for sub_component in self._sub_components:
            y_limits.extend(
                (sub_component.start_y_position, sub_component.end_y_position)
            )
        y_min = min(y_limits)
        y_max = max(y_limits)
        del y_limits
        # Now do some label placement magic...
        # from reportlab.pdfbase import pdfmetrics
        # font = pdfmetrics.getFont('Helvetica')
        # h = (font.face.ascent + font.face.descent) * 0.90
        h = self.label_size
        for x1, x2, labels, anchor in [
            (
                segment_x,
                segment_x - label_sep,
                _place_labels(left_labels, y_min, y_max, h),
                "end",
            ),
            (
                segment_x + segment_width,
                segment_x + segment_width + label_sep,
                _place_labels(right_labels, y_min, y_max, h),
                "start",
            ),
        ]:
            for (y1, y2, color, back_color, name) in labels:
                cur_drawing.add(
                    Line(x1, y1, x2, y2, strokeColor=color, strokeWidth=0.25)
                )
                label_string = String(x2, y2, name, textAnchor=anchor)
                label_string.fontName = "Helvetica"
                label_string.fontSize = h
                if color_label:
                    label_string.fillColor = color
                if back_color:
                    w = stringWidth(name, label_string.fontName, label_string.fontSize)
                    if x1 > x2:
                        w = w * -1.0
                    cur_drawing.add(
                        Rect(
                            x2,
                            y2 - 0.1 * h,
                            w,
                            h,
                            strokeColor=back_color,
                            fillColor=back_color,
                        )
                    )
                cur_drawing.add(label_string)


class ChromosomeSegment(_ChromosomeComponent):
    """Draw a segment of a chromosome.

    This class provides the important configurable functionality of drawing
    a Chromosome. Each segment has some customization available here, or can
    be subclassed to define additional functionality. Most of the interesting
    drawing stuff is likely to happen at the ChromosomeSegment level.
    """

    def __init__(self):
        """Initialize a ChromosomeSegment.

        Attributes:
         - start_x_position, end_x_position - Defines the x range we have
           to draw things in.
         - start_y_position, end_y_position - Defines the y range we have
           to draw things in.

        Configuration Attributes:
         - scale - A scaling value for the component. By default this is
           set at 1 (ie -- has the same scale as everything else). Higher
           values give more size to the component, smaller values give less.
         - fill_color - A color to fill in the segment with. Colors are
           available in reportlab.lib.colors
         - label - A label to place on the chromosome segment. This should
           be a text string specifying what is to be included in the label.
         - label_size - The size of the label.
         - chr_percent - The percentage of area that the chromosome
           segment takes up.

        """
        _ChromosomeComponent.__init__(self)

        self.start_x_position = -1
        self.end_x_position = -1
        self.start_y_position = -1
        self.end_y_position = -1

        # --- attributes for configuration
        self.scale = 1
        self.fill_color = None
        self.label = None
        self.label_size = 6
        self.chr_percent = 0.25

    def draw(self, cur_drawing):
        """Draw a chromosome segment.

        Before drawing, the range we are drawing in needs to be set.
        """
        for position in (
            self.start_x_position,
            self.end_x_position,
            self.start_y_position,
            self.end_y_position,
        ):
            assert position != -1, "Need to set drawing coordinates."

        self._draw_subcomponents(cur_drawing)  # Anything behind
        self._draw_segment(cur_drawing)
        self._overdraw_subcomponents(cur_drawing)  # Anything on top
        self._draw_label(cur_drawing)

    def _draw_subcomponents(self, cur_drawing):
        """Draw any subcomponents of the chromosome segment (PRIVATE).

        This should be overridden in derived classes if there are
        subcomponents to be drawn.
        """
        pass

    def _draw_segment(self, cur_drawing):
        """Draw the current chromosome segment (PRIVATE)."""
        # set the coordinates of the segment -- it'll take up the MIDDLE part
        # of the space we have.
        segment_y = self.end_y_position
        segment_width = (self.end_x_position - self.start_x_position) * self.chr_percent
        segment_height = self.start_y_position - self.end_y_position
        segment_x = self.start_x_position + 0.5 * (
            self.end_x_position - self.start_x_position - segment_width
        )

        # first draw the sides of the segment
        right_line = Line(segment_x, segment_y, segment_x, segment_y + segment_height)
        left_line = Line(
            segment_x + segment_width,
            segment_y,
            segment_x + segment_width,
            segment_y + segment_height,
        )

        cur_drawing.add(right_line)
        cur_drawing.add(left_line)

        # now draw the box, if it is filled in
        if self.fill_color is not None:
            fill_rectangle = Rect(segment_x, segment_y, segment_width, segment_height)
            fill_rectangle.fillColor = self.fill_color
            fill_rectangle.strokeColor = None

            cur_drawing.add(fill_rectangle)

    def _overdraw_subcomponents(self, cur_drawing):
        """Draw any subcomponents of the chromosome segment over the main part (PRIVATE).

        This should be overridden in derived classes if there are
        subcomponents to be drawn.
        """
        pass

    def _draw_label(self, cur_drawing):
        """Add a label to the chromosome segment (PRIVATE).

        The label will be applied to the right of the segment.

        This may be overlapped by any sub-feature labels on other segments!
        """
        if self.label is not None:

            label_x = 0.5 * (self.start_x_position + self.end_x_position) + (
                self.chr_percent + 0.05
            ) * (self.end_x_position - self.start_x_position)
            label_y = (
                self.start_y_position - self.end_y_position
            ) / 2 + self.end_y_position

            label_string = String(label_x, label_y, self.label)
            label_string.fontName = "Helvetica"
            label_string.fontSize = self.label_size

            cur_drawing.add(label_string)


def _spring_layout(desired, minimum, maximum, gap=0):
    """Try to layout label co-ordinates or other floats (PRIVATE).

    Originally written for the y-axis vertical positioning of labels on a
    chromosome diagram (where the minimum gap between y-axis co-ordinates is
    the label height), it could also potentially be used for x-axis placement,
    or indeed radial placement for circular chromosomes within GenomeDiagram.

    In essence this is an optimisation problem, balancing the desire to have
    each label as close as possible to its data point, but also to spread out
    the labels to avoid overlaps. This could be described with a cost function
    (modelling the label distance from the desired placement, and the inter-
    label separations as springs) and solved as a multi-variable minimization
    problem - perhaps with NumPy or SciPy.

    For now however, the implementation is a somewhat crude ad hoc algorithm.

    NOTE - This expects the input data to have been sorted!
    """
    count = len(desired)
    if count <= 1:
        return desired  # Easy!
    if minimum >= maximum:
        raise ValueError("Bad min/max %f and %f" % (minimum, maximum))
    if min(desired) < minimum or max(desired) > maximum:
        raise ValueError(
            "Data %f to %f out of bounds (%f to %f)"
            % (min(desired), max(desired), minimum, maximum)
        )
    equal_step = float(maximum - minimum) / (count - 1)

    if equal_step < gap:
        import warnings
        from Bio import BiopythonWarning

        warnings.warn("Too many labels to avoid overlap", BiopythonWarning)
        # Crudest solution
        return [minimum + i * equal_step for i in range(count)]

    good = True
    if gap:
        prev = desired[0]
        for next in desired[1:]:
            if prev - next < gap:
                good = False
                break
    if good:
        return desired

    span = maximum - minimum
    for split in [0.5 * span, span / 3.0, 2 * span / 3.0, 0.25 * span, 0.75 * span]:
        midpoint = minimum + split
        low = [x for x in desired if x <= midpoint - 0.5 * gap]
        high = [x for x in desired if x > midpoint + 0.5 * gap]
        if len(low) + len(high) < count:
            # Bad split point, points right on boundary
            continue
        elif not low and len(high) * gap <= (span - split) + 0.5 * gap:
            # Give a little of the unused low space to the high points
            return _spring_layout(high, midpoint + 0.5 * gap, maximum, gap)
        elif not high and len(low) * gap <= split + 0.5 * gap:
            # Give a little of the unused highspace to the low points
            return _spring_layout(low, minimum, midpoint - 0.5 * gap, gap)
        elif (
            len(low) * gap <= split - 0.5 * gap
            and len(high) * gap <= (span - split) - 0.5 * gap
        ):
            return _spring_layout(
                low, minimum, midpoint - 0.5 * gap, gap
            ) + _spring_layout(high, midpoint + 0.5 * gap, maximum, gap)

    # This can be count-productive now we can split out into the telomere or
    # spacer-segment's vertical space...
    # Try not to spread out as far as the min/max unless needed
    low = min(desired)
    high = max(desired)
    if (high - low) / (count - 1) >= gap:
        # Good, we don't need the full range, and can position the
        # min and max exactly as well :)
        equal_step = (high - low) / (count - 1)
        return [low + i * equal_step for i in range(count)]

    low = 0.5 * (minimum + min(desired))
    high = 0.5 * (max(desired) + maximum)
    if (high - low) / (count - 1) >= gap:
        # Good, we don't need the full range
        equal_step = (high - low) / (count - 1)
        return [low + i * equal_step for i in range(count)]

    # Crudest solution
    return [minimum + i * equal_step for i in range(count)]


# assert False, _spring_layout([0.10,0.12,0.13,0.14,0.5,0.75, 1.0], 0, 1, 0.1)
# assert _spring_layout([0.10,0.12,0.13,0.14,0.5,0.75, 1.0], 0, 1, 0.1) == \
#     [0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0]
# assert _spring_layout([0.10,0.12,0.13,0.14,0.5,0.75, 1.0], 0, 1, 0.1) == \
#     [0.0, 0.16666666666666666, 0.33333333333333331, 0.5,
#      0.66666666666666663, 0.83333333333333326, 1.0]


def _place_labels(desired_etc, minimum, maximum, gap=0):
    # Want a list of lists/tuples for desired_etc
    desired_etc.sort()
    placed = _spring_layout([row[0] for row in desired_etc], minimum, maximum, gap)
    for old, y2 in zip(desired_etc, placed):
        # (y1, a, b, c, ..., z) --> (y1, y2, a, b, c, ..., z)
        yield (old[0], y2) + tuple(old[1:])


class AnnotatedChromosomeSegment(ChromosomeSegment):
    """Annotated chromosome segment.

    This is like the ChromosomeSegment, but accepts a list of features.
    """

    def __init__(
        self,
        bp_length,
        features,
        default_feature_color=colors.blue,
        name_qualifiers=("gene", "label", "name", "locus_tag", "product"),
    ):
        """Initialize.

        The features can either be SeqFeature objects, or tuples of values:
        start (int), end (int), strand (+1, -1, O or None), label (string),
        ReportLab color (string or object), and optional ReportLab fill color.

        Note we require 0 <= start <= end <= bp_length, and within the vertical
        space allocated to this segmenet lines will be places according to the
        start/end coordinates (starting from the top).

        Positive stand features are drawn on the right, negative on the left,
        otherwise all the way across.

        We recommend using consisent units for all the segment's scale values
        (e.g. their length in base pairs).

        When providing features as SeqFeature objects, the default color
        is used, unless the feature's qualifiers include an Artemis colour
        string (functionality also in GenomeDiagram). The caption also follows
        the GenomeDiagram approach and takes the first qualifier from the list
        or tuple specified in name_qualifiers.

        Note additional attribute label_sep_percent controls the percentage of
        area that the chromosome segment takes up, by default half of the
        chr_percent attribute (half of 25%, thus 12.5%)

        """
        ChromosomeSegment.__init__(self)
        self.bp_length = bp_length
        self.features = features
        self.default_feature_color = default_feature_color
        self.name_qualifiers = name_qualifiers
        self.label_sep_percent = self.chr_percent * 0.5

    def _overdraw_subcomponents(self, cur_drawing):
        """Draw any annotated features on the chromosome segment (PRIVATE).

        Assumes _draw_segment already called to fill out the basic shape,
        and assmes that uses the same boundaries.
        """
        # set the coordinates of the segment -- it'll take up the MIDDLE part
        # of the space we have.
        segment_y = self.end_y_position
        segment_width = (self.end_x_position - self.start_x_position) * self.chr_percent
        label_sep = (
            self.end_x_position - self.start_x_position
        ) * self.label_sep_percent
        segment_height = self.start_y_position - self.end_y_position
        segment_x = self.start_x_position + 0.5 * (
            self.end_x_position - self.start_x_position - segment_width
        )

        left_labels = []
        right_labels = []
        for f in self.features:
            try:
                # Assume SeqFeature objects
                start = f.location.start
                end = f.location.end
                strand = f.strand
                try:
                    # Handles Artemis colour integers, HTML colors, etc
                    color = _color_trans.translate(f.qualifiers["color"][0])
                except Exception:  # TODO: ValueError?
                    color = self.default_feature_color
                fill_color = color
                name = ""
                for qualifier in self.name_qualifiers:
                    if qualifier in f.qualifiers:
                        name = f.qualifiers[qualifier][0]
                        break
            except AttributeError:
                # Assume tuple of ints, string, and color
                start, end, strand, name, color = f[:5]
                color = _color_trans.translate(color)
                if len(f) > 5:
                    fill_color = _color_trans.translate(f[5])
                else:
                    fill_color = color
            assert 0 <= start <= end <= self.bp_length
            if strand == +1:
                # Right side only
                x = segment_x + segment_width * 0.6
                w = segment_width * 0.4
            elif strand == -1:
                # Left side only
                x = segment_x
                w = segment_width * 0.4
            else:
                # Both or neither - full width
                x = segment_x
                w = segment_width
            local_scale = segment_height / self.bp_length
            fill_rectangle = Rect(
                x,
                segment_y + segment_height - local_scale * start,
                w,
                local_scale * (start - end),
            )
            fill_rectangle.fillColor = fill_color
            fill_rectangle.strokeColor = color
            cur_drawing.add(fill_rectangle)
            if name:
                if fill_color == color:
                    back_color = None
                else:
                    back_color = fill_color
                value = (
                    segment_y + segment_height - local_scale * start,
                    color,
                    back_color,
                    name,
                )
                if strand == -1:
                    self._left_labels.append(value)
                else:
                    self._right_labels.append(value)


class TelomereSegment(ChromosomeSegment):
    """A segment that is located at the end of a linear chromosome.

    This is just like a regular segment, but it draws the end of a chromosome
    which is represented by a half circle. This just overrides the
    _draw_segment class of ChromosomeSegment to provide that specialized
    drawing.
    """

    def __init__(self, inverted=0):
        """Initialize a segment at the end of a chromosome.

        See ChromosomeSegment for all of the attributes that can be
        customized in a TelomereSegments.

        Arguments:
         - inverted -- Whether or not the telomere should be inverted
           (ie. drawn on the bottom of a chromosome)

        """
        ChromosomeSegment.__init__(self)

        self._inverted = inverted

    def _draw_segment(self, cur_drawing):
        """Draw a half circle representing the end of a linear chromosome (PRIVATE)."""
        # set the coordinates of the segment -- it'll take up the MIDDLE part
        # of the space we have.
        width = (self.end_x_position - self.start_x_position) * self.chr_percent
        height = self.start_y_position - self.end_y_position
        center_x = 0.5 * (self.end_x_position + self.start_x_position)
        start_x = center_x - 0.5 * width
        if self._inverted:
            center_y = self.start_y_position
            start_angle = 180
            end_angle = 360
        else:
            center_y = self.end_y_position
            start_angle = 0
            end_angle = 180

        cap_wedge = Wedge(center_x, center_y, width / 2, start_angle, end_angle, height)
        cap_wedge.strokeColor = None
        cap_wedge.fillColor = self.fill_color
        cur_drawing.add(cap_wedge)

        # Now draw an arc for the curved edge of the wedge,
        # omitting the flat end.
        cap_arc = ArcPath()
        cap_arc.addArc(center_x, center_y, width / 2, start_angle, end_angle, height)
        cur_drawing.add(cap_arc)


class SpacerSegment(ChromosomeSegment):
    """A segment that is located at the end of a linear chromosome.

    Doesn't draw anything, just empty space which can be helpful
    for layout purposes (e.g. making room for feature labels).
    """

    def draw(self, cur_diagram):
        """Draw nothing to the current diagram (dummy method).

        The segment spacer has no actual image in the diagram,
        so this method therefore does nothing, but is defined
        to match the expected API of the other segment objects.
        """
        pass
