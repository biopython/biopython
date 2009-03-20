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
# standard library
import os

# reportlab
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors

from reportlab.graphics.shapes import Drawing, String, Line, Rect, Wedge
from reportlab.graphics import renderPDF, renderPS
from reportlab.graphics.widgetbase import Widget

from Bio.Graphics import _write

class _ChromosomeComponent(Widget):
    """Base class specifying the interface for a component of the system.

    This class should not be instantiated directly, but should be used
    from derived classes.
    """
    def __init__(self):
        """Initialize a chromosome component.

        Attributes:

        o _sub_components -- Any components which are contained under
        this parent component. This attribute should be accessed through
        the add() and remove() functions.
        """
        self._sub_components = []

    def add(self, component):
        """Add a sub_component to the list of components under this item.
        """
        assert isinstance(component, _ChromosomeComponent), \
               "Expected a _ChromosomeComponent object, got %s" % component
        
        self._sub_components.append(component)

    def remove(self, component):
        """Remove the specified component from the subcomponents.

        Raises a ValueError if the component is not registered as a
        sub_component.
        """
        try:
            self._sub_components.remove(component)
        except ValueError:
            raise ValueError("Component %s not found in sub_components." %
                             component)

    def draw(self):
        """Draw the specified component.
        """
        raise AssertionError("Subclasses must implement.")
    
class Organism(_ChromosomeComponent):
    """Top level class for drawing chromosomes.

    This class holds information about an organism and all of it's
    chromosomes, and provides the top level object which could be used
    for drawing a chromosome representation of an organism.

    Chromosomes should be added and removed from the Organism via the
    add and remove functions.
    """
    def __init__(self, output_format = 'pdf'):
        _ChromosomeComponent.__init__(self)

        # customizable attributes
        self.page_size = letter
        self.title_size = 20

        self.output_format = output_format

    def draw(self, output_file, title):
        """Draw out the information for the Organism.

        Arguments:

        o output_file -- The name of a file specifying where the
        document should be saved, or a handle to be written to.
        The output format is set when creating the Organism object.

        o title -- The output title of the produced document.
        """
        width, height = self.page_size
        cur_drawing = Drawing(width, height)

        self._draw_title(cur_drawing, title, width, height)

        cur_x_pos = inch * .5
        if len(self._sub_components) > 0:
            x_pos_change = (width - inch) / len(self._sub_components)
        # no sub_components
        else:
            pass
        
        for sub_component in self._sub_components:
            # set the drawing location of the chromosome
            sub_component.start_x_position = cur_x_pos
            sub_component.end_x_position = cur_x_pos + .9 * x_pos_change
            sub_component.start_y_position = height - 1.5 * inch
            sub_component.end_y_position = 3 * inch

            # do the drawing
            sub_component.draw(cur_drawing)

            # update the locations for the next chromosome
            cur_x_pos += x_pos_change

        self._draw_legend(cur_drawing, 2.5 * inch, width)

        return _write(cur_drawing, output_file, self.output_format)

    def _draw_title(self, cur_drawing, title, width, height):
        """Write out the title of the organism figure.
        """
        title_string = String(width / 2, height - inch, title)
        title_string.fontName = 'Helvetica-Bold'
        title_string.fontSize = self.title_size
        title_string.textAnchor = "middle"

        cur_drawing.add(title_string)

    def _draw_legend(self, cur_drawing, start_y, width):
        """Draw a legend for the figure.

        Subclasses should implement this to provide specialized legends.
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

        o chromosome_name - The label for the chromosome.

        Attributes:

        o start_x_position, end_x_position - The x positions on the page
        where the chromosome should be drawn. This allows multiple
        chromosomes to be drawn on a single page.

        o start_y_position, end_y_position - The y positions on the page
        where the chromosome should be contained.

        Configuration Attributes:

        o title_size - The size of the chromosome title.

        o scale_num - A number of scale the drawing by. This is useful if
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

    def subcomponent_size(self):
        """Return the scaled size of all subcomponents of this component.
        """
        total_sub = 0
        for sub_component in self._sub_components:
            total_sub += sub_component.scale

        return total_sub

    def draw(self, cur_drawing):
        """Draw a chromosome on the specified template.

        Ideally, the x_position and y_*_position attributes should be
        set prior to drawing -- otherwise we're going to have some problems.
        """
        for position in (self.start_x_position, self.end_x_position,
                         self.start_y_position, self.end_y_position):
            assert position != -1, "Need to set drawing coordinates."

        # first draw all of the sub-sections of the chromosome -- this
        # will actually be the picture of the chromosome
        cur_y_pos = self.start_y_position
        if self.scale_num:
            y_pos_change = ((self.start_y_position * .95 - self.end_y_position)
                            / self.scale_num)
        elif len(self._sub_components) > 0:
            y_pos_change = ((self.start_y_position * .95 - self.end_y_position)
                            / self.subcomponent_size())
        # no sub_components to draw
        else:
            pass
            
        for sub_component in self._sub_components:
            this_y_pos_change = sub_component.scale * y_pos_change
            
            # set the location of the component to draw
            sub_component.start_x_position = self.start_x_position
            sub_component.end_x_position = self.end_x_position
            sub_component.start_y_position = cur_y_pos
            sub_component.end_y_position = cur_y_pos - this_y_pos_change

            # draw the sub component
            sub_component.draw(cur_drawing)

            # update the position for the next component
            cur_y_pos -= this_y_pos_change

        self._draw_label(cur_drawing, self._name)

    def _draw_label(self, cur_drawing, label_name):
        """Draw a label for the chromosome.
        """
        x_position = self.start_x_position
        y_position = self.end_y_position

        label_string = String(x_position, y_position, label_name)
        label_string.fontName = 'Times-BoldItalic'
        label_string.fontSize = self.title_size
        label_string.textAnchor = 'start'

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
        o start_x_position, end_x_position - Defines the x range we have
        to draw things in.

        o start_y_position, end_y_position - Defines the y range we have
        to draw things in.

        Configuration Attributes:

        o scale - A scaling value for the component. By default this is
        set at 1 (ie -- has the same scale as everything else). Higher
        values give more size to the component, smaller values give less.

        o fill_color - A color to fill in the segment with. Colors are
        available in reportlab.lib.colors

        o label - A label to place on the chromosome segment. This should
        be a text string specifying what is to be included in the label.

        o label_size - The size of the label.

        o chr_percent - The percentage of area that the chromosome
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
        self.chr_percent = .25

    def draw(self, cur_drawing):
        """Draw a chromosome segment.

        Before drawing, the range we are drawing in needs to be set.
        """
        for position in (self.start_x_position, self.end_x_position,
                         self.start_y_position, self.end_y_position):
            assert position != -1, "Need to set drawing coordinates."

        self._draw_subcomponents(cur_drawing)
        self._draw_segment(cur_drawing)
        self._draw_label(cur_drawing)

    def _draw_subcomponents(self, cur_drawing):
        """Draw any subcomponents of the chromosome segment.

        This should be overridden in derived classes if there are
        subcomponents to be drawn.
        """
        pass

    def _draw_segment(self, cur_drawing):
        """Draw the current chromosome segment.
        """
        # set the coordinates of the segment -- it'll take up the left part
        # of the space we have.
        segment_x = self.start_x_position
        segment_y = self.end_y_position
        segment_width = (self.end_x_position - self.start_x_position) \
                        * self.chr_percent
        segment_height = self.start_y_position - self.end_y_position
        
        # first draw the sides of the segment
        right_line = Line(segment_x, segment_y,
                          segment_x, segment_y + segment_height)
        left_line = Line(segment_x + segment_width, segment_y,
                         segment_x + segment_width, segment_y + segment_height)
        
        cur_drawing.add(right_line)
        cur_drawing.add(left_line)
        
        # now draw the box, if it is filled in
        if self.fill_color is not None:
            fill_rectangle = Rect(segment_x, segment_y,
                                  segment_width, segment_height)
            fill_rectangle.fillColor = self.fill_color
            fill_rectangle.strokeColor = None

            cur_drawing.add(fill_rectangle)
            
    def _draw_label(self, cur_drawing):
        """Add a label to the chromosome segment.
        """
        # the label will be applied to the right of the segment
        if self.label is not None:

            label_x = self.start_x_position + \
                      (self.chr_percent + 0.05) * (self.end_x_position -
                                                   self.start_x_position)
            label_y = ((self.start_y_position - self.end_y_position) / 2 +
                       self.end_y_position)

            label_string = String(label_x, label_y, self.label)
            label_string.fontName = 'Helvetica'
            label_string.fontSize = self.label_size

            cur_drawing.add(label_string)
        
class TelomereSegment(ChromosomeSegment):
    """A segment that is located at the end of a linear chromosome.

    This is just like a regular segment, but it draws the end of a chromosome
    which is represented by a half circle. This just overrides the
    _draw_segment class of ChromosomeSegment to provide that specialized
    drawing.
    """
    def __init__(self, inverted = 0):
        """Initialize a segment at the end of a chromosome.

        See ChromosomeSegment for all of the attributes that can be
        customized in a TelomereSegments.

        Arguments:

        o inverted -- Whether or not the telomere should be inverted
        (ie. drawn on the bottom of a chromosome)
        """
        ChromosomeSegment.__init__(self)

        self._inverted = inverted

    def _draw_segment(self, cur_drawing):
        """Draw a half circle representing the end of a linear chromosome.
        """
        # set the coordinates of the segment -- it'll take up the left part
        # of the space we have.
        width = (self.end_x_position - self.start_x_position) \
                * self.chr_percent
        height = self.start_y_position - self.end_y_position

        center_x = self.start_x_position + width / 2
        if self._inverted:
            center_y = self.start_y_position
            start_angle = 180
            end_angle = 360
        else:
            center_y = self.end_y_position
            start_angle = 0
            end_angle = 180
        
        cap_wedge = Wedge(center_x, center_y, width / 2,
                          start_angle, end_angle, height / 2)

        cap_wedge.fillColor = self.fill_color
        cur_drawing.add(cap_wedge)

        # draw a line to cover up the the bottom part of the wedge
        if self._inverted:
            cover_line = Line(self.start_x_position, self.start_y_position,
                              self.start_x_position + width,
                              self.start_y_position)
        else:
            cover_line = Line(self.start_x_position, self.end_y_position,
                              self.start_x_position + width,
                              self.end_y_position)

        if self.fill_color is not None:
            cover_color = self.fill_color
        else:
            cover_color = colors.white

        cover_line.strokeColor = cover_color
        cur_drawing.add(cover_line)


        
        

        
        

        
                   
