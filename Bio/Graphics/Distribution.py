# Copyright 2001 by Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Display information distributed across a Chromosome-like object.

These classes are meant to show the distribution of some kind of information
as it changes across any kind of segment. It was designed with chromosome
distributions in mind, but could also work for chromosome regions, BAC clones
or anything similar.

Reportlab is used for producing the graphical output.
"""
# standard library
import math

# reportlab
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors

from reportlab.graphics.shapes import Drawing, String
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.barcharts import BarChartProperties
from reportlab.graphics.widgetbase import TypedPropertyCollection

from Bio.Graphics import _write


class DistributionPage:
    """Display a grouping of distributions on a page.

    This organizes Distributions, and will display them nicely
    on a single page.
    """

    def __init__(self, output_format="pdf"):
        """Initialize."""
        self.distributions = []

        # customizable attributes
        self.number_of_columns = 1
        self.page_size = letter
        self.title_size = 20

        self.output_format = output_format

    def draw(self, output_file, title):
        """Draw out the distribution information.

        Arguments:
         - output_file - The name of the file to output the information to,
           or a handle to write to.
         - title - A title to display on the graphic.

        """
        width, height = self.page_size
        cur_drawing = Drawing(width, height)

        self._draw_title(cur_drawing, title, width, height)

        # calculate the x and y position changes for each distribution
        cur_x_pos = inch * 0.5
        end_x_pos = width - inch * 0.5
        cur_y_pos = height - 1.5 * inch
        end_y_pos = 0.5 * inch
        x_pos_change = (end_x_pos - cur_x_pos) / float(self.number_of_columns)
        num_y_rows = math.ceil(
            float(len(self.distributions)) / float(self.number_of_columns)
        )
        y_pos_change = (cur_y_pos - end_y_pos) / num_y_rows

        self._draw_distributions(
            cur_drawing, cur_x_pos, x_pos_change, cur_y_pos, y_pos_change, num_y_rows
        )
        self._draw_legend(cur_drawing, 2.5 * inch, width)

        return _write(cur_drawing, output_file, self.output_format)

    def _draw_title(self, cur_drawing, title, width, height):
        """Add the title of the figure to the drawing (PRIVATE)."""
        title_string = String(width / 2, height - inch, title)
        title_string.fontName = "Helvetica-Bold"
        title_string.fontSize = self.title_size
        title_string.textAnchor = "middle"

        cur_drawing.add(title_string)

    def _draw_distributions(
        self,
        cur_drawing,
        start_x_pos,
        x_pos_change,
        start_y_pos,
        y_pos_change,
        num_y_drawings,
    ):
        """Draw all of the distributions on the page (PRIVATE).

        Arguments:
         - cur_drawing - The drawing we are working with.
         - start_x_pos - The x position on the page to start drawing at.
         - x_pos_change - The change in x position between each figure.
         - start_y_pos - The y position on the page to start drawing at.
         - y_pos_change - The change in y position between each figure.
         - num_y_drawings - The number of drawings we'll have in the y
           (up/down) direction.

        """
        for y_drawing in range(int(num_y_drawings)):
            # if we are on the last y position, we may not be able
            # to fill all of the x columns
            if (y_drawing + 1) * self.number_of_columns > len(self.distributions):
                num_x_drawings = (
                    len(self.distributions) - y_drawing * self.number_of_columns
                )
            else:
                num_x_drawings = self.number_of_columns
            for x_drawing in range(num_x_drawings):
                dist_num = y_drawing * self.number_of_columns + x_drawing
                cur_distribution = self.distributions[dist_num]

                # find the x and y boundaries of the distribution
                x_pos = start_x_pos + x_drawing * x_pos_change
                end_x_pos = x_pos + x_pos_change
                end_y_pos = start_y_pos - y_drawing * y_pos_change
                y_pos = end_y_pos - y_pos_change

                # draw the distribution
                cur_distribution.draw(cur_drawing, x_pos, y_pos, end_x_pos, end_y_pos)

    def _draw_legend(self, cur_drawing, start_y, width):
        """Add a legend to the figure (PRIVATE).

        Subclasses can implement to provide a specialized legend.
        """
        pass


class BarChartDistribution:
    """Display the distribution of values as a bunch of bars."""

    def __init__(self, display_info=None):
        """Initialize a Bar Chart display of distribution info.

        Attributes:
         - display_info - the information to be displayed in the distribution.
           This should be ordered as a list of lists, where each internal list
           is a data set to display in the bar chart.

        """
        if display_info is None:
            display_info = []
        self.display_info = display_info

        self.x_axis_title = ""
        self.y_axis_title = ""
        self.chart_title = ""
        self.chart_title_size = 10

        self.padding_percent = 0.15

    def draw(self, cur_drawing, start_x, start_y, end_x, end_y):
        """Draw a bar chart with the info in the specified range."""
        bar_chart = VerticalBarChart()
        if self.chart_title:
            self._draw_title(
                cur_drawing, self.chart_title, start_x, start_y, end_x, end_y
            )
        # set the position of the bar chart
        x_start, x_end, y_start, y_end = self._determine_position(
            start_x, start_y, end_x, end_y
        )

        bar_chart.x = x_start
        bar_chart.y = y_start
        bar_chart.width = abs(x_start - x_end)
        bar_chart.height = abs(y_start - y_end)

        # set the information in the bar chart
        bar_chart.data = self.display_info
        bar_chart.valueAxis.valueMin = min(self.display_info[0])
        bar_chart.valueAxis.valueMax = max(self.display_info[0])
        for data_set in self.display_info[1:]:
            if min(data_set) < bar_chart.valueAxis.valueMin:
                bar_chart.valueAxis.valueMin = min(data_set)
            if max(data_set) > bar_chart.valueAxis.valueMax:
                bar_chart.valueAxis.valueMax = max(data_set)

        # set other formatting options
        if len(self.display_info) == 1:
            bar_chart.groupSpacing = 0
            style = TypedPropertyCollection(BarChartProperties)
            style.strokeWidth = 0
            style.strokeColor = colors.green
            style[0].fillColor = colors.green

            bar_chart.bars = style

        # set the labels
        # XXX labels don't work yet
        # bar_chart.valueAxis.title = self.x_axis_title
        # bar_chart.categoryAxis.title = self.y_axis_title

        cur_drawing.add(bar_chart)

    def _draw_title(self, cur_drawing, title, start_x, start_y, end_x, end_y):
        """Add the title of the figure to the drawing (PRIVATE)."""
        x_center = start_x + (end_x - start_x) / 2
        y_pos = end_y + (self.padding_percent * (start_y - end_y)) / 2
        title_string = String(x_center, y_pos, title)
        title_string.fontName = "Helvetica-Bold"
        title_string.fontSize = self.chart_title_size
        title_string.textAnchor = "middle"

        cur_drawing.add(title_string)

    def _determine_position(self, start_x, start_y, end_x, end_y):
        """Calculate the position of the chart with blank space (PRIVATE).

        This uses some padding around the chart, and takes into account
        whether the chart has a title. It returns 4 values, which are,
        in order, the x_start, x_end, y_start and y_end of the chart
        itself.
        """
        x_padding = self.padding_percent * (end_x - start_x)
        y_padding = self.padding_percent * (start_y - end_y)

        new_x_start = start_x + x_padding
        new_x_end = end_x - x_padding

        if self.chart_title:
            new_y_start = start_y - y_padding - self.chart_title_size
        else:
            new_y_start = start_y - y_padding

        new_y_end = end_y + y_padding

        return new_x_start, new_x_end, new_y_start, new_y_end


class LineDistribution:
    """Display the distribution of values as connected lines.

    This distribution displays the change in values across the object as
    lines. This also allows multiple distributions to be displayed on a
    single graph.
    """

    def __init__(self):
        """Initialize."""
        pass

    def draw(self, cur_drawing, start_x, start_y, end_x, end_y):
        """Draw a line distribution into the current drawing."""
        pass
