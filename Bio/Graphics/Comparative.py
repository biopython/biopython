# Copyright 2001 by Brad Chapman.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Plots to compare information between different sources.

This file contains high level plots which are designed to be used to
compare different types of information. The most basic example is comparing
two variables in a traditional scatter plot.
"""
# reportlab
from reportlab.lib import colors
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

from reportlab.graphics.shapes import Drawing, String
from reportlab.graphics.charts.markers import makeEmptySquare, makeFilledSquare
from reportlab.graphics.charts.markers import makeFilledDiamond, makeSmiley
from reportlab.graphics.charts.markers import makeFilledCircle, makeEmptyCircle

from Bio.Graphics import _write


class ComparativeScatterPlot:
    """Display a scatter-type plot comparing two different kinds of info.

    Attributes;
     - display_info - a 2D list of the information we'll be outputting. Each
       top level list is a different data type, and each data point is a
       two-tuple of the coordinates of a point.

    So if you had two distributions of points, it should look like::

       display_info = [[(1, 2), (3, 4)],
                       [(5, 6), (7, 8)]]

    If everything is just one set of points, display_info can look like::

        display_info = [[(1, 2), (3, 4), (5, 6)]]

    """

    def __init__(self, output_format="pdf"):
        """Initialize the class."""
        # customizable attributes
        self.number_of_columns = 1
        self.page_size = letter
        self.title_size = 20

        self.output_format = output_format

        # the information we'll be writing
        self.display_info = []

        # initial colors and shapes used for drawing points
        self.color_choices = [
            colors.red,
            colors.green,
            colors.blue,
            colors.yellow,
            colors.orange,
            colors.black,
        ]
        self.shape_choices = [
            makeFilledCircle,
            makeEmptySquare,
            makeFilledDiamond,
            makeFilledSquare,
            makeEmptyCircle,
            makeSmiley,
        ]

    def draw_to_file(self, output_file, title):
        """Write the comparative plot to a file.

        Arguments:
         - output_file - The name of the file to output the information to,
           or a handle to write to.
         - title - A title to display on the graphic.

        """
        width, height = self.page_size
        cur_drawing = Drawing(width, height)

        self._draw_title(cur_drawing, title, width, height)

        start_x = inch * 0.5
        end_x = width - inch * 0.5
        end_y = height - 1.5 * inch
        start_y = 0.5 * inch
        self._draw_scatter_plot(cur_drawing, start_x, start_y, end_x, end_y)

        return _write(cur_drawing, output_file, self.output_format)

    def _draw_title(self, cur_drawing, title, width, height):
        """Add a title to the page we are outputting (PRIVATE)."""
        title_string = String(width / 2, height - inch, title)
        title_string.fontName = "Helvetica-Bold"
        title_string.fontSize = self.title_size
        title_string.textAnchor = "middle"

        cur_drawing.add(title_string)

    def _draw_scatter_plot(self, cur_drawing, x_start, y_start, x_end, y_end):
        """Draw a scatter plot on the drawing with the given coordinates (PRIVATE)."""
        scatter_plot = LinePlot()

        # set the dimensions of the scatter plot
        scatter_plot.x = x_start
        scatter_plot.y = y_start
        scatter_plot.width = abs(x_start - x_end)
        scatter_plot.height = abs(y_start - y_end)

        scatter_plot.data = self.display_info

        scatter_plot.joinedLines = 0

        # set the axes of the plot
        x_min, x_max, y_min, y_max = self._find_min_max(self.display_info)
        scatter_plot.xValueAxis.valueMin = x_min
        scatter_plot.xValueAxis.valueMax = x_max
        scatter_plot.xValueAxis.valueStep = (x_max - x_min) / 10.0

        scatter_plot.yValueAxis.valueMin = y_min
        scatter_plot.yValueAxis.valueMax = y_max
        scatter_plot.yValueAxis.valueStep = (y_max - y_min) / 10.0

        self._set_colors_and_shapes(scatter_plot, self.display_info)

        cur_drawing.add(scatter_plot)

    def _set_colors_and_shapes(self, scatter_plot, display_info):
        """Set the colors and shapes of the points displayed (PRIVATE).

        By default this just sets all of the points according to the order
        of colors and shapes defined in self.color_choices and
        self.shape_choices. The first 5 shapes and colors are unique, the
        rest of them are just set to the same color and shape (since I
        ran out of shapes!).

        You can change how this function works by either changing the
        values of the color_choices and shape_choices attributes, or
        by inheriting from this class and overriding this function.
        """
        for value_num in range(len(display_info)):
            # if we have unique colors, add them
            if (value_num + 1) < len(self.color_choices):
                scatter_plot.lines[value_num].strokeColor = self.color_choices[
                    value_num
                ]
                scatter_plot.lines[value_num].symbol = self.shape_choices[value_num]
            # otherwise just use the last number
            else:
                scatter_plot.lines[value_num].strokeColor = self.color_choices[-1]
                scatter_plot.lines[value_num].symbol = self.shape_choices[-1]

    def _find_min_max(self, info):
        """Find min and max for x and y coordinates in the given data (PRIVATE)."""
        x_min = info[0][0][0]
        x_max = info[0][0][0]
        y_min = info[0][0][1]
        y_max = info[0][0][1]

        for two_d_list in info:
            for x, y in two_d_list:
                if x > x_max:
                    x_max = x
                if x < x_min:
                    x_min = x
                if y > y_max:
                    y_max = y
                if y < y_min:
                    y_min = y

        return x_min, x_max, y_min, y_max
