#!/usr/bin/env python
# Copyright 2013 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for general functionality of the ColorSpiral utility."""

# Builtins
import colorsys
from math import pi
import os
import unittest
from cmath import rect

# Do we have ReportLab?  Raise error if not present.
from Bio import MissingPythonDependencyError
try:
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.lib.pagesizes import A4
except ImportError:
    raise MissingPythonDependencyError(
            "Install reportlab if you want to use Bio.Graphics.")


# Biopython Bio.Graphics.ColorSpiral
from Bio.Graphics.ColorSpiral import ColorSpiral, get_colors, get_color_dict


class SpiralTest(unittest.TestCase):
    """Construct and draw ColorSpiral colours placed on HSV spiral."""
    def setUp(self):
        """Set up canvas for drawing"""
        output_filename = os.path.join("Graphics", "spiral_test.pdf")
        self.c = Canvas(output_filename, pagesize=A4)
        # co-ordinates of the centre of the canvas
        self.x_0, self.y_0 = 0.5 * A4[0], 0.5 * A4[1]

    def test_colorlist(self):
        """Get set of eight colours, no jitter, using ColorSpiral."""
        cs = ColorSpiral(a=4, b=0.33, jitter=0)
        colours = list(cs.get_colors(8))
        cstr = ["(%.2f, %.2f, %.2f)" % (r, g, b)
                for r, g, b in colours]
        expected = \
            ['(0.64, 0.74, 0.81)', '(0.68, 0.52, 0.76)', '(0.72, 0.41, 0.55)',
             '(0.68, 0.39, 0.31)', '(0.63, 0.54, 0.22)', '(0.48, 0.59, 0.13)',
             '(0.24, 0.54, 0.06)', '(0.01, 0.50, -0.00)']
        self.assertEqual(cstr, expected)

    def test_colorspiral(self):
        """Get set of 16 colours, no jitter, using ColorSpiral."""
        cs = ColorSpiral(a=4, b=0.33, jitter=0)
        radius = A4[0] * 0.025
        for r, g, b in cs.get_colors(16):
            self.c.setFillColor((r, g, b))
            # Convert HSV colour to rectangular coordinates on HSV disc
            h, s, v = colorsys.rgb_to_hsv(r, g, b)
            coords = rect(s * A4[0] * 0.45, h * 2 * pi)
            x, y = self.x_0 + coords.real, self.y_0 + coords.imag
            self.c.ellipse(x - radius, y - radius, x + radius, y + radius,
                           stroke=0, fill=1)
        self.finish()

    def finish(self):
        """Clean up and save image."""
        self.c.save()


class SquareTest(unittest.TestCase):
    """Construct and draw ColorSpiral colours placed in a square, with jitter."""
    def setUp(self):
        """Set up canvas for drawing"""
        output_filename = os.path.join("Graphics", "square_test.pdf")
        self.c = Canvas(output_filename, pagesize=(500, 500))

    def test_colorspiral(self):
        """Set of 625 colours, with jitter, using get_colors()."""
        boxedge = 20
        boxes_per_row = 25
        rows = 0
        for i, c in enumerate(get_colors(625)):
            self.c.setFillColor(c)
            x1 = boxedge * (i % boxes_per_row)
            y1 = rows * boxedge
            self.c.rect(x1, y1, boxedge, boxedge, fill=1, stroke=0)
            if not (i + 1) % boxes_per_row:
                rows += 1
        self.finish()

    def finish(self):
        """Clean up and save image."""
        self.c.save()


class DictTest(unittest.TestCase):
    """Generate set of colours on the basis of an iterable."""
    def test_dict(self):
        """get_color_dict() for classes A-D, no jitter."""
        classes = ['A', 'B', 'C', 'D']
        colors = get_color_dict(classes, jitter=0)
        cstr = ["%s: (%.2f, %.2f, %.2f)" % (c, r, g, b)
                for c, (r, g, b) in sorted(colors.items())]
        expected = ['A: (0.52, 0.76, 0.69)',
                    'B: (0.40, 0.31, 0.68)',
                    'C: (0.59, 0.13, 0.47)',
                    'D: (0.50, 0.00, 0.00)']
        self.assertEqual(cstr, expected)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
