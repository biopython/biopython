# Copyright 2008 by Brad Chapman. All rights reserved.
# Copyright 2008 by Michiel de Hoon.  All rights reserved.
# Copyright 2009-2017 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Graphics offers several graphical outputs, all using ReportLab."""

# Check if ReportLab is installed.
try:
    import reportlab as r
    del r
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Please install ReportLab if you want "
        "to use Bio.Graphics. You can find ReportLab at "
        "http://www.reportlab.com/software/opensource/")


# The following code is to allow all the Bio.Graphics
# code to deal with the different ReportLab renderers
# and the API quirks consistently.


def _write(drawing, output_file, format, dpi=72):
    """Standardize output to files (PRIVATE).

    Writes the provided drawing out to a file in a prescribed format.

      - drawing - suitable ReportLab drawing object.
      - output_file - a handle to write to, or a filename to write to.
      - format - String indicating output format, one of PS, PDF, SVG,
        or provided the ReportLab renderPM module is installed,
        one of the bitmap formats JPG, BMP, GIF, PNG, TIFF or TIFF.
        The format can be given in any case.
      - dpi - Resolution (dots per inch) for bitmap formats.

    No return value.
    """
    from reportlab.graphics import renderPS, renderPDF, renderSVG
    try:
        from reportlab.graphics import renderPM
    except ImportError:
        # This is an optional part of ReportLab, so may not be installed.
        # We'll raise a missing dependency error if rendering to a
        # bitmap format is attempted.
        renderPM = None

    formatdict = {'PS': renderPS, 'EPS': renderPS,
                  # not sure which you actually get, PS or EPS, but
                  # GenomeDiagram used PS while other modules used EPS.
                  'PDF': renderPDF,
                  'SVG': renderSVG,
                  'JPG': renderPM,
                  'BMP': renderPM,
                  'GIF': renderPM,
                  'PNG': renderPM,
                  'TIFF': renderPM,
                  'TIF': renderPM
                  }
    try:
        # If output is not a string, then .upper() will trigger
        # an attribute error...
        drawmethod = formatdict[format.upper()]  # select drawing method
    except (KeyError, AttributeError):
        raise ValueError("Output format should be one of %s"
                         % ", ".join(formatdict))

    if drawmethod is None:
        # i.e. We wanted renderPM but it isn't installed
        # See the import at the top of the function.
        from Bio import MissingPythonDependencyError
        raise MissingPythonDependencyError(
            "Please install ReportLab's renderPM module")

    if drawmethod == renderPM:
        # This has a different API to the other render objects
        return drawmethod.drawToFile(drawing, output_file,
                                     format, dpi=dpi)
    else:
        return drawmethod.drawToFile(drawing, output_file)
