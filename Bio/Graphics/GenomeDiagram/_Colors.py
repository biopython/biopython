# Copyright 2003-2008 by Leighton Pritchard.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Contact:       Leighton Pritchard, Scottish Crop Research Institute,
#                Invergowrie, Dundee, Scotland, DD2 5DA, UK
#                L.Pritchard@scri.ac.uk
################################################################################

""" Colors module

    Provides:

    o ColorTranslator -  class to convert tuples of integers and floats into
                            colors.Color objects

    For drawing capabilities, this module uses reportlab to define colors:

    http://www.reportlab.com
"""

# ReportLab imports
from reportlab.lib import colors

class ColorTranslator(object):
    """ Class providing methods for translating representations of color into
    """
    def __init__(self, filename=None):
        """ __init__(self, filename)

            o filename      Location of a file containing colorscheme
                            information

            Optional parameters set the color scheme 
        """
        self._artemis_colorscheme = {0: (colors.Color(1, 1, 1,), "pathogenicity, adaptation, chaperones"),
               1: (colors.Color(0.39, 0.39, 0.39), "energy metabolism"),
               2: (colors.Color(1, 0, 0), "information transfer"),
               3: (colors.Color(0, 1, 0), "surface"),
               4: (colors.Color(0, 0, 1), "stable RNA"),
               5: (colors.Color(0, 1, 1), "degradation of large molecules"),
               6: (colors.Color(1, 0, 1), "degradation of small molecules"),
               7: (colors.Color(1, 1, 0), "central/intermediary/miscellaneous metabolism"),
               8: (colors.Color(0.60, 0.98, 0.60), "unknown"),
               9: (colors.Color(0.53, 0.81, 0.98), "regulators"),
               10: (colors.Color(1, 0.65, 0), "conserved hypotheticals"),
               11: (colors.Color(0.78, 0.59, 0.39), "pseudogenes and partial genes"),
               12: (colors.Color(1, 0.78, 0.78), "phage/IS elements"),
               13: (colors.Color(0.70, 0.70, 0.70), "some miscellaneous information"),
               14: (colors.Color(0, 0, 0), ""),
               15: (colors.Color(1, 0.25, 0.25), "secondary metabolism"),
               16: (colors.Color(1, 0.5, 0.5), ""),
               17: (colors.Color(1, 0.75, 0.75), "")
                }      # Hardwired Artemis color scheme                                                                        
        self._colorscheme = {}
        if filename is not None:
            self.read_colorscheme(filename)# Imported color scheme
        else:
            self._colorscheme = self._artemis_colorscheme


    def translate(self, color=None, colour=None):
        """ translate(self, color)

            o color    Color defined as an int, a tuple of three ints 0->255
                       or a tuple of three floats 0 -> 1, or a string giving
                       one of the named colors defined by ReportLab, or a
                       ReportLab color object (returned as is).

                       (This argument is overridden by a backwards compatible
                       argument with UK spelling, colour).

            Returns a colors.Color object, determined semi-intelligently
            depending on the input values
        """
        #Let the UK spelling (colour) override the USA spelling (color)
        if colour is not None:
            color = colour
        
        if color is None:
            raise ValueError, "Passed color (or colour) must be a valid color type"
        elif isinstance(color, int):
            color = self.scheme_color(color)
        elif isinstance(color, colors.Color):
            return color
        elif isinstance(color, basestring):
            #Assume its a named reportlab color like "red".
            color = colors.toColor(color)
        elif type(color) == type((1., 2., 3.)) and type(color[0]) == type(1.):
            color = self.float1_color(color)        
        elif type(color) == type((1, 2, 3)) and type(color[0]) == type(1):
            color = self.int255_color(color)
        return color
        

    def read_colorscheme(self, filename):
        """ read_colorscheme(self, filename)

            o filename      The location of a file defining colors in tab-separated
                            format plaintext as:
                            INT \t RED \t GREEN \t BLUE \t Comment
                            Where RED, GREEN and BLUE are intensities in the range
                            0 -> 255
                            e.g.
                            2 \t 255 \t 0 \t 0 \t Red: Information transfer

            Reads information from a file containing color information and
            stores it internally
        """
        lines = open(filename, 'r').readlines()
        for line in lines:
            data = line.strip().split('\t')
            try:
                label = int(data[0])
                red, green, blue = int(data[1]), int(data[2]), int(data[3])
                if len(data) > 4:
                    comment = data[4]
                else:
                    comment = ""
                self._colorscheme[label] = (self.int255_color((red, green, blue)),
                                             comment)
            except:
                raise IOError, "Expected INT \t INT \t INT \t INT \t string input"

    def get_artemis_colorscheme(self):
        """ get_artemis_colorscheme(self)

            Return the Artemis color scheme as a dictionary
        """
        return self._artemis_colorscheme
        

    def artemis_color(self, value):
        """ artemis_color(self, value)

            o value     An int representing a functional class in the Artemis
                        color scheme (see www.sanger.ac.uk for a description),
                        or a string from a GenBank feature annotation for the
                        color which may be dot delimited (in which case the
                        first value is used).

            Takes an int representing a functional class in the Artemis color
            scheme, and returns the appropriate colors.Color object
        """
        try:
            value = int(value)
        except ValueError:
            if value.count('.'):                           # dot-delimited
                value = int(artemis_color.split('.',1)[0]) # Use only first integer
            else:
                raise
        if value in self._artemis_colorscheme:
            return self._artemis_colorscheme[value][0]
        else:
            raise ValueError, "Artemis color out of range: %d" % value


    def get_colorscheme(self):
        """ get_colorscheme(self)

            Return the user-defined color scheme as a dictionary
        """
        return self._colorscheme


    def scheme_color(self, value):
        """ scheme_color(self, value)

            o value     An int representing a single color in the user-defined
                        color scheme

            Takes an int representing a user-defined color and returns the
            appropriate colors.Color object
        """
        if value in self._colorscheme:
            return self._colorscheme[value][0]
        else:
            raise ValueError, "Scheme color out of range: %d" % value


    def int255_color(self, values):
        """ int255_color(self, values)

            o values        A tuple of (red, green, blue) intensities as
                            integers in the range 0->255

            Takes a tuple of (red, green, blue) intensity values in the range
            0 -> 255 and returns an appropriate colors.Color object
        """
        red, green, blue = values
        factor = 1/255.
        red, green, blue = red * factor, green * factor, blue * factor
        return colors.Color(red, green, blue)


    def float1_color(self, values):
        """ float1_color(self, values)

            o values        A tuple of (red, green, blue) intensities as floats
                            in the range 0 -> 1

            Takes a tuple of (red, green, blue) intensity values in the range
            0 -> 1 and returns an appropriate colors.Color object
        """
        red, green, blue = values
        return colors.Color(red, green, blue)


################################################################################
# RUN AS SCRIPT
################################################################################

if __name__ == '__main__':

    # Test code
    gdct = ColorTranslator()
    print gdct.float1_color((0.5, 0.5, 0.5))
    print gdct.int255_color((1, 75, 240))
    print gdct.artemis_color(7)
    print gdct.scheme_color(2)

    print gdct.translate((0.5, 0.5, 0.5))
    print gdct.translate((1, 75, 240))
    print gdct.translate(7)
    print gdct.translate(2)

                
                        

