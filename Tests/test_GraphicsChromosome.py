#!/usr/bin/env python
# Copyright 2001 by Brad Chapman.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test the creation of Chromosome Graphics.

Provide tests for creating and manipulating pdf pictures of chromosomes.

This tests the Graphics.BasicChromosome classes and the
Graphics.DisplayRepresentation classes.
"""
# standard library
import os
import sys
import random
import cStringIO
import unittest

from Bio import MissingPythonDependencyError
try:
    # reportlab
    from reportlab.lib import colors
except:
    raise MissingPythonDependencyError(
        "Install reportlab if you want to use Bio.Graphics.")

# local stuff
from Bio.Graphics import BasicChromosome
from Bio.Graphics.DisplayRepresentation import ChromosomeCounts


# hold the chromosome info for testing
# the info is structured as (label, color, scale)
chr1_info = (("", None, 1),
             ("AC30001", None, 1),
             ("", colors.blue, 2),
             ("", colors.blue, 1.5),
             ("", colors.red, 2),
             ("AC12345", colors.blue, 1),
             ("", colors.blue, 1),
             ("", colors.blue, 1),
             ("", None, 1))

chr2_info = (("", None, 2),
             ("AC23456", None, 2),
             ("", colors.blue, 1),
             ("", colors.blue, 1),
             ("", colors.red, 1),
             ("AC00002", colors.blue, 1.5),
             ("", colors.blue, 1.5),
             ("", colors.blue, 2),
             ("", None, 1.5),
             ("", colors.red, 2),
             ("", colors.red, 1.5),
             ("AB0001", None, 1),
             ("", colors.red, 1))

chr3_info = (("", colors.green, 2),
             ("", colors.green, 1),
             ("AD00002", colors.blue, 1),
             ("", colors.blue, 1),
             ("", colors.red, 1))

chr4_info = (("", None, 1.5),
             ("", None, 1),
             ("AD11111", colors.blue, 2),
             ("", colors.blue, 1),
             ("", colors.red, 1.5),
             ("", colors.blue, 1),
             ("", colors.blue, 1),
             ("AC22222", colors.blue, 1))

all_chr_info = {"I" : chr1_info,
                "II" : chr2_info,
                "III" : chr3_info,
                "IV" : chr4_info}

def load_chromosome(chr_name):
    """Load a chromosome and all of its segments.
    """
    cur_chromosome = BasicChromosome.Chromosome(chr_name)

    chr_segment_info = all_chr_info[chr_name]

    for seg_info_num in range(len(chr_segment_info)):
        label, fill_color, scale = chr_segment_info[seg_info_num]
        # make the top and bottom telomeres
        if seg_info_num == 0:
            cur_segment = BasicChromosome.TelomereSegment()
        elif seg_info_num == len(chr_segment_info) - 1:
            cur_segment = BasicChromosome.TelomereSegment(1)
        # otherwise, they are just regular segments
        else:
            cur_segment = BasicChromosome.ChromosomeSegment()
        if label != "":
            cur_segment.label = label
        if fill_color is not None:
            cur_segment.fill_color = fill_color

        cur_segment.scale = scale

        cur_chromosome.add(cur_segment)

    # scale by the size of chromosome 2
    cur_chromosome.scale_num = 19

    return cur_chromosome

# --- stuff for generating random organisms
color_choices = (colors.red, colors.blue)
letter_choices = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
num_possible_segments = 500
color_prob = .3
id_prob = .025

def get_random_id():
    """Generate a random id number.
    """
    id = ''
    for n in range(6):
        letter = random.choice(letter_choices)
        id += letter

    return id

def load_random_chromosome(chr_name):
    """Generate a chromosome with random information about it.
    """
    cur_chromosome = BasicChromosome.Chromosome(chr_name)

    num_segments = random.randrange(num_possible_segments)
    for seg in range(num_segments):
        # make the top and bottom telomeres
        if seg == 0:
            cur_segment = BasicChromosome.TelomereSegment()
        elif seg == num_segments - 1:
            cur_segment = BasicChromosome.TelomereSegment(1)
        # otherwise, they are just regular segments
        else:
            cur_segment = BasicChromosome.ChromosomeSegment()
            
        color_chance = random.random()
        if color_chance <= color_prob:
            fill_color = random.choice(color_choices)
            cur_segment.fill_color = fill_color

        id_chance = random.random()
        if id_chance <= id_prob:
            id = get_random_id()
            cur_segment.label = id

        cur_chromosome.add(cur_segment)
        
    return cur_chromosome, num_segments

class OrganismGraphicTest(unittest.TestCase):
    """Test the creation of all chromosomes of an organism.
    """
    def setUp(self):
        self.test_file = os.path.join("Graphics", "organism.pdf")

    def test_simple_organism(self):
        """Test the basic functionality of drawing an organism.
        """
        pdf_organism = BasicChromosome.Organism()

        # add chromosomes
        for chr_name in ["I", "II", "III", "IV"]:
            cur_chromosome = load_chromosome(chr_name)
            pdf_organism.add(cur_chromosome)

        pdf_organism.draw(self.test_file, "Test organism")


    def _simple_organism(self, filename, format):
        """Output a simple organism to given format."""
        test_organism = BasicChromosome.Organism(format)
        test_file = os.path.join("Graphics", filename)

        # add chromosomes
        for chr_name in ["I", "II", "III", "IV"]:
            cur_chromosome = load_chromosome(chr_name)
            test_organism.add(cur_chromosome)
        test_organism.draw(test_file, "Test organism")
        
    def test_simple_organism_ps(self):
        """Output a simple organism to a postscript file.
        """
        self._simple_organism("organism.eps", "ps")

    def test_simple_organism_pdf(self):
        """Output a simple organism to a PDF file.
        """
        self._simple_organism("organism.pdf", "pdf")

    def test_simple_organism_svg(self):
        """Output a simple organism to an SVG file.
        """
        self._simple_organism("organism.svg", "svg")
     
    def test_random_organism(self):
        """Generate an organism with random chromosome info.
        """
        random_file = os.path.join("Graphics", "random_organism.pdf")
        pdf_organism = BasicChromosome.Organism()

        all_segs = []
        all_chrs = []
        
        num_chrs = random.randrange(1, 15)
        for chr_name in range(num_chrs):
            cur_chromosome, num_segs = load_random_chromosome(str(chr_name))
            all_chrs.append(cur_chromosome)
            all_segs.append(num_segs)

        # scale all of the chromosomes by the maximum number of segments
        max_segs = max(all_segs)
        for chr in all_chrs:
            chr.scale_num = max_segs
            pdf_organism.add(chr)

        pdf_organism.draw(random_file, "Randomly generated Organism")

    def test_widget(self):
        """Try widget derived functionality.
        """
        test_widget = BasicChromosome.ChromosomeSegment()

        expected_string = "chr_percent = 0.25"

        # trick to write the properties to a string
        save_stdout = sys.stdout
        new_stdout = cStringIO.StringIO()
        sys.stdout = new_stdout
        
        test_widget.dumpProperties()

        properties = new_stdout.getvalue()
        sys.stdout = save_stdout

        self.assertTrue(properties.find(expected_string) >= 0,
               "Unexpected results from dumpProperties: \n %s" % properties)

        properties = test_widget.getProperties()
        self.assertEqual(properties["label_size"], 6, 
               "Unexpected results from getProperties: %s" % properties)

        test_widget.setProperties({"start_x_position" : 12})
        self.assertEqual(test_widget.start_x_position, 12,
               "setProperties doesn't seem to work right: %s" \
               % test_widget.start_x_position)
        
class ChromosomeCountTest(unittest.TestCase):
    """Test the display representation for simple counts on a chromosome.
    """
    def setUp(self):
        self.names = ["Bob", "Dylan", "Doesn't", "Like", "Spam"]
        self.count_display = ChromosomeCounts(self.names)

    def test_add_count(self):
        """Add counts to specific chromosome segments.
        """
        self.count_display.add_count(self.names[1])
        self.count_display.add_count(self.names[2], 5)

        try:
            self.count_display.add_count("Non-existent")
            raise AssertionError("Didn't raise a KeyError on a fake key")
        except KeyError:
            pass

    def test_add_label(self):
        """Add labels to chromosome segments.
        """
        self.count_display.add_label(self.names[1], "Rules")

        try:
            self.count_display.add_label("Non-existent", "Elephant")
            raise AssertionError("Didn't raise a KeyError on a fake key")
        except KeyError:
            pass

    def test_set_scale(self):
        """Set the scale for a chromosome segment.
        """
        self.count_display.set_scale(self.names[1], 1.5)

        try:
            self.count_display.set_scale("Non-existant", 5)
            raise AssertionError("Didn't raise a KeyError on a fake key.")
        except KeyError:
            pass

    def test_color_from_count(self):
        """Retrieve a color from a count number with the default color scheme.
        """
        test_color = self.count_display._color_from_count(3)
        assert test_color == colors.blue, "Unexpected color %s" % test_color

        test_color = self.count_display._color_from_count(9)
        assert test_color == colors.red, "Unexpected color %s" % test_color

        try:
            self.count_display._color_from_count(200)
            raise AssertionError("Didn't raise an error for bad count.")
        except ValueError:
            pass

    def test_fill_chromosome(self):
        """Test filling out the information on a chromosome.
        """
        test_chr = BasicChromosome.Chromosome("1")
        self.count_display.add_count(self.names[2], 5)
        self.count_display.add_count(self.names[1], 2)
        self.count_display.add_label(self.names[3], "Test-Label")

        new_chr = self.count_display.fill_chromosome(test_chr)

    def test_get_segment_info(self):
        """Test retrieval of segment information.
        """
        test_count_num = 1
        test_count_value = 5
        
        test_label_num = 3
        test_label_value = "BigBird"
        
        self.count_display.add_count(self.names[test_count_num],
                                     test_count_value)
        self.count_display.add_label(self.names[test_label_num],
                                     test_label_value)

        seg_info = self.count_display.get_segment_info()
        
        assert seg_info[test_count_num][0] == test_count_value, \
               "Did not set and retrieve counts correctly."
        assert seg_info[test_label_num][1] == test_label_value, \
               "Did not set and retrieve label correctly."

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
