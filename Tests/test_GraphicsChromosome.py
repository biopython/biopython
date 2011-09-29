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


class OrganismSubAnnotationsTest(unittest.TestCase):
    """Test sub-annotations on a segment."""
    def test_simple_tRNA(self):
        """Test sub-annotations on a genome segment, tRNA for cress"""
        f1 = [(111889, 111961, -1, 'AT1G01270'), (306383, 306456, 1, 'AT1G01870'), (309274, 309347, -1, 'AT1G01890'), (515493, 515566, 1, 'AT1G02480'), (552639, 552711, 1, 'AT1G02600'), (604401, 604474, 1, 'AT1G02760'), (877648, 877720, 1, 'AT1G03515'), (892513, 892585, 1, 'AT1G03570'), (909809, 909882, -1, 'AT1G03640'), (1159021, 1159092, 1, 'AT1G04320'), (1324921, 1324959, 1, 'AT1G04720'), (1583770, 1583844, -1, 'AT1G05390'), (1817398, 1817470, 1, 'AT1G05980'), (1978082, 1978156, 1, 'AT1G06480'), (2025354, 2025427, 1, 'AT1G06610'), (2107396, 2107467, -1, 'AT1G06860'), (2111146, 2111217, -1, 'AT1G06880'), (2177883, 2177957, 1, 'AT1G07100'), (2334818, 2334891, 1, 'AT1G07580'), (2406830, 2406902, -1, 'AT1G07760'), (2588521, 2588593, 1, 'AT1G08240'), (2846538, 2846611, -1, 'AT1G08870'), (2879305, 2879377, 1, 'AT1G08950'), (2939418, 2939490, 1, 'AT1G09110'), (3431185, 3431257, -1, 'AT1G10440'), (3676606, 3676644, 1, 'AT1G11010'), (3678774, 3678848, -1, 'AT1G11030'), (3881528, 3881608, 1, 'AT1G11550'), (3914628, 3914700, -1, 'AT1G11640'), (4266985, 4267059, -1, 'AT1G12510'), (4285884, 4285956, -1, 'AT1G12590'), (4440211, 4440284, 1, 'AT1G13010'), (4522705, 4522779, -1, 'AT1G13240'), (4709631, 4709703, 1, 'AT1G13720'), (4741995, 4742068, 1, 'AT1G13840'), (4743091, 4743164, 1, 'AT1G13850'), (5189681, 5189755, -1, 'AT1G15090'), (5309641, 5309713, -1, 'AT1G15450'), (5380901, 5380983, 1, 'AT1G15650'), (5518055, 5518128, -1, 'AT1G16100'), (5619464, 5619537, -1, 'AT1G16450'), (6038749, 6038831, 1, 'AT1G17570'), (6075812, 6075884, 1, 'AT1G17660'), (6075937, 6076011, -1, 'AT1G17670'), (6345756, 6345828, 1, 'AT1G18430'), (6488645, 6488726, 1, 'AT1G18820'), (6948850, 6948934, -1, 'AT1G20040'), (6995272, 6995344, -1, 'AT1G20170'), (7004504, 7004576, 1, 'AT1G20210'), (7016506, 7016579, 1, 'AT1G20250'), (7082657, 7082729, 1, 'AT1G20420'), (7242749, 7242821, -1, 'AT1G20820'), (7499721, 7499793, -1, 'AT1G21420'), (7656108, 7656180, -1, 'AT1G21800'), (7884405, 7884443, -1, 'AT1G22320'), (8520278, 8520352, -1, 'AT1G24080'), (9143796, 9143870, 1, 'AT1G26430'), (9158169, 9158242, 1, 'AT1G26490'), (10089422, 10089494, 1, 'AT1G28720'), (10089883, 10089955, 1, 'AT1G28730'), (10090353, 10090425, 1, 'AT1G28740'), (10090754, 10090826, 1, 'AT1G28750'), (10092310, 10092382, 1, 'AT1G28770'), (10092786, 10092858, 1, 'AT1G28780'), (10093294, 10093366, 1, 'AT1G28790'), (10093731, 10093803, 1, 'AT1G28800'), (10094158, 10094230, 1, 'AT1G28810'), (10096936, 10097008, 1, 'AT1G28820'), (10097099, 10097171, 1, 'AT1G28830'), (10097703, 10097775, 1, 'AT1G28840'), (10098638, 10098710, 1, 'AT1G28850'), (10099064, 10099136, 1, 'AT1G28860'), (10099410, 10099482, 1, 'AT1G28870'), (10099812, 10099884, 1, 'AT1G28880'), (10100258, 10100330, 1, 'AT1G28890'), (10101013, 10101085, 1, 'AT1G28900'), (10101585, 10101657, 1, 'AT1G28910'), (10101978, 10102050, 1, 'AT1G28920'), (10106075, 10106147, 1, 'AT1G28930'), (10106513, 10106585, 1, 'AT1G28940'), (10106883, 10106955, 1, 'AT1G28950'), (10107634, 10107706, 1, 'AT1G28970'), (10108374, 10108446, 1, 'AT1G28980'), (10108695, 10108767, 1, 'AT1G28990'), (10207291, 10207364, -1, 'AT1G29210'), (10756703, 10756776, 1, 'AT1G30430'), (10963553, 10963627, -1, 'AT1G30830'), (11104093, 11104167, 1, 'AT1G31110'), (11797227, 11797265, -1, 'AT1G32620'), (12097258, 12097327, -1, 'AT1G33370'), (13687637, 13687710, 1, 'AT1G36350'), (15733055, 15733127, -1, 'AT1G42120'), (16588144, 16588216, -1, 'AT1G43820'), (17159046, 17159118, 1, 'AT1G45234'), (17159799, 17159871, 1, 'AT1G45236'), (17160970, 17161042, 1, 'AT1G45238'), (17161418, 17161490, 1, 'AT1G45240'), (17162967, 17163039, 1, 'AT1G45242'), (17163408, 17163480, 1, 'AT1G45244'), (17164461, 17164533, 1, 'AT1G45246'), (17735509, 17735582, 1, 'AT1G48080'), (18139265, 18139337, -1, 'AT1G49020'), (18234146, 18234220, -1, 'AT1G49280'), (18312570, 18312607, 1, 'AT1G49460'), (18391469, 18391542, 1, 'AT1G49690'), (18556666, 18556746, 1, 'AT1G50070'), (18561567, 18561647, 1, 'AT1G50100'), (19428223, 19428297, 1, 'AT1G52170'), (19502087, 19502161, -1, 'AT1G52350'), (19688850, 19688887, -1, 'AT1G52860'), (19851640, 19851714, 1, 'AT1G53220'), (19929506, 19929578, -1, 'AT1G53410'), (20416594, 20416667, -1, 'AT1G54670'), (20794976, 20795058, 1, 'AT1G55625'), (21272451, 21272533, 1, 'AT1G56730'), (21272786, 21272823, 1, 'AT1G56740'), (21273216, 21273253, 1, 'AT1G56750'), (21273960, 21274042, 1, 'AT1G56760'), (21274295, 21274332, 1, 'AT1G56770'), (21274725, 21274762, 1, 'AT1G56780'), (21275469, 21275551, 1, 'AT1G56790'), (21275804, 21275841, 1, 'AT1G56800'), (21276234, 21276271, 1, 'AT1G56810'), (21276978, 21277060, 1, 'AT1G56820'), (21277313, 21277350, 1, 'AT1G56830'), (21277743, 21277780, 1, 'AT1G56840'), (21278487, 21278569, 1, 'AT1G56850'), (21278822, 21278859, 1, 'AT1G56860'), (21279273, 21279310, 1, 'AT1G56870'), (21280016, 21280098, 1, 'AT1G56880'), (21280351, 21280388, 1, 'AT1G56890'), (21280781, 21280818, 1, 'AT1G56900'), (21281525, 21281607, 1, 'AT1G56910'), (21281860, 21281897, 1, 'AT1G56920'), (21282311, 21282348, 1, 'AT1G56930'), (21283054, 21283136, 1, 'AT1G56940'), (21283384, 21283421, 1, 'AT1G56950'), (21283842, 21283879, 1, 'AT1G56960'), (21284586, 21284668, 1, 'AT1G56970'), (21284916, 21284953, 1, 'AT1G56980'), (21285374, 21285411, 1, 'AT1G56990'), (21286118, 21286200, 1, 'AT1G57000'), (21286448, 21286485, 1, 'AT1G57010'), (21286906, 21286943, 1, 'AT1G57020'), (21287650, 21287732, 1, 'AT1G57030'), (21287980, 21288017, 1, 'AT1G57040'), (21288438, 21288475, 1, 'AT1G57050'), (21289183, 21289265, 1, 'AT1G57060'), (21289513, 21289550, 1, 'AT1G57070'), (21289970, 21290007, 1, 'AT1G57080'), (21290714, 21290796, 1, 'AT1G57090'), (21291044, 21291081, 1, 'AT1G57100'), (21291501, 21291538, 1, 'AT1G57110'), (21292245, 21292327, 1, 'AT1G57120'), (21292574, 21292611, 1, 'AT1G57130'), (21293032, 21293069, 1, 'AT1G57140'), (21293776, 21293858, 1, 'AT1G57150'), (21294109, 21294146, 1, 'AT1G57160'), (21294567, 21294604, 1, 'AT1G57170'), (21295125, 21295207, 1, 'AT1G57180'), (21295455, 21295492, 1, 'AT1G57190'), (21295912, 21295949, 1, 'AT1G57200'), (21296656, 21296738, 1, 'AT1G57210'), (21296989, 21297026, 1, 'AT1G57220'), (21297447, 21297484, 1, 'AT1G57230'), (21298005, 21298087, 1, 'AT1G57240'), (21298335, 21298372, 1, 'AT1G57250'), (21298792, 21298829, 1, 'AT1G57260'), (21299536, 21299618, 1, 'AT1G57270'), (21299869, 21299906, 1, 'AT1G57280'), (21300327, 21300364, 1, 'AT1G57290'), (21300885, 21300967, 1, 'AT1G57300'), (21301215, 21301252, 1, 'AT1G57310'), (21301673, 21301710, 1, 'AT1G57320'), (21302417, 21302499, 1, 'AT1G57330'), (21302750, 21302787, 1, 'AT1G57340'), (21303208, 21303245, 1, 'AT1G57350'), (21303766, 21303848, 1, 'AT1G57360'), (21304096, 21304133, 1, 'AT1G57370'), (21304554, 21304591, 1, 'AT1G57380'), (21305298, 21305380, 1, 'AT1G57390'), (21305631, 21305668, 1, 'AT1G57400'), (21306089, 21306126, 1, 'AT1G57410'), (21306647, 21306729, 1, 'AT1G57420'), (21306981, 21307018, 1, 'AT1G57430'), (21307441, 21307478, 1, 'AT1G57440'), (21308184, 21308268, 1, 'AT1G57450'), (21308520, 21308557, 1, 'AT1G57460'), (21308975, 21309012, 1, 'AT1G57470'), (21309719, 21309801, 1, 'AT1G57480'), (21310053, 21310090, 1, 'AT1G57490'), (21310513, 21310550, 1, 'AT1G57500'), (21311256, 21311340, 1, 'AT1G57510'), (21311592, 21311629, 1, 'AT1G57520'), (21312051, 21312088, 1, 'AT1G57530'), (21377983, 21378054, -1, 'AT1G57710'), (21887507, 21887589, -1, 'AT1G59570'), (22044276, 22044348, -1, 'AT1G59880'), (22317078, 22317149, -1, 'AT1G60580'), (22398301, 22398372, -1, 'AT1G60820'), (22401256, 22401327, -1, 'AT1G60840'), (22431831, 22431902, 1, 'AT1G60910'), (22481437, 22481511, -1, 'AT1G61020'), (22870422, 22870494, -1, 'AT1G61880'), (22890754, 22890834, 1, 'AT1G61910'), (23562849, 23562921, -1, 'AT1G63510'), (23671147, 23671219, -1, 'AT1G63790'), (23806215, 23806299, 1, 'AT1G64120'), (23936799, 23936872, 1, 'AT1G64420'), (24490654, 24490736, -1, 'AT1G65830'), (25833316, 25833388, 1, 'AT1G68770'), (25890198, 25890272, 1, 'AT1G68860'), (25931858, 25931931, 1, 'AT1G68950'), (25935739, 25935812, -1, 'AT1G68970'), (25944826, 25944898, 1, 'AT1G69000'), (25993392, 25993466, 1, 'AT1G69130'), (26053140, 26053214, 1, 'AT1G69300'), (26385816, 26385888, -1, 'AT1G70050'), (26977050, 26977121, 1, 'AT1G71700'), (27397046, 27397128, 1, 'AT1G72780'), (27792643, 27792715, 1, 'AT1G73900'), (28024043, 28024124, -1, 'AT1G74570'), (28031620, 28031701, 1, 'AT1G74610'), (28188192, 28188264, 1, 'AT1G75070'), (28377149, 28377222, -1, 'AT1G75570'), (28411644, 28411717, 1, 'AT1G75650'), (28444549, 28444621, 1, 'AT1G75740'), (28523645, 28523717, -1, 'AT1G75970'), (28531427, 28531499, 1, 'AT1G76000'), (28639585, 28639667, 1, 'AT1G76330'), (28952447, 28952519, -1, 'AT1G77040'), (29007098, 29007180, -1, 'AT1G77190'), (29147983, 29148055, -1, 'AT1G77560'), (29448865, 29448903, -1, 'AT1G78250'), (29809015, 29809088, 1, 'AT1G79240'), (29838009, 29838081, 1, 'AT1G79290'), (29838610, 29838682, 1, 'AT1G79300'), (30088888, 30088962, -1, 'AT1G79980'), (30178905, 30178977, -1, 'AT1G80250'), (30242675, 30242757, 1, 'AT1G80430')]
        f2 = [(102063, 102137, 1, 'AT2G01160'), (706794, 706867, 1, 'AT2G02600'), (846853, 846926, -1, 'AT2G02900'), (1054714, 1054787, -1, 'AT2G03490'), (1113980, 1114052, -1, 'AT2G03660'), (1123386, 1123458, -1, 'AT2G03700'), (1154381, 1154454, 1, 'AT2G03790'), (3239653, 3239725, -1, 'AT2G07742'), (3255828, 3255902, -1, 'AT2G07743'), (3268803, 3268883, 1, 'AT2G07745'), (3276436, 3276508, 1, 'AT2G07746'), (3280859, 3280933, 1, 'AT2G07748'), (3290962, 3291034, 1, 'AT2G07778'), (3303240, 3303312, -1, 'AT2G07752'), (3303350, 3303425, -1, 'AT2G07753'), (3303781, 3303819, -1, 'AT2G07754'), (3328666, 3328739, -1, 'AT2G07755'), (3332674, 3332756, 1, 'AT2G07792'), (3369350, 3369437, 1, 'AT2G07793'), (3383400, 3383474, -1, 'AT2G07794'), (3444359, 3444431, -1, 'AT2G07756'), (3452973, 3453060, 1, 'AT2G07757'), (3462074, 3462148, 1, 'AT2G07758'), (3494378, 3494416, 1, 'AT2G07759'), (3494772, 3494847, 1, 'AT2G07761'), (3495008, 3495083, 1, 'AT2G07762'), (3495438, 3495509, 1, 'AT2G07763'), (3496436, 3496508, 1, 'AT2G07764'), (3497354, 3497437, 1, 'AT2G07765'), (3503518, 3503605, 1, 'AT2G07766'), (6953924, 6953961, -1, 'AT2G15950'), (7046175, 7046247, 1, 'AT2G16240'), (7749793, 7749867, 1, 'AT2G17810'), (7962758, 7962832, -1, 'AT2G18310'), (9144435, 9144507, 1, 'AT2G21360'), (9241319, 9241356, -1, 'AT2G21570'), (9273888, 9273969, -1, 'AT2G21670'), (9277742, 9277814, -1, 'AT2G21700'), (9291113, 9291185, 1, 'AT2G21760'), (9400749, 9400823, 1, 'AT2G22110'), (9456888, 9456962, -1, 'AT2G22220'), (9472660, 9472733, -1, 'AT2G22280'), (9509359, 9509433, 1, 'AT2G22380'), (9598106, 9598179, 1, 'AT2G22580'), (9810296, 9810368, -1, 'AT2G23020'), (10066525, 10066597, -1, 'AT2G23650'), (10380655, 10380728, 1, 'AT2G24380'), (10820917, 10820990, 1, 'AT2G25400'), (11122756, 11122837, -1, 'AT2G26090'), (11781928, 11782000, -1, 'AT2G27560'), (11871230, 11871302, -1, 'AT2G27850'), (12336079, 12336151, 1, 'AT2G28730'), (12346827, 12346899, 1, 'AT2G28770'), (12478849, 12478921, -1, 'AT2G29030'), (12645232, 12645305, -1, 'AT2G29520'), (12888667, 12888738, 1, 'AT2G30180'), (12889810, 12889881, 1, 'AT2G30190'), (12983024, 12983095, -1, 'AT2G30450'), (13144312, 13144385, -1, 'AT2G30850'), (13658350, 13658425, 1, 'AT2G32110'), (14054465, 14054503, -1, 'AT2G33140'), (14250206, 14250278, 1, 'AT2G33650'), (14251774, 14251846, 1, 'AT2G33660'), (14357464, 14357536, 1, 'AT2G33890'), (14358437, 14358509, 1, 'AT2G33900'), (14359269, 14359341, 1, 'AT2G33910'), (14360221, 14360293, 1, 'AT2G33920'), (14360734, 14360806, 1, 'AT2G33930'), (14361176, 14361248, 1, 'AT2G33940'), (14362215, 14362287, 1, 'AT2G33950'), (14363133, 14363205, 1, 'AT2G33960'), (14363599, 14363671, 1, 'AT2G33970'), (14750553, 14750627, -1, 'AT2G34950'), (14757142, 14757213, 1, 'AT2G34985'), (14847685, 14847723, 1, 'AT2G35220'), (15175940, 15176014, 1, 'AT2G36140'), (15176656, 15176736, 1, 'AT2G36150'), (15215480, 15215517, -1, 'AT2G36280'), (15327312, 15327395, 1, 'AT2G36510'), (15327463, 15327546, -1, 'AT2G36520'), (15353238, 15353311, 1, 'AT2G36600'), (15477287, 15477324, -1, 'AT2G36860'), (15923894, 15923967, 1, 'AT2G38030'), (16525641, 16525713, -1, 'AT2G39600'), (16525846, 16525918, 1, 'AT2G39610'), (16646857, 16646929, -1, 'AT2G39860'), (17545780, 17545862, -1, 'AT2G42020'), (17667855, 17667926, 1, 'AT2G42420'), (17880766, 17880839, 1, 'AT2G42970'), (18002649, 18002721, -1, 'AT2G43300'), (18317052, 18317134, -1, 'AT2G44320'), (18576985, 18577058, 1, 'AT2G45020'), (18710751, 18710824, 1, 'AT2G45390'), (18963713, 18963786, 1, 'AT2G46120'), (19351496, 19351569, 1, 'AT2G47100'), (19566924, 19566995, -1, 'AT2G47740')]
        f3 = [(259640, 259712, 1, 'AT3G01705'), (469666, 469740, 1, 'AT3G02315'), (476808, 476880, 1, 'AT3G02335'), (586092, 586174, 1, 'AT3G02715'), (981975, 982047, 1, 'AT3G03845'), (984105, 984177, 1, 'AT3G03852'), (1220234, 1220307, 1, 'AT3G04525'), (1601343, 1601415, -1, 'AT3G05525'), (1707743, 1707815, -1, 'AT3G05755'), (1738796, 1738870, 1, 'AT3G05835'), (1843329, 1843400, -1, 'AT3G06105'), (1920038, 1920110, -1, 'AT3G06335'), (2104961, 2105033, -1, 'AT3G06665'), (2222251, 2222324, 1, 'AT3G07025'), (2232470, 2232506, -1, 'AT3G07055'), (2253680, 2253762, -1, 'AT3G07115'), (2285607, 2285679, 1, 'AT3G07185'), (2918418, 2918492, -1, 'AT3G09505'), (2944616, 2944698, 1, 'AT3G09585'), (2945700, 2945782, -1, 'AT3G09595'), (3090548, 3090631, 1, 'AT3G10015'), (3096220, 3096293, 1, 'AT3G10035'), (3238371, 3238407, -1, 'AT3G10415'), (3535151, 3535224, 1, 'AT3G11285'), (3575849, 3575923, 1, 'AT3G11395'), (3622697, 3622769, -1, 'AT3G11505'), (3942012, 3942084, 1, 'AT3G12385'), (3995103, 3995176, -1, 'AT3G12585'), (4254534, 4254615, 1, 'AT3G13223'), (4330778, 4330850, 1, 'AT3G13335'), (4998147, 4998219, 1, 'AT3G14855'), (5068300, 5068374, -1, 'AT3G15055'), (5275155, 5275228, 1, 'AT3G15585'), (5632857, 5632930, 1, 'AT3G16552'), (6483945, 6484019, -1, 'AT3G18815'), (6540636, 6540673, 1, 'AT3G18952'), (6663713, 6663786, 1, 'AT3G19235'), (7104314, 7104398, 1, 'AT3G20365'), (7224223, 7224296, -1, 'AT3G20655'), (7319582, 7319664, -1, 'AT3G20885'), (7567399, 7567471, -1, 'AT3G21475'), (9373610, 9373684, -1, 'AT3G25715'), (9840420, 9840494, 1, 'AT3G26747'), (10211564, 10211636, 1, 'AT3G27555'), (10319498, 10319570, 1, 'AT3G27825'), (10325875, 10325947, 1, 'AT3G27845'), (10753667, 10753740, 1, 'AT3G28685'), (10760629, 10760702, -1, 'AT3G28695'), (11076814, 11076886, 1, 'AT3G29095'), (11961645, 11961718, 1, 'AT3G30345'), (16438025, 16438097, -1, 'AT3G44955'), (16896875, 16896949, 1, 'AT3G45935'), (16902623, 16902697, 1, 'AT3G45955'), (16905147, 16905221, 1, 'AT3G45965'), (17160736, 17160808, 1, 'AT3G46585'), (17275564, 17275646, 1, 'AT3G46875'), (17905395, 17905467, 1, 'AT3G48275'), (17985575, 17985611, -1, 'AT3G48515'), (18080062, 18080134, 1, 'AT3G48745'), (18518796, 18518870, 1, 'AT3G49925'), (18755788, 18755860, -1, 'AT3G50505'), (18837020, 18837092, 1, 'AT3G50665'), (18907851, 18907924, 1, 'AT3G50835'), (18928413, 18928487, 1, 'AT3G50895'), (19008621, 19008694, -1, 'AT3G51135'), (19044371, 19044443, -1, 'AT3G51265'), (19403651, 19403723, -1, 'AT3G52285'), (19420345, 19420417, -1, 'AT3G52345'), (19511965, 19512045, 1, 'AT3G52565'), (19566013, 19566085, 1, 'AT3G52765'), (19648105, 19648188, 1, 'AT3G52955'), (19935354, 19935426, 1, 'AT3G53775'), (19995918, 19995989, 1, 'AT3G53965'), (20704664, 20704736, 1, 'AT3G55735'), (20720151, 20720223, 1, 'AT3G55795'), (20824495, 20824568, -1, 'AT3G56085'), (21498293, 21498375, 1, 'AT3G58035'), (21553258, 21553329, 1, 'AT3G58165'), (21970486, 21970557, 1, 'AT3G59415'), (22149699, 22149773, 1, 'AT3G59923'), (22149823, 22149895, -1, 'AT3G59926'), (22197810, 22197892, -1, 'AT3G60075'), (22481215, 22481288, -1, 'AT3G60805'), (22622384, 22622465, 1, 'AT3G61105'), (22786896, 22786969, 1, 'AT3G61545'), (22853496, 22853567, 1, 'AT3G61715'), (22871101, 22871174, 1, 'AT3G61755'), (22892781, 22892853, 1, 'AT3G61825'), (23047854, 23047927, 1, 'AT3G62245'), (23062444, 23062517, -1, 'AT3G62285'), (23221682, 23221753, 1, 'AT3G62735'), (23296567, 23296640, -1, 'AT3G63003'), (23296728, 23296801, -1, 'AT3G63006')]
        f4 = [(33799, 33872, 1, 'AT4G00085'), (424716, 424788, -1, 'AT4G00985'), (562560, 562634, -1, 'AT4G01355'), (611865, 611932, -1, 'AT4G01455'), (808269, 808342, -1, 'AT4G01865'), (901175, 901247, 1, 'AT4G02055'), (1390894, 1390966, 1, 'AT4G03135'), (1442004, 1442076, 1, 'AT4G03285'), (1501605, 1501677, 1, 'AT4G03405'), (1520781, 1520854, -1, 'AT4G03435'), (5268124, 5268210, -1, 'AT4G08345'), (6646425, 6646496, 1, 'AT4G10815'), (6819287, 6819324, 1, 'AT4G11177'), (6837555, 6837639, -1, 'AT4G11213'), (6837769, 6837853, -1, 'AT4G11216'), (6905479, 6905552, -1, 'AT4G11355'), (6944721, 6944793, 1, 'AT4G11405'), (7185697, 7185771, 1, 'AT4G11985'), (7232792, 7232865, -1, 'AT4G12065'), (7256408, 7256481, 1, 'AT4G12115'), (7341420, 7341494, -1, 'AT4G12405'), (7730956, 7731037, 1, 'AT4G13265'), (7814197, 7814270, 1, 'AT4G13445'), (8255695, 8255767, 1, 'AT4G14345'), (8301720, 8301794, -1, 'AT4G14415'), (8979656, 8979729, 1, 'AT4G15775'), (9108317, 9108391, 1, 'AT4G16105'), (9191590, 9191663, 1, 'AT4G16235'), (9287230, 9287304, 1, 'AT4G16465'), (9289706, 9289787, 1, 'AT4G16475'), (9815215, 9815287, -1, 'AT4G17612'), (9873524, 9873596, -1, 'AT4G17765'), (9978117, 9978189, -1, 'AT4G17975'), (10093077, 10093157, -1, 'AT4G18255'), (10302011, 10302084, 1, 'AT4G18725'), (10325975, 10326047, -1, 'AT4G18815'), (10878733, 10878807, -1, 'AT4G20115'), (11774472, 11774508, -1, 'AT4G22265'), (11910299, 11910373, 1, 'AT4G22635'), (11954751, 11954824, -1, 'AT4G22754'), (11974951, 11975032, 1, 'AT4G22785'), (12320119, 12320203, 1, 'AT4G23635'), (12429608, 12429681, 1, 'AT4G23915'), (12486211, 12486282, -1, 'AT4G24025'), (12686148, 12686230, 1, 'AT4G24565'), (13006243, 13006316, -1, 'AT4G25435'), (13058840, 13058922, -1, 'AT4G25585'), (13076582, 13076666, -1, 'AT4G25635'), (13285431, 13285503, -1, 'AT4G26225'), (13336345, 13336419, -1, 'AT4G26375'), (13341501, 13341575, -1, 'AT4G26385'), (13454562, 13454635, 1, 'AT4G26675'), (13704787, 13704860, 1, 'AT4G27395'), (13882922, 13882994, -1, 'AT4G27875'), (13885196, 13885269, -1, 'AT4G27885'), (14032495, 14032567, 1, 'AT4G28362'), (14267286, 14267368, 1, 'AT4G28915'), (14470283, 14470355, 1, 'AT4G29415'), (15120655, 15120728, 1, 'AT4G31075'), (15183089, 15183162, 1, 'AT4G31265'), (15345717, 15345753, -1, 'AT4G31695'), (15430229, 15430303, -1, 'AT4G31895'), (15576655, 15576728, 1, 'AT4G32265'), (15671398, 15671469, 1, 'AT4G32475'), (15804553, 15804635, 1, 'AT4G32765'), (16304128, 16304201, 1, 'AT4G34035'), (16454700, 16454773, -1, 'AT4G34415'), (16556627, 16556700, 1, 'AT4G34695'), (16655290, 16655364, 1, 'AT4G34975'), (17130054, 17130127, 1, 'AT4G36197'), (17149473, 17149545, 1, 'AT4G36245'), (17276705, 17276779, -1, 'AT4G36635'), (17500800, 17500872, -1, 'AT4G37175'), (18254982, 18255018, -1, 'AT4G39195'), (18293773, 18293845, 1, 'AT4G39345'), (18395021, 18395093, 1, 'AT4G39615'), (18411258, 18411332, 1, 'AT4G39672'), (18501705, 18501778, -1, 'AT4G39865'), (18542164, 18542238, 1, 'AT4G39985')]
        f5 = [(150353, 150426, -1, 'AT5G01365'), (389889, 389960, -1, 'AT5G02025'), (508427, 508500, -1, 'AT5G02385'), (530819, 530893, 1, 'AT5G02435'), (559327, 559399, -1, 'AT5G02505'), (588890, 588964, -1, 'AT5G02615'), (614641, 614723, 1, 'AT5G02725'), (642397, 642479, -1, 'AT5G02815'), (858534, 858571, 1, 'AT5G03445'), (862395, 862468, -1, 'AT5G03452'), (970797, 970878, -1, 'AT5G03705'), (984365, 984448, 1, 'AT5G03745'), (998940, 999013, 1, 'AT5G03775'), (1742692, 1742765, 1, 'AT5G05795'), (1788651, 1788723, 1, 'AT5G05945'), (1804616, 1804690, 1, 'AT5G05985'), (1853302, 1853382, -1, 'AT5G06125'), (2060153, 2060235, -1, 'AT5G06685'), (2212678, 2212749, -1, 'AT5G07135'), (2309512, 2309549, -1, 'AT5G07315'), (2411148, 2411232, 1, 'AT5G07625'), (2432263, 2432336, -1, 'AT5G07675'), (2587826, 2587899, -1, 'AT5G08075'), (2898867, 2898951, -1, 'AT5G09345'), (2993327, 2993401, 1, 'AT5G09655'), (3030817, 3030890, -1, 'AT5G09755'), (3118377, 3118458, 1, 'AT5G09975'), (3212351, 3212424, -1, 'AT5G10235'), (3287553, 3287635, -1, 'AT5G10455'), (3324702, 3324775, 1, 'AT5G10525'), (3578295, 3578367, -1, 'AT5G11225'), (3617058, 3617130, 1, 'AT5G11325'), (3669000, 3669073, -1, 'AT5G11475'), (4471050, 4471122, 1, 'AT5G13845'), (4530475, 4530548, 1, 'AT5G14035'), (4673902, 4673974, 1, 'AT5G14495'), (4929562, 4929636, 1, 'AT5G15175'), (5157641, 5157715, 1, 'AT5G15805'), (5161514, 5161586, 1, 'AT5G15815'), (5358918, 5359000, 1, 'AT5G16375'), (5962699, 5962771, -1, 'AT5G18005'), (5965972, 5966044, -1, 'AT5G18015'), (5984378, 5984450, 1, 'AT5G18085'), (6258146, 6258218, 1, 'AT5G18755'), (6401240, 6401311, 1, 'AT5G19095'), (7073531, 7073603, -1, 'AT5G20852'), (7073944, 7074016, -1, 'AT5G20854'), (7074357, 7074429, -1, 'AT5G20856'), (7074773, 7074845, -1, 'AT5G20858'), (7222059, 7222131, -1, 'AT5G21378'), (7387890, 7387962, 1, 'AT5G22315'), (7981400, 7981472, 1, 'AT5G23665'), (8906418, 8906502, 1, 'AT5G25585'), (8946826, 8946899, -1, 'AT5G25625'), (9815405, 9815477, -1, 'AT5G27715'), (11802284, 11802356, 1, 'AT5G32017'), (13823211, 13823284, -1, 'AT5G35605'), (15049737, 15049811, -1, 'AT5G37795'), (15242547, 15242621, 1, 'AT5G38155'), (15593086, 15593160, 1, 'AT5G38905'), (15844253, 15844325, -1, 'AT5G39535'), (15993514, 15993587, 1, 'AT5G39895'), (16256865, 16256937, -1, 'AT5G40545'), (16427812, 16427893, 1, 'AT5G40945'), (16524760, 16524832, -1, 'AT5G41265'), (16655393, 16655477, 1, 'AT5G41605'), (16684663, 16684735, -1, 'AT5G41675'), (17476402, 17476475, -1, 'AT5G43455'), (17512768, 17512839, -1, 'AT5G43535'), (17856811, 17856883, -1, 'AT5G44283'), (17894906, 17894979, -1, 'AT5G44375'), (18058014, 18058088, 1, 'AT5G44705'), (18560206, 18560278, -1, 'AT5G45715'), (18576071, 18576143, 1, 'AT5G45745'), (18715888, 18715960, -1, 'AT5G46105'), (18807534, 18807614, 1, 'AT5G46325'), (18924749, 18924821, 1, 'AT5G46595'), (19658828, 19658900, 1, 'AT5G48465'), (19761400, 19761472, -1, 'AT5G48675'), (19820360, 19820398, 1, 'AT5G48835'), (20064048, 20064120, 1, 'AT5G49435'), (20692447, 20692519, 1, 'AT5G50805'), (20758903, 20758940, -1, 'AT5G50995'), (20773555, 20773637, 1, 'AT5G51055'), (21275059, 21275141, -1, 'AT5G52355'), (21318105, 21318189, -1, 'AT5G52495'), (21418369, 21418441, 1, 'AT5G52815'), (21740339, 21740410, -1, 'AT5G53487'), (22091631, 22091704, 1, 'AT5G54365'), (22094087, 22094160, 1, 'AT5G54375'), (22304851, 22304923, -1, 'AT5G54865'), (22355897, 22355970, -1, 'AT5G55045'), (22357726, 22357799, -1, 'AT5G55055'), (22501995, 22502068, -1, 'AT5G55505'), (22845356, 22845430, 1, 'AT5G56365'), (22973066, 22973138, 1, 'AT5G56745'), (23071996, 23072070, -1, 'AT5G56975'), (23463219, 23463291, 1, 'AT5G57885'), (23661936, 23662018, 1, 'AT5G58495'), (23861431, 23861503, 1, 'AT5G59055'), (23971167, 23971239, 1, 'AT5G59385'), (23974655, 23974727, 1, 'AT5G59395'), (24157171, 24157245, -1, 'AT5G59945'), (24279805, 24279886, 1, 'AT5G60285'), (24547401, 24547474, 1, 'AT5G60963'), (24548892, 24548964, 1, 'AT5G60966'), (24684507, 24684579, 1, 'AT5G61345'), (24726891, 24726964, 1, 'AT5G61445'), (24856205, 24856242, 1, 'AT5G61835'), (25347261, 25347333, 1, 'AT5G63145'), (25801340, 25801414, 1, 'AT5G64505'), (25892619, 25892691, -1, 'AT5G64735'), (25942291, 25942372, 1, 'AT5G64855'), (25989903, 25989976, 1, 'AT5G65015'), (26114755, 26114793, -1, 'AT5G65305'), (26174414, 26174496, -1, 'AT5G65445'), (26212684, 26212757, 1, 'AT5G65535'), (26238859, 26238933, -1, 'AT5G65615'), (26573248, 26573322, -1, 'AT5G66535'), (26585622, 26585696, 1, 'AT5G66568'), (26670495, 26670567, -1, 'AT5G66755'), (26699933, 26700004, -1, 'AT5G66817'), (26938897, 26938969, 1, 'AT5G67455')]
        entries = [("Chr I", "NC_003070", 30432563, f1, colors.red),
                   ("Chr II", "NC_003071", 19705359, f2, colors.green),
                   ("Chr III", "NC_003074", 23470805, f3, colors.blue),
                   ("Chr IV", "NC_003075", 18585042, f4, colors.orange),
                   ("Chr V", "NC_003076", 26992728, f5, colors.purple)]
        max_length = max([row[2] for row in entries])
 
        chr_diagram = BasicChromosome.Organism()
        chr_diagram._legend_height = 0
        for name, acc, length, features, color in entries:
            if False:
                #How I generated the values above... and tested passing in SeqFeatures
                filename = "/Users/pjcock/Documents/comp_genomics/seed/%s.gbk" % acc
                import os
                if not os.path.isfile(filename):
                    continue
                from Bio import SeqIO
                record = SeqIO.read(filename, "gb")
                assert length == len(record)
                features = [f for f in record.features if f.type=="tRNA"]
                print name
                print [(int(f.location.start), int(f.location.end), \
                        f.strand, f.qualifiers['locus_tag'][0]) \
                       for f in features]
                #Output was copy and pasted to the script, see above.
                #Continue test using SeqFeature objects!
                #To test colours from the qualifiers,
                for i,f in enumerate(features):
                    f.qualifiers['color'] = [str(i % 16)]
            else:
                features = [(start,end,strand,label,color) \
                            for (start,end,strand,label) in features]
            cur_chromosome = BasicChromosome.Chromosome(name)
            #Set the length, adding and extra 20 percent for the tolomeres:
            cur_chromosome.scale_num = max_length * 1.2
            #Add an opening telomere
            start = BasicChromosome.TelomereSegment()
            start.scale = 0.1 * max_length
            cur_chromosome.add(start)
            #Add a body - using bp as the scale length here.
            body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
            body.scale = length
            cur_chromosome.add(body)
            #Add a closing telomere
            end = BasicChromosome.TelomereSegment(inverted=True)
            end.scale = 0.1 * max_length
            cur_chromosome.add(end)
            #This chromosome is done
            chr_diagram.add(cur_chromosome)
        chr_diagram.draw("Graphics/tRNA_chrom.pdf", "Arabidopsis thaliana tRNA")

        
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
