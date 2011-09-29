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
        f1 = [(111889, 111961, -1), (306383, 306456, 1), (309274, 309347, -1), (515493, 515566, 1), (552639, 552711, 1), (604401, 604474, 1), (877648, 877720, 1), (892513, 892585, 1), (909809, 909882, -1), (1159021, 1159092, 1), (1324921, 1324959, 1), (1583770, 1583844, -1), (1817398, 1817470, 1), (1978082, 1978156, 1), (2025354, 2025427, 1), (2107396, 2107467, -1), (2111146, 2111217, -1), (2177883, 2177957, 1), (2334818, 2334891, 1), (2406830, 2406902, -1), (2588521, 2588593, 1), (2846538, 2846611, -1), (2879305, 2879377, 1), (2939418, 2939490, 1), (3431185, 3431257, -1), (3676606, 3676644, 1), (3678774, 3678848, -1), (3881528, 3881608, 1), (3914628, 3914700, -1), (4266985, 4267059, -1), (4285884, 4285956, -1), (4440211, 4440284, 1), (4522705, 4522779, -1), (4709631, 4709703, 1), (4741995, 4742068, 1), (4743091, 4743164, 1), (5189681, 5189755, -1), (5309641, 5309713, -1), (5380901, 5380983, 1), (5518055, 5518128, -1), (5619464, 5619537, -1), (6038749, 6038831, 1), (6075812, 6075884, 1), (6075937, 6076011, -1), (6345756, 6345828, 1), (6488645, 6488726, 1), (6948850, 6948934, -1), (6995272, 6995344, -1), (7004504, 7004576, 1), (7016506, 7016579, 1), (7082657, 7082729, 1), (7242749, 7242821, -1), (7499721, 7499793, -1), (7656108, 7656180, -1), (7884405, 7884443, -1), (8520278, 8520352, -1), (9143796, 9143870, 1), (9158169, 9158242, 1), (10089422, 10089494, 1), (10089883, 10089955, 1), (10090353, 10090425, 1), (10090754, 10090826, 1), (10092310, 10092382, 1), (10092786, 10092858, 1), (10093294, 10093366, 1), (10093731, 10093803, 1), (10094158, 10094230, 1), (10096936, 10097008, 1), (10097099, 10097171, 1), (10097703, 10097775, 1), (10098638, 10098710, 1), (10099064, 10099136, 1), (10099410, 10099482, 1), (10099812, 10099884, 1), (10100258, 10100330, 1), (10101013, 10101085, 1), (10101585, 10101657, 1), (10101978, 10102050, 1), (10106075, 10106147, 1), (10106513, 10106585, 1), (10106883, 10106955, 1), (10107634, 10107706, 1), (10108374, 10108446, 1), (10108695, 10108767, 1), (10207291, 10207364, -1), (10756703, 10756776, 1), (10963553, 10963627, -1), (11104093, 11104167, 1), (11797227, 11797265, -1), (12097258, 12097327, -1), (13687637, 13687710, 1), (15733055, 15733127, -1), (16588144, 16588216, -1), (17159046, 17159118, 1), (17159799, 17159871, 1), (17160970, 17161042, 1), (17161418, 17161490, 1), (17162967, 17163039, 1), (17163408, 17163480, 1), (17164461, 17164533, 1), (17735509, 17735582, 1), (18139265, 18139337, -1), (18234146, 18234220, -1), (18312570, 18312607, 1), (18391469, 18391542, 1), (18556666, 18556746, 1), (18561567, 18561647, 1), (19428223, 19428297, 1), (19502087, 19502161, -1), (19688850, 19688887, -1), (19851640, 19851714, 1), (19929506, 19929578, -1), (20416594, 20416667, -1), (20794976, 20795058, 1), (21272451, 21272533, 1), (21272786, 21272823, 1), (21273216, 21273253, 1), (21273960, 21274042, 1), (21274295, 21274332, 1), (21274725, 21274762, 1), (21275469, 21275551, 1), (21275804, 21275841, 1), (21276234, 21276271, 1), (21276978, 21277060, 1), (21277313, 21277350, 1), (21277743, 21277780, 1), (21278487, 21278569, 1), (21278822, 21278859, 1), (21279273, 21279310, 1), (21280016, 21280098, 1), (21280351, 21280388, 1), (21280781, 21280818, 1), (21281525, 21281607, 1), (21281860, 21281897, 1), (21282311, 21282348, 1), (21283054, 21283136, 1), (21283384, 21283421, 1), (21283842, 21283879, 1), (21284586, 21284668, 1), (21284916, 21284953, 1), (21285374, 21285411, 1), (21286118, 21286200, 1), (21286448, 21286485, 1), (21286906, 21286943, 1), (21287650, 21287732, 1), (21287980, 21288017, 1), (21288438, 21288475, 1), (21289183, 21289265, 1), (21289513, 21289550, 1), (21289970, 21290007, 1), (21290714, 21290796, 1), (21291044, 21291081, 1), (21291501, 21291538, 1), (21292245, 21292327, 1), (21292574, 21292611, 1), (21293032, 21293069, 1), (21293776, 21293858, 1), (21294109, 21294146, 1), (21294567, 21294604, 1), (21295125, 21295207, 1), (21295455, 21295492, 1), (21295912, 21295949, 1), (21296656, 21296738, 1), (21296989, 21297026, 1), (21297447, 21297484, 1), (21298005, 21298087, 1), (21298335, 21298372, 1), (21298792, 21298829, 1), (21299536, 21299618, 1), (21299869, 21299906, 1), (21300327, 21300364, 1), (21300885, 21300967, 1), (21301215, 21301252, 1), (21301673, 21301710, 1), (21302417, 21302499, 1), (21302750, 21302787, 1), (21303208, 21303245, 1), (21303766, 21303848, 1), (21304096, 21304133, 1), (21304554, 21304591, 1), (21305298, 21305380, 1), (21305631, 21305668, 1), (21306089, 21306126, 1), (21306647, 21306729, 1), (21306981, 21307018, 1), (21307441, 21307478, 1), (21308184, 21308268, 1), (21308520, 21308557, 1), (21308975, 21309012, 1), (21309719, 21309801, 1), (21310053, 21310090, 1), (21310513, 21310550, 1), (21311256, 21311340, 1), (21311592, 21311629, 1), (21312051, 21312088, 1), (21377983, 21378054, -1), (21887507, 21887589, -1), (22044276, 22044348, -1), (22317078, 22317149, -1), (22398301, 22398372, -1), (22401256, 22401327, -1), (22431831, 22431902, 1), (22481437, 22481511, -1), (22870422, 22870494, -1), (22890754, 22890834, 1), (23562849, 23562921, -1), (23671147, 23671219, -1), (23806215, 23806299, 1), (23936799, 23936872, 1), (24490654, 24490736, -1), (25833316, 25833388, 1), (25890198, 25890272, 1), (25931858, 25931931, 1), (25935739, 25935812, -1), (25944826, 25944898, 1), (25993392, 25993466, 1), (26053140, 26053214, 1), (26385816, 26385888, -1), (26977050, 26977121, 1), (27397046, 27397128, 1), (27792643, 27792715, 1), (28024043, 28024124, -1), (28031620, 28031701, 1), (28188192, 28188264, 1), (28377149, 28377222, -1), (28411644, 28411717, 1), (28444549, 28444621, 1), (28523645, 28523717, -1), (28531427, 28531499, 1), (28639585, 28639667, 1), (28952447, 28952519, -1), (29007098, 29007180, -1), (29147983, 29148055, -1), (29448865, 29448903, -1), (29809015, 29809088, 1), (29838009, 29838081, 1), (29838610, 29838682, 1), (30088888, 30088962, -1), (30178905, 30178977, -1), (30242675, 30242757, 1,)]
        f2 = [(102063, 102137, 1), (706794, 706867, 1), (846853, 846926, -1), (1054714, 1054787, -1), (1113980, 1114052, -1), (1123386, 1123458, -1), (1154381, 1154454, 1), (3239653, 3239725, -1), (3255828, 3255902, -1), (3268803, 3268883, 1), (3276436, 3276508, 1), (3280859, 3280933, 1), (3290962, 3291034, 1), (3303240, 3303312, -1), (3303350, 3303425, -1), (3303781, 3303819, -1), (3328666, 3328739, -1), (3332674, 3332756, 1), (3369350, 3369437, 1), (3383400, 3383474, -1), (3444359, 3444431, -1), (3452973, 3453060, 1), (3462074, 3462148, 1), (3494378, 3494416, 1), (3494772, 3494847, 1), (3495008, 3495083, 1), (3495438, 3495509, 1), (3496436, 3496508, 1), (3497354, 3497437, 1), (3503518, 3503605, 1), (6953924, 6953961, -1), (7046175, 7046247, 1), (7749793, 7749867, 1), (7962758, 7962832, -1), (9144435, 9144507, 1), (9241319, 9241356, -1), (9273888, 9273969, -1), (9277742, 9277814, -1), (9291113, 9291185, 1), (9400749, 9400823, 1), (9456888, 9456962, -1), (9472660, 9472733, -1), (9509359, 9509433, 1), (9598106, 9598179, 1), (9810296, 9810368, -1), (10066525, 10066597, -1), (10380655, 10380728, 1), (10820917, 10820990, 1), (11122756, 11122837, -1), (11781928, 11782000, -1), (11871230, 11871302, -1), (12336079, 12336151, 1), (12346827, 12346899, 1), (12478849, 12478921, -1), (12645232, 12645305, -1), (12888667, 12888738, 1), (12889810, 12889881, 1), (12983024, 12983095, -1), (13144312, 13144385, -1), (13658350, 13658425, 1), (14054465, 14054503, -1), (14250206, 14250278, 1), (14251774, 14251846, 1), (14357464, 14357536, 1), (14358437, 14358509, 1), (14359269, 14359341, 1), (14360221, 14360293, 1), (14360734, 14360806, 1), (14361176, 14361248, 1), (14362215, 14362287, 1), (14363133, 14363205, 1), (14363599, 14363671, 1), (14750553, 14750627, -1), (14757142, 14757213, 1), (14847685, 14847723, 1), (15175940, 15176014, 1), (15176656, 15176736, 1), (15215480, 15215517, -1), (15327312, 15327395, 1), (15327463, 15327546, -1), (15353238, 15353311, 1), (15477287, 15477324, -1), (15923894, 15923967, 1), (16525641, 16525713, -1), (16525846, 16525918, 1), (16646857, 16646929, -1), (17545780, 17545862, -1), (17667855, 17667926, 1), (17880766, 17880839, 1), (18002649, 18002721, -1), (18317052, 18317134, -1), (18576985, 18577058, 1), (18710751, 18710824, 1), (18963713, 18963786, 1), (19351496, 19351569, 1), (19566924, 19566995, -1)]
        f3 = [(259640, 259712, 1), (469666, 469740, 1), (476808, 476880, 1), (586092, 586174, 1), (981975, 982047, 1), (984105, 984177, 1), (1220234, 1220307, 1), (1601343, 1601415, -1), (1707743, 1707815, -1), (1738796, 1738870, 1), (1843329, 1843400, -1), (1920038, 1920110, -1), (2104961, 2105033, -1), (2222251, 2222324, 1), (2232470, 2232506, -1), (2253680, 2253762, -1), (2285607, 2285679, 1), (2918418, 2918492, -1), (2944616, 2944698, 1), (2945700, 2945782, -1), (3090548, 3090631, 1), (3096220, 3096293, 1), (3238371, 3238407, -1), (3535151, 3535224, 1), (3575849, 3575923, 1), (3622697, 3622769, -1), (3942012, 3942084, 1), (3995103, 3995176, -1), (4254534, 4254615, 1), (4330778, 4330850, 1), (4998147, 4998219, 1), (5068300, 5068374, -1), (5275155, 5275228, 1), (5632857, 5632930, 1), (6483945, 6484019, -1), (6540636, 6540673, 1), (6663713, 6663786, 1), (7104314, 7104398, 1), (7224223, 7224296, -1), (7319582, 7319664, -1), (7567399, 7567471, -1), (9373610, 9373684, -1), (9840420, 9840494, 1), (10211564, 10211636, 1), (10319498, 10319570, 1), (10325875, 10325947, 1), (10753667, 10753740, 1), (10760629, 10760702, -1), (11076814, 11076886, 1), (11961645, 11961718, 1), (16438025, 16438097, -1), (16896875, 16896949, 1), (16902623, 16902697, 1), (16905147, 16905221, 1), (17160736, 17160808, 1), (17275564, 17275646, 1), (17905395, 17905467, 1), (17985575, 17985611, -1), (18080062, 18080134, 1), (18518796, 18518870, 1), (18755788, 18755860, -1), (18837020, 18837092, 1), (18907851, 18907924, 1), (18928413, 18928487, 1), (19008621, 19008694, -1), (19044371, 19044443, -1), (19403651, 19403723, -1), (19420345, 19420417, -1), (19511965, 19512045, 1), (19566013, 19566085, 1), (19648105, 19648188, 1), (19935354, 19935426, 1), (19995918, 19995989, 1), (20704664, 20704736, 1), (20720151, 20720223, 1), (20824495, 20824568, -1), (21498293, 21498375, 1), (21553258, 21553329, 1), (21970486, 21970557, 1), (22149699, 22149773, 1), (22149823, 22149895, -1), (22197810, 22197892, -1), (22481215, 22481288, -1), (22622384, 22622465, 1), (22786896, 22786969, 1), (22853496, 22853567, 1), (22871101, 22871174, 1), (22892781, 22892853, 1), (23047854, 23047927, 1), (23062444, 23062517, -1), (23221682, 23221753, 1), (23296567, 23296640, -1), (23296728, 23296801, -1)]
        f4 = [(33799, 33872, 1), (424716, 424788, -1), (562560, 562634, -1), (611865, 611932, -1), (808269, 808342, -1), (901175, 901247, 1), (1390894, 1390966, 1), (1442004, 1442076, 1), (1501605, 1501677, 1), (1520781, 1520854, -1), (5268124, 5268210, -1), (6646425, 6646496, 1), (6819287, 6819324, 1), (6837555, 6837639, -1), (6837769, 6837853, -1), (6905479, 6905552, -1), (6944721, 6944793, 1), (7185697, 7185771, 1), (7232792, 7232865, -1), (7256408, 7256481, 1), (7341420, 7341494, -1), (7730956, 7731037, 1), (7814197, 7814270, 1), (8255695, 8255767, 1), (8301720, 8301794, -1), (8979656, 8979729, 1), (9108317, 9108391, 1), (9191590, 9191663, 1), (9287230, 9287304, 1), (9289706, 9289787, 1), (9815215, 9815287, -1), (9873524, 9873596, -1), (9978117, 9978189, -1), (10093077, 10093157, -1), (10302011, 10302084, 1), (10325975, 10326047, -1), (10878733, 10878807, -1), (11774472, 11774508, -1), (11910299, 11910373, 1), (11954751, 11954824, -1), (11974951, 11975032, 1), (12320119, 12320203, 1), (12429608, 12429681, 1), (12486211, 12486282, -1), (12686148, 12686230, 1), (13006243, 13006316, -1), (13058840, 13058922, -1), (13076582, 13076666, -1), (13285431, 13285503, -1), (13336345, 13336419, -1), (13341501, 13341575, -1), (13454562, 13454635, 1), (13704787, 13704860, 1), (13882922, 13882994, -1), (13885196, 13885269, -1), (14032495, 14032567, 1), (14267286, 14267368, 1), (14470283, 14470355, 1), (15120655, 15120728, 1), (15183089, 15183162, 1), (15345717, 15345753, -1), (15430229, 15430303, -1), (15576655, 15576728, 1), (15671398, 15671469, 1), (15804553, 15804635, 1), (16304128, 16304201, 1), (16454700, 16454773, -1), (16556627, 16556700, 1), (16655290, 16655364, 1), (17130054, 17130127, 1), (17149473, 17149545, 1), (17276705, 17276779, -1), (17500800, 17500872, -1), (18254982, 18255018, -1), (18293773, 18293845, 1), (18395021, 18395093, 1), (18411258, 18411332, 1), (18501705, 18501778, -1), (18542164, 18542238, 1)]
        f5 = [(150353, 150426, -1), (389889, 389960, -1), (508427, 508500, -1), (530819, 530893, 1), (559327, 559399, -1), (588890, 588964, -1), (614641, 614723, 1), (642397, 642479, -1), (858534, 858571, 1), (862395, 862468, -1), (970797, 970878, -1), (984365, 984448, 1), (998940, 999013, 1), (1742692, 1742765, 1), (1788651, 1788723, 1), (1804616, 1804690, 1), (1853302, 1853382, -1), (2060153, 2060235, -1), (2212678, 2212749, -1), (2309512, 2309549, -1), (2411148, 2411232, 1), (2432263, 2432336, -1), (2587826, 2587899, -1), (2898867, 2898951, -1), (2993327, 2993401, 1), (3030817, 3030890, -1), (3118377, 3118458, 1), (3212351, 3212424, -1), (3287553, 3287635, -1), (3324702, 3324775, 1), (3578295, 3578367, -1), (3617058, 3617130, 1), (3669000, 3669073, -1), (4471050, 4471122, 1), (4530475, 4530548, 1), (4673902, 4673974, 1), (4929562, 4929636, 1), (5157641, 5157715, 1), (5161514, 5161586, 1), (5358918, 5359000, 1), (5962699, 5962771, -1), (5965972, 5966044, -1), (5984378, 5984450, 1), (6258146, 6258218, 1), (6401240, 6401311, 1), (7073531, 7073603, -1), (7073944, 7074016, -1), (7074357, 7074429, -1), (7074773, 7074845, -1), (7222059, 7222131, -1), (7387890, 7387962, 1), (7981400, 7981472, 1), (8906418, 8906502, 1), (8946826, 8946899, -1), (9815405, 9815477, -1), (11802284, 11802356, 1), (13823211, 13823284, -1), (15049737, 15049811, -1), (15242547, 15242621, 1), (15593086, 15593160, 1), (15844253, 15844325, -1), (15993514, 15993587, 1), (16256865, 16256937, -1), (16427812, 16427893, 1), (16524760, 16524832, -1), (16655393, 16655477, 1), (16684663, 16684735, -1), (17476402, 17476475, -1), (17512768, 17512839, -1), (17856811, 17856883, -1), (17894906, 17894979, -1), (18058014, 18058088, 1), (18560206, 18560278, -1), (18576071, 18576143, 1), (18715888, 18715960, -1), (18807534, 18807614, 1), (18924749, 18924821, 1), (19658828, 19658900, 1), (19761400, 19761472, -1), (19820360, 19820398, 1), (20064048, 20064120, 1), (20692447, 20692519, 1), (20758903, 20758940, -1), (20773555, 20773637, 1), (21275059, 21275141, -1), (21318105, 21318189, -1), (21418369, 21418441, 1), (21740339, 21740410, -1), (22091631, 22091704, 1), (22094087, 22094160, 1), (22304851, 22304923, -1), (22355897, 22355970, -1), (22357726, 22357799, -1), (22501995, 22502068, -1), (22845356, 22845430, 1), (22973066, 22973138, 1), (23071996, 23072070, -1), (23463219, 23463291, 1), (23661936, 23662018, 1), (23861431, 23861503, 1), (23971167, 23971239, 1), (23974655, 23974727, 1), (24157171, 24157245, -1), (24279805, 24279886, 1), (24547401, 24547474, 1), (24548892, 24548964, 1), (24684507, 24684579, 1), (24726891, 24726964, 1), (24856205, 24856242, 1), (25347261, 25347333, 1), (25801340, 25801414, 1), (25892619, 25892691, -1), (25942291, 25942372, 1), (25989903, 25989976, 1), (26114755, 26114793, -1), (26174414, 26174496, -1), (26212684, 26212757, 1), (26238859, 26238933, -1), (26573248, 26573322, -1), (26585622, 26585696, 1), (26670495, 26670567, -1), (26699933, 26700004, -1), (26938897, 26938969, 1)]
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
                print [(int(f.location.start), int(f.location.end), f.strand) for f in features]
                #Output was copy and pasted to the script, see above.
                #Continue test using SeqFeature objects!
                #To test colours from the qualifiers,
                for i,f in enumerate(features):
                    f.qualifiers['color'] = [str(i % 16)]
            else:
                features = [(start,end,strand,color) for (start,end,strand) in features]
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
