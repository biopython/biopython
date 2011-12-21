#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for GenomeDiagram general functionality.
"""

##########
# IMPORTS

# Builtins
import os
import unittest
import math

# Do we have ReportLab?  Raise error if not present.
from Bio import MissingPythonDependencyError
try:
    from reportlab.lib import colors
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.lib.units import cm
except ImportError:
    raise MissingPythonDependencyError(
            "Install reportlab if you want to use Bio.Graphics.")

# Biopython core
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqUtils

# Bio.Graphics.GenomeDiagram
from Bio.Graphics.GenomeDiagram import FeatureSet, GraphSet, Track, Diagram
from Bio.Graphics.GenomeDiagram import CrossLink
#from Bio.Graphics.GenomeDiagram.Utilities import *

#Currently private, but we test them here:
from Bio.Graphics.GenomeDiagram._Graph import GraphData
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator

def fill_and_border(base_color, alpha=0.5):
    try:
        c = base_color.clone()
        c.alpha = alpha
        return c, base_color
    except AttributeError:
        #Old ReportLab, no transparency and/or no clone
        return base_color, base_color

###############################################################################
# Utility functions for graph plotting, originally in GenomeDiagram.Utilities #
# See Bug 2705 for discussion on where to put these functions in Biopython... #
###############################################################################
def apply_to_window(sequence, window_size, function, step=None):
    """ apply_to_window(sequence, window_size, function) -> [(int, float),(int, float),...]

        o sequence      Bio.Seq.Seq object

        o window_size   Int describing the length of sequence to consider

        o step          Int describing the step to take between windows
                        (default = window_size//2)

        o function      Method or function that accepts a Bio.Seq.Seq object
                        as its sole argument and returns a single value

        Returns a list of (position, value) tuples for fragments of the passed
        sequence of length window_size (stepped by step), calculated by the
        passed function.  Returned positions are the midpoint of each window.
    """
    seqlen = len(sequence)      # Total length of sequence to be used
    if step is None:    # No step specified, so use half window-width or 1 if larger
        step = max(window_size//2, 1)
    else:               # Use specified step, or 1 if greater
        step = max(step, 1)

    results = []    # Holds (position, value) results

    # Perform the passed function on as many windows as possible, short of
    # overrunning the sequence
    pos = 0
    while pos < seqlen-window_size+1:
        # Obtain sequence fragment
        start, middle, end = pos, (pos+window_size+pos)//2, pos+window_size
        fragment = sequence[start:end]
        # Apply function to the sequence fragment
        value = function(fragment)
        results.append((middle, value)) # Add results to list
        # Advance to next fragment
        pos += step

    # Use the last available window on the sequence, even if it means
    # re-covering old ground
    if pos != seqlen - window_size:
        # Obtain sequence fragment
        pos = seqlen - window_size
        start, middle, end = pos, (pos+window_size+pos)//2, pos+window_size
        fragment = sequence[start:end]
        # Apply function to sequence fragment
        value = function(fragment)
        results.append((middle, value)) # Add results to list
        
    # Check on last sequence
    #print fragment
    #print seq[-100:]
    return results      # Return the list of (position, value) results

def calc_gc_content(sequence):
    """ calc_gc_content(sequence)

        o sequence  A Bio.Seq.Seq object

        Returns the % G+C content in a passed sequence
    """
    d = {}
    for nt in ['A','T','G','C']:
        d[nt] = sequence.count(nt) + sequence.count(nt.lower())
    gc = d.get('G',0) + d.get('C',0)

    if gc == 0: return 0
    #print gc*100.0/(d['A'] +d['T'] + gc)
    return gc*1./(d['A'] +d['T'] + gc)


def calc_at_content(sequence):
    """ calc_at_content(sequence)

        o sequence  A Bio.Seq.Seq object

        Returns the % A+T content in a passed sequence
    """
    d = {}
    for nt in ['A','T','G','C']:
        d[nt] = sequence.count(nt) + sequence.count(nt.lower())
    at = d.get('A',0) + d.get('T',0)

    if at == 0: return 0
    return at*1./(d['G'] +d['G'] + at)


def calc_gc_skew(sequence):
    """ calc_gc_skew(sequence)

        o sequence   A Bio.Seq.Seq object

        Returns the (G-C)/(G+C) GC skew in a passed sequence
    """
    g = sequence.count('G') + sequence.count('g')
    c = sequence.count('C') + sequence.count('c')
    if g+c == 0:
        return 0.0 #TODO - return NaN or None here?
    else:
        return (g-c)/float(g+c)


def calc_at_skew(sequence):
    """ calc_at_skew(sequence)

        o sequence   A Bio.Seq.Seq object

        Returns the (A-T)/(A+T) AT skew in a passed sequence
    """
    a = sequence.count('A') + sequence.count('a')
    t = sequence.count('T') + sequence.count('t')
    if a+t == 0:
        return 0.0 #TODO - return NaN or None here?
    else:
        return (a-t)/float(a+t)

def calc_dinucleotide_counts(sequence):
    """Returns the total count of di-nucleotides repeats (e.g. "AA", "CC").

    This is purely for the sake of generating some non-random sequence
    based score for plotting, with no expected biological meaning.

    NOTE - Only considers same case pairs.
    NOTE - "AA" scores 1, "AAA" scores 2, "AAAA" scores 3 etc.
    """
    total = 0
    for letter in "ACTGUactgu":
        total += sequence.count(letter+letter)
    return total
    

###############################################################################
# End of utility functions for graph plotting                                 #
###############################################################################

# Tests
class TrackTest(unittest.TestCase):
    # TODO Bring code from Track.py, unsure about what test does
    pass

class ColorsTest(unittest.TestCase):
    def test_color_conversions(self):
        """Test color translations.
        """
        translator = ColorTranslator()
        
        # Does the translate method correctly convert the passed argument?
        assert translator.float1_color((0.5, 0.5, 0.5)) == translator.translate((0.5, 0.5, 0.5)), \
            "Did not correctly translate colour from floating point RGB tuple"
        assert translator.int255_color((1, 75, 240)) == translator.translate((1, 75, 240)), \
            "Did not correctly translate colour from integer RGB tuple"
        assert translator.artemis_color(7) == translator.translate(7), \
            "Did not correctly translate colour from Artemis colour scheme"                        
        assert translator.scheme_color(2) == translator.translate(2), \
            "Did not correctly translate colour from user-defined colour scheme"

            
class GraphTest(unittest.TestCase):
    def test_limits(self):
        """Check line graphs."""
        #TODO - Fix GD so that the same min/max is used for all three lines?
        points = 1000
        scale = math.pi * 2.0 / points
        data1 = [math.sin(x*scale) for x in range(points)]
        data2 = [math.cos(x*scale) for x in range(points)]
        data3 = [2*math.sin(2*x*scale) for x in range(points)]
        
        gdd = Diagram('Test Diagram', circular=False,
                      y=0.01, yt=0.01, yb=0.01,
                      x=0.01, xl=0.01, xr=0.01)
        gdt_data = gdd.new_track(1, greytrack=False)
        gds_data = gdt_data.new_set("graph")
        for data_values, name, color in zip([data1,data2,data3],
                                            ["sin", "cos", "2sin2"],
                                            ["red","green","blue"]):
            data = zip(range(points), data_values)
            gds_data.new_graph(data, "", style="line",
                               color = color, altcolor = color,
                               center = 0)

        gdd.draw(format='linear',
                 tracklines=False,
                 pagesize=(15*cm,15*cm),
                 fragments=1,
                 start=0, end=points)
        gdd.write(os.path.join('Graphics', "line_graph.pdf"), "pdf")
        #Circular diagram
        gdd.draw(tracklines=False,
                 pagesize=(15*cm,15*cm),
                 circular=True, #Data designed to be periodic
                 start=0, end=points, circle_core=0.5)
        gdd.write(os.path.join('Graphics', "line_graph_c.pdf"), "pdf")
        
    def test_slicing(self):
        """Check GraphData slicing."""
        gd = GraphData()
        gd.set_data([(1, 10), (5, 15), (20, 40)])
        gd.add_point((10, 20))
        
        assert gd[4:16] == [(5, 15), (10, 20)], \
                "Unable to insert and retrieve points correctly"


class LabelTest(unittest.TestCase):
    """Check label positioning."""
    def setUp(self):
        self.gdd = Diagram('Test Diagram', circular=False,
                           y=0.01, yt=0.01, yb=0.01,
                           x=0.01, xl=0.01, xr=0.01)

    def finish(self, name, circular=True):
        #And draw it...
        tracks = len(self.gdd.tracks)
        #Work arround the page orientation code being too clever
        #and flipping the h & w round:
        if tracks <= 3:
            orient = "landscape"
        else:
            orient = "portrait"
        self.gdd.draw(format='linear', orientation=orient,
                      tracklines=False,
                      pagesize=(15*cm,5*cm*tracks),
                      fragments=1,
                      start=0, end=400)
        self.gdd.write(os.path.join('Graphics', name+".pdf"), "pdf")
        #For the tutorial this might be useful:
        #self.gdd.write(os.path.join('Graphics', name+".png"), "png")
        if circular:
            #Circular diagram
            self.gdd.draw(tracklines=False,
                          pagesize=(15*cm,15*cm),
                          fragments=1,
                          circle_core=0.5,
                          start=0, end=400)
            self.gdd.write(os.path.join('Graphics', name+"_c.pdf"), "pdf")
    
    def add_track_with_sigils(self, **kwargs):
        self.gdt_features = self.gdd.new_track(1, greytrack=False)
        self.gds_features = self.gdt_features.new_set()
        for i in range(18):
            start = int((400 * i)/18.0)
            end = start + 17
            if i % 3 == 0:
                strand=None
                name = "Strandless"
                color=colors.orange
            elif i % 3 == 1:
                strand=+1
                name="Forward"
                color=colors.red
            else:
                strand = -1
                name="Reverse"
                color=colors.blue
            feature = SeqFeature(FeatureLocation(start, end), strand=strand)
            self.gds_features.add_feature(feature, name=name,
                                          color=color, label=True, **kwargs)

    def test_label_default(self):
        """Feature labels - default."""
        self.add_track_with_sigils()
        self.finish("labels_default")

class SigilsTest(unittest.TestCase):
    """Check the different feature sigils.

    These figures are intended to be used in the Tutorial..."""
    def setUp(self):
        self.gdd = Diagram('Test Diagram', circular=False,
                           y=0.01, yt=0.01, yb=0.01,
                           x=0.01, xl=0.01, xr=0.01)

    def add_track_with_sigils(self, **kwargs):
        #Add a track of features,
        self.gdt_features = self.gdd.new_track(1, greytrack=False)
        #We'll just use one feature set for these features,
        self.gds_features = self.gdt_features.new_set()
        #Add three features to show the strand options,
        feature = SeqFeature(FeatureLocation(25, 125), strand=+1)
        self.gds_features.add_feature(feature, name="Forward", **kwargs)
        feature = SeqFeature(FeatureLocation(150, 250), strand=None)
        self.gds_features.add_feature(feature, name="Strandless", **kwargs)
        feature = SeqFeature(FeatureLocation(275, 375), strand=-1)
        self.gds_features.add_feature(feature, name="Reverse", **kwargs)

    def finish(self, name, circular=True):
        #And draw it...
        tracks = len(self.gdd.tracks)
        #Work arround the page orientation code being too clever
        #and flipping the h & w round:
        if tracks <= 3:
            orient = "landscape"
        else:
            orient = "portrait"
        self.gdd.draw(format='linear', orientation=orient,
                      tracklines=False,
                      pagesize=(15*cm,5*cm*tracks),
                      fragments=1,
                      start=0, end=400)
        self.gdd.write(os.path.join('Graphics', name+".pdf"), "pdf")
        #For the tutorial this might be useful:
        #self.gdd.write(os.path.join('Graphics', name+".png"), "png")
        if circular:
            #Circular diagram
            self.gdd.draw(tracklines=False,
                          pagesize=(15*cm,15*cm),
                          fragments=1,
                          circle_core=0.5,
                          start=0, end=400)
            self.gdd.write(os.path.join('Graphics', name+"_c.pdf"), "pdf")

    def test_labels(self):
        """Feature labels."""
        self.add_track_with_sigils(label=True)
        self.add_track_with_sigils(label=True, color="green",
                                   label_size=25, label_angle=0)
        self.add_track_with_sigils(label=True, color="purple",
                                   label_position="end",
                                   label_size=4, label_angle=90)
        self.add_track_with_sigils(label=True, color="blue",
                                   label_position="middle",
                                   label_size=6, label_angle=-90)
        self.assertEqual(len(self.gdd.tracks), 4)
        self.finish("GD_sigil_labels", circular=False)

    def test_arrow_shafts(self):
        """Feature arrow sigils, varying shafts."""
        self.add_track_with_sigils(sigil="ARROW")
        self.add_track_with_sigils(sigil="ARROW", color="brown",
                                   arrowshaft_height=1.0)
        self.add_track_with_sigils(sigil="ARROW", color="teal",
                                   arrowshaft_height=0.2)
        self.add_track_with_sigils(sigil="ARROW", color="darkgreen",
                                   arrowshaft_height=0.1)
        self.assertEqual(len(self.gdd.tracks), 4)
        self.finish("GD_sigil_arrow_shafts")        

    def test_arrow_heads(self):
        """Feature arrow sigils, varying heads."""
        self.add_track_with_sigils(sigil="ARROW")
        self.add_track_with_sigils(sigil="ARROW", color="blue",
                                   arrowhead_length=0.25)
        self.add_track_with_sigils(sigil="ARROW", color="orange",
                                   arrowhead_length=1)
        self.add_track_with_sigils(sigil="ARROW", color="red",
                                   arrowhead_length=10000) #Triangles
        self.assertEqual(len(self.gdd.tracks), 4)
        self.finish("GD_sigil_arrows")

    def test_small_arrow_heads(self):
        """Feature arrow sigil heads within bounding box."""
        #Add a track of features, bigger height to emphasise any sigil errors
        self.gdt_features = self.gdd.new_track(1, greytrack=True, height=3)
        #We'll just use one feature set for these features,
        self.gds_features = self.gdt_features.new_set()
        #Green arrows just have small heads (meaning if there is a mitre
        #it will escape the bounding box).  Red arrows are small triangles.
        feature = SeqFeature(FeatureLocation(15, 30), strand=+1)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Forward", sigil="ARROW",
                                      arrowhead_length=0.05)
        feature = SeqFeature(FeatureLocation(55, 60), strand=+1)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Forward", sigil="ARROW",
                                      arrowhead_length=1000, color="red")
        feature = SeqFeature(FeatureLocation(75, 125), strand=+1)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Forward", sigil="ARROW",
                                      arrowhead_length=0.05)
        feature = SeqFeature(FeatureLocation(140, 155), strand=None)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Strandless", sigil="ARROW",
                                      arrowhead_length=0.05)
        feature = SeqFeature(FeatureLocation(180, 185), strand=None)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Strandless", sigil="ARROW",
                                      arrowhead_length=1000, color="red")
        feature = SeqFeature(FeatureLocation(200, 250), strand=None)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Strandless", sigil="ARROW",
                                      arrowhead_length=0.05)
        feature = SeqFeature(FeatureLocation(265, 280), strand=-1)
        self.gds_features.add_feature(feature, name="Reverse", sigil="ARROW",
                                      arrowhead_length=0.05)
        feature = SeqFeature(FeatureLocation(305, 310), strand=-1)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Reverse", sigil="ARROW",
                                      arrowhead_length=1000, color="red")
        feature = SeqFeature(FeatureLocation(325, 375), strand=-1)
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(feature, name="Reverse", sigil="ARROW",
                                      arrowhead_length=0.05)
        self.finish("GD_sigil_arrows_small")

    def test_long_arrow_heads(self):
        """Feature arrow sigil heads within bounding box."""
        #Add a track of features, bigger height to emphasise any sigil errors
        self.gdt_features = self.gdd.new_track(1, greytrack=True, height=3)
        #We'll just use one feature set for these features,
        self.gds_features = self.gdt_features.new_set()
        feature = SeqFeature(FeatureLocation(25, 375), strand=+1)
        self.gds_features.add_feature(feature, color="lightblue")
        self.gds_features.add_feature(feature, name="Forward", sigil="ARROW",
                                      color="blue", arrowhead_length=2.0)
        feature = SeqFeature(FeatureLocation(25, 375), strand=-1)
        self.gds_features.add_feature(feature, color="pink")
        self.gds_features.add_feature(feature, name="Reverse", sigil="ARROW",
                                      color="red", arrowhead_length=2.0)
        #Add another track of features, bigger height to emphasise any sigil errors
        self.gdt_features = self.gdd.new_track(1, greytrack=True, height=3)
        #We'll just use one feature set for these features,
        self.gds_features = self.gdt_features.new_set()
        feature = SeqFeature(FeatureLocation(25, 375), strand=None)
        self.gds_features.add_feature(feature, color="lightgreen")
        self.gds_features.add_feature(feature, name="Standless", sigil="ARROW",
                                      color="green", arrowhead_length=2.0)
        self.finish("GD_sigil_arrows_long")

class DiagramTest(unittest.TestCase):
    """Creating feature sets, graph sets, tracks etc individually for the diagram."""
    def setUp(self):
        """Test setup, just loads a GenBank file as a SeqRecord."""
        handle = open(os.path.join("GenBank","NC_005816.gb"), 'r')
        self.record = SeqIO.read(handle, "genbank")
        handle.close()

    def test_write_arguments(self):
        """Check how the write methods respond to output format arguments."""
        gdd = Diagram('Test Diagram')
        gdd.drawing = None #Hack - need the ReportLab drawing object to be created.
        filename = os.path.join("Graphics","error.txt")
        #We (now) allow valid formats in any case.
        for output in ["XXX","xxx",None,123,5.9]:
            try:
                gdd.write(filename, output)
                assert False, \
                       "Should have rejected %s as an output format" % output
            except ValueError, e:
                #Good!
                pass
            try:
                gdd.write_to_string(output)
                assert False, \
                       "Should have rejected %s as an output format" % output
            except ValueError, e:
                #Good!
                pass

    def test_partial_diagram(self):
        """construct and draw SVG and PDF for just part of a SeqRecord."""
        genbank_entry = self.record
        start = 6500
        end = 8750
        
        gdd = Diagram('Test Diagram',
                      #For the circular diagram we don't want a closed cirle:
                      circular=False,
                      )
        #Add a track of features,
        gdt_features = gdd.new_track(1, greytrack=True,
                                     name="CDS Features",
                                     scale_largetick_interval=1000,
                                     scale_smalltick_interval=100,
                                     scale_format = "SInt",
                                     greytrack_labels=False,
                                     height=0.5)
        #We'll just use one feature set for these features,
        gds_features = gdt_features.new_set()
        for feature in genbank_entry.features:
            if feature.type != "CDS":
                #We're going to ignore these.
                continue
            if feature.location.end.position < start:
                #Out of frame (too far left)
                continue
            if feature.location.start.position > end:
                #Out of frame (too far right)
                continue

            #This URL should work in SVG output from recent versions
            #of ReportLab.  You need ReportLab 2.4 or later
            try :
                url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi"+\
                      "?db=protein&id=%s" % feature.qualifiers["protein_id"][0]
            except KeyError :
                url = None
                
            #Note that I am using strings for color names, instead
            #of passing in color objects.  This should also work!
            if len(gds_features) % 2 == 0:
                color = "white" #for testing the automatic black border!
            else:
                color = "red"
            #Checking it can cope with the old UK spelling colour.
            #Also show the labels perpendicular to the track.
            gds_features.add_feature(feature, colour=color,
                                     url = url,
                                     sigil="ARROW",
                                     label_position = "start",
                                     label_size = 8,
                                     label_angle = 90,
                                     label=True)

        #And draw it...
        gdd.draw(format='linear', orientation='landscape',
                 tracklines=False, pagesize=(10*cm,6*cm), fragments=1,
                 start=start, end=end)
        output_filename = os.path.join('Graphics', 'GD_region_linear.pdf')
        gdd.write(output_filename, 'PDF')

        #Also check the write_to_string method matches,
        #(Note the possible confusion over new lines on Windows)
        assert open(output_filename).read().replace("\r\n","\n") \
               == gdd.write_to_string('PDF').replace("\r\n","\n")

        output_filename = os.path.join('Graphics', 'GD_region_linear.svg')
        gdd.write(output_filename, 'SVG')

        #Circular with a particular start/end is a bit odd, but by setting
        #circular=False (above) a sweep of 90% is used (a wedge is left out)
        gdd.draw(format='circular',
                 tracklines=False, pagesize=(10*cm,10*cm),
                 start=start, end=end)
        output_filename = os.path.join('Graphics', 'GD_region_circular.pdf')
        gdd.write(output_filename, 'PDF')
        output_filename = os.path.join('Graphics', 'GD_region_circular.svg')
        gdd.write(output_filename, 'SVG')

    def test_diagram_via_methods_pdf(self):
        """Construct and draw PDF using method approach."""
        genbank_entry = self.record
        gdd = Diagram('Test Diagram')

        #Add a track of features,
        gdt_features = gdd.new_track(1, greytrack=True,
                                     name="CDS Features", greytrack_labels=0,
                                     height=0.5)
        #We'll just use one feature set for the genes and misc_features,
        gds_features = gdt_features.new_set()
        for feature in genbank_entry.features:
            if feature.type == "gene":
                if len(gds_features) % 2 == 0:
                    color = "blue"
                else:
                    color = "lightblue"
                gds_features.add_feature(feature, color=color,
                                            #label_position = "middle",
                                            #label_position = "end",
                                            label_position = "start",
                                            label_size = 11,
                                            #label_angle = 90,
                                            sigil="ARROW",
                                            label=True)

        #I want to include some strandless features, so for an example
        #will use EcoRI recognition sites etc.
        for site, name, color in [("GAATTC","EcoRI","green"),
                                  ("CCCGGG","SmaI","orange"),
                                  ("AAGCTT","HindIII","red"),
                                  ("GGATCC","BamHI","purple")]:
            index = 0
            while True:
                index  = genbank_entry.seq.find(site, start=index)
                if index == -1 : break
                feature = SeqFeature(FeatureLocation(index, index+6), strand=None)

                #This URL should work in SVG output from recent versions
                #of ReportLab.  You need ReportLab 2.4 or later
                try :
                    url = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi"+\
                          "?db=protein&id=%s" % feature.qualifiers["protein_id"][0]
                except KeyError :
                    url = None

                gds_features.add_feature(feature, color = color,
                                         url = url,
                                         #label_position = "middle",
                                         label_size = 10,
                                         label_color = color,
                                         #label_angle = 90,
                                         name = name,
                                         label = True)
                index += len(site)
            del index

        #Now add a graph track...
        gdt_at_gc = gdd.new_track(2, greytrack=True,
                                  name="AT and GC content",
                                  greytrack_labels=True)
        gds_at_gc = gdt_at_gc.new_set(type="graph")

        step = len(genbank_entry)//200
        gds_at_gc.new_graph(apply_to_window(genbank_entry.seq, step, calc_gc_content, step),
                        'GC content', style='line', 
                        color=colors.lightgreen,
                        altcolor=colors.darkseagreen)
        gds_at_gc.new_graph(apply_to_window(genbank_entry.seq, step, calc_at_content, step),
                        'AT content', style='line', 
                        color=colors.orange,
                        altcolor=colors.red)
        
        #Finally draw it in both formats,
        gdd.draw(format='linear', orientation='landscape',
             tracklines=0, pagesize='A4', fragments=3)
        output_filename = os.path.join('Graphics', 'GD_by_meth_linear.pdf')
        gdd.write(output_filename, 'PDF')

        gdd.draw(format='circular', tracklines=False, circle_core=0.8,
                 pagesize=(20*cm,20*cm), circular=True)
        output_filename = os.path.join('Graphics', 'GD_by_meth_circular.pdf')
        gdd.write(output_filename, 'PDF')

    def test_diagram_via_object_pdf(self):
        """Construct and draw PDF using object approach."""
        genbank_entry = self.record
        gdd = Diagram('Test Diagram')

        gdt1 = Track('CDS features', greytrack=True,
                     scale_largetick_interval=1e4,
                     scale_smalltick_interval=1e3,
                     greytrack_labels=10,
                     greytrack_font_color="red",
                     scale_format = "SInt")
        gdt2 = Track('gene features', greytrack=1,
                   scale_largetick_interval=1e4)

        #First add some feature sets:
        gdfsA = FeatureSet(name='CDS backgrounds')
        gdfsB = FeatureSet(name='gene background')


        gdfs1 = FeatureSet(name='CDS features')
        gdfs2 = FeatureSet(name='gene features')
        gdfs3 = FeatureSet(name='misc_features')
        gdfs4 = FeatureSet(name='repeat regions')

        prev_gene = None
        cds_count = 0
        for feature in genbank_entry.features:
            if feature.type == 'CDS':
                cds_count += 1
                if prev_gene:
                    #Assuming it goes with this CDS!
                    if cds_count % 2 == 0:
                        dark, light = colors.peru, colors.tan
                    else:
                        dark, light = colors.burlywood, colors.bisque
                    #Background for CDS,
                    a = gdfsA.add_feature(SeqFeature(FeatureLocation(feature.location.start, feature.location.end, strand=0)),
                                         color=dark)
                    #Background for gene,
                    b = gdfsB.add_feature(SeqFeature(FeatureLocation(prev_gene.location.start, prev_gene.location.end, strand=0)),
                                          color=dark)
                    #Cross link,
                    gdd.cross_track_links.append(CrossLink(a, b, light, dark))
                    prev_gene = None
            if feature.type == 'gene':
                prev_gene = feature

        #Some cross links on the same linear diagram fragment,
        f, c = fill_and_border(colors.red)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(2220,2230)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(2200,2210)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c))

        f, c = fill_and_border(colors.blue)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(2150,2200)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(2220,2290)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c, flip=True))

        f, c = fill_and_border(colors.green)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(2250,2560)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(2300,2860)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c))

        #Some cross links where both parts are saddling the linear diagram fragment boundary,
        f, c = fill_and_border(colors.red)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(3155,3250)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(3130,3300)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c))
        #Nestled within that (drawn on top),
        f, c = fill_and_border(colors.blue)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(3160,3275)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(3180,3225)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c, flip=True))

        #Some cross links where two features are on either side of the linear diagram fragment boundary,
        f, c = fill_and_border(colors.green)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(6450,6550)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(6265,6365)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c))
        f, c = fill_and_border(colors.gold)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(6265,6365)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(6450,6550)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c))
        f, c = fill_and_border(colors.red)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(6275,6375)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(6430,6530)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c, flip=True))
        f, c = fill_and_border(colors.blue)
        a = gdfsA.add_feature(SeqFeature(FeatureLocation(6430,6530)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(FeatureLocation(6275,6375)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c, flip=True))


        cds_count = 0
        for feature in genbank_entry.features:
            if feature.type == 'CDS':
                cds_count += 1
                if cds_count % 2 == 0:
                    gdfs1.add_feature(feature, color=colors.pink, sigil="ARROW")
                else:
                    gdfs1.add_feature(feature, color=colors.red, sigil="ARROW")

            if feature.type == 'gene':
                #Note we set the colour of ALL the genes later on as a test,
                gdfs2.add_feature(feature, sigil="ARROW")

            if feature.type == 'misc_feature':
                gdfs3.add_feature(feature, color=colors.orange)

            if feature.type == 'repeat_region':
                gdfs4.add_feature(feature, color=colors.purple)

        #gdd.cross_track_links = gdd.cross_track_links[:1]

        gdfs1.set_all_features('label', 1)
        gdfs2.set_all_features('label', 1)
        gdfs3.set_all_features('label', 1)
        gdfs4.set_all_features('label', 1)

        gdfs3.set_all_features('hide', 0)
        gdfs4.set_all_features('hide', 0)

        #gdfs1.set_all_features('color', colors.red)
        gdfs2.set_all_features('color', colors.blue)

        gdt1.add_set(gdfsA) #Before CDS so under them!
        gdt1.add_set(gdfs1)

        gdt2.add_set(gdfsB) #Before genes so under them!
        gdt2.add_set(gdfs2)
                
        gdt3 = Track('misc features and repeats', greytrack=1,
                   scale_largetick_interval=1e4)
        gdt3.add_set(gdfs3)
        gdt3.add_set(gdfs4)

        #Now add some graph sets:

        #Use a fairly large step so we can easily tell the difference
        #between the bar and line graphs.
        step = len(genbank_entry)//200
        gdgs1 = GraphSet('GC skew')
        
        graphdata1 = apply_to_window(genbank_entry.seq, step, calc_gc_skew, step)
        gdgs1.new_graph(graphdata1, 'GC Skew', style='bar',
                color=colors.violet,
                altcolor=colors.purple)
        
        gdt4 = Track(\
                'GC Skew (bar)',
                height=1.94, greytrack=1,
                scale_largetick_interval=1e4)
        gdt4.add_set(gdgs1)


        gdgs2 = GraphSet('GC and AT Content')
        gdgs2.new_graph(apply_to_window(genbank_entry.seq, step, calc_gc_content, step),
                        'GC content', style='line', 
                        color=colors.lightgreen,
                        altcolor=colors.darkseagreen)

        gdgs2.new_graph(apply_to_window(genbank_entry.seq, step, calc_at_content, step),
                        'AT content', style='line', 
                        color=colors.orange,
                        altcolor=colors.red)    

        gdt5 = Track(\
                'GC Content(green line), AT Content(red line)',
                height=1.94, greytrack=1,
                scale_largetick_interval=1e4)
        gdt5.add_set(gdgs2)

        gdgs3 = GraphSet('Di-nucleotide count')
        step = len(genbank_entry)//400 #smaller step
        gdgs3.new_graph(apply_to_window(genbank_entry.seq, step, calc_dinucleotide_counts, step),
                        'Di-nucleotide count', style='heat', 
                        color=colors.red, altcolor=colors.orange)
        gdt6 = Track('Di-nucleotide count', height=0.5, greytrack=False, scale=False)
        gdt6.add_set(gdgs3)

        #Add the tracks (from both features and graphs)
        #Leave some white space in the middle/bottom
        gdd.add_track(gdt4, 3) # GC skew
        gdd.add_track(gdt5, 4) # GC and AT content
        gdd.add_track(gdt1, 5) # CDS features
        gdd.add_track(gdt2, 6) # Gene features
        gdd.add_track(gdt3, 7) # Misc features and repeat feature
        gdd.add_track(gdt6, 8) # Feature depth

        #Finally draw it in both formats, and full view and partial
        gdd.draw(format='circular', orientation='landscape',
             tracklines=0, pagesize='A0')
        output_filename = os.path.join('Graphics', 'GD_by_obj_circular.pdf')
        gdd.write(output_filename, 'PDF')

        gdd.circular=False
        gdd.draw(format='circular', orientation='landscape',
             tracklines=0, pagesize='A0', start=3000, end=6300)
        output_filename = os.path.join('Graphics', 'GD_by_obj_frag_circular.pdf')
        gdd.write(output_filename, 'PDF')

        gdd.draw(format='linear', orientation='landscape',
             tracklines=0, pagesize='A0', fragments=3)
        output_filename = os.path.join('Graphics', 'GD_by_obj_linear.pdf')
        gdd.write(output_filename, 'PDF')

        gdd.set_all_tracks("greytrack_labels", 2)
        gdd.draw(format='linear', orientation='landscape',
             tracklines=0, pagesize=(30*cm,10*cm), fragments=1,
             start=3000, end=6300)
        output_filename = os.path.join('Graphics', 'GD_by_obj_frag_linear.pdf')
        gdd.write(output_filename, 'PDF')

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
