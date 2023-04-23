# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for GenomeDiagram general functionality."""


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
        "Install reportlab if you want to use Bio.Graphics."
    ) from None

try:
    # The preferred PIL import has changed over time...
    try:
        from PIL import Image
    except ImportError:
        import Image
    from reportlab.graphics import renderPM
except ImportError:
    # This is an optional part of ReportLab, so may not be installed.
    # We'll raise a missing dependency error if rendering to a
    # bitmap format is attempted.
    renderPM = None

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation

from Bio.Graphics.GenomeDiagram import FeatureSet, GraphSet, Track, Diagram
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.Graphics.GenomeDiagram._Graph import GraphData
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator

from reportlab import rl_config

rl_config.invariant = True


def fill_and_border(base_color, alpha=0.5):
    """Return fill and border colors given a base color."""
    try:
        c = base_color.clone()
        c.alpha = alpha
        return c, base_color
    except AttributeError:
        # Old ReportLab, no transparency and/or no clone
        return base_color, base_color


###############################################################################
# Utility functions for graph plotting, originally in GenomeDiagram.Utilities #
# See Bug 2705 for discussion on where to put these functions in Biopython... #
###############################################################################


def apply_to_window(sequence, window_size, function, step=None):
    """Apply function to windows of the given sequence.

    Returns a list of (position, value) tuples for fragments of the passed
    sequence of length window_size (stepped by step), calculated by the passed
    function.  Returned positions are the midpoint of each window.

    - sequence - Bio.Seq.Seq object.
    - window_size - an integer describing the length of sequence to consider.
    - step - an integer describing the step to take between windows
      (default = window_size//2).
    - function - Method or function that accepts a Bio.Seq.Seq object
      as its sole argument and returns a single value.
    """
    seqlen = len(sequence)  # Total length of sequence to be used
    if step is None:  # No step specified, so use half window-width or 1 if larger
        step = max(window_size // 2, 1)
    else:  # Use specified step, or 1 if greater
        step = max(step, 1)

    results = []  # Holds (position, value) results

    # Perform the passed function on as many windows as possible, short of
    # overrunning the sequence
    pos = 0
    while pos < seqlen - window_size + 1:
        # Obtain sequence fragment
        start, middle, end = pos, (pos + window_size + pos) // 2, pos + window_size
        fragment = sequence[start:end]
        # Apply function to the sequence fragment
        value = function(fragment)
        results.append((middle, value))  # Add results to list
        # Advance to next fragment
        pos += step

    # Use the last available window on the sequence, even if it means
    # re-covering old ground
    if pos != seqlen - window_size:
        # Obtain sequence fragment
        pos = seqlen - window_size
        start, middle, end = pos, (pos + window_size + pos) // 2, pos + window_size
        fragment = sequence[start:end]
        # Apply function to sequence fragment
        value = function(fragment)
        results.append((middle, value))  # Add results to list

    return results  # Return the list of (position, value) results


def calc_gc_content(sequence):
    """Return the % G+C content in a passed sequence.

    Arguments:
        - sequence  - a Bio.Seq.Seq object.

    calc_gc_content(sequence)

    """
    d = {}
    for nt in ["A", "T", "G", "C"]:
        d[nt] = sequence.count(nt) + sequence.count(nt.lower())
    gc = d.get("G", 0) + d.get("C", 0)

    if gc == 0:
        return 0
    return gc / (d["A"] + d["T"] + gc)


def calc_at_content(sequence):
    """Return the % A+T content in a passed sequence.

    Arguments:
        - sequence  - a Bio.Seq.Seq object.

    calc_at_content(sequence)

    """
    d = {}
    for nt in ["A", "T", "G", "C"]:
        d[nt] = sequence.count(nt) + sequence.count(nt.lower())
    at = d.get("A", 0) + d.get("T", 0)

    if at == 0:
        return 0
    return at / (d["G"] + d["G"] + at)


def calc_gc_skew(sequence):
    """Return the (G-C)/(G+C) GC skew in a passed sequence.

    Arguments:
        - sequence   - a Bio.Seq.Seq object.

    calc_gc_skew(sequence)

    """
    g = sequence.count("G") + sequence.count("g")
    c = sequence.count("C") + sequence.count("c")
    if g + c == 0:
        return 0.0  # TODO - return NaN or None here?
    else:
        return (g - c) / (g + c)


def calc_at_skew(sequence):
    """Return the (A-T)/(A+T) AT skew in a passed sequence.

    Arguments:
        - sequence   - a Bio.Seq.Seq object.

    calc_at_skew(sequence)

    """
    a = sequence.count("A") + sequence.count("a")
    t = sequence.count("T") + sequence.count("t")
    if a + t == 0:
        return 0.0  # TODO - return NaN or None here?
    else:
        return (a - t) / (a + t)


def calc_dinucleotide_counts(sequence):
    """Return the total count of di-nucleotides repeats (e.g. "AA", "CC").

    This is purely for the sake of generating some non-random sequence
    based score for plotting, with no expected biological meaning.

    NOTE - Only considers same case pairs.
    NOTE - "AA" scores 1, "AAA" scores 2, "AAAA" scores 3 etc.
    """
    total = 0
    for letter in "ACTGUactgu":
        total += sequence.count(letter + letter)
    return total


###############################################################################
# End of utility functions for graph plotting                                 #
###############################################################################


# Tests
# class TrackTest(unittest.TestCase):
#    # TODO Bring code from Track.py, unsure about what test does
#    pass


class ColorsTest(unittest.TestCase):
    """Color tests."""

    def test_color_conversions(self):
        """Test color translations."""
        translator = ColorTranslator()

        # Does the translate method correctly convert the passed argument?
        self.assertEqual(
            translator.float1_color((0.5, 0.5, 0.5)),
            translator.translate((0.5, 0.5, 0.5)),
            "Did not correctly translate colour from floating point RGB tuple",
        )
        self.assertEqual(
            translator.int255_color((1, 75, 240)),
            translator.translate((1, 75, 240)),
            "Did not correctly translate colour from integer RGB tuple",
        )
        self.assertEqual(
            translator.artemis_color(7),
            translator.translate(7),
            "Did not correctly translate colour from Artemis colour scheme",
        )
        self.assertEqual(
            translator.scheme_color(2),
            translator.translate(2),
            "Did not correctly translate colour from user-defined colour scheme",
        )


class GraphTest(unittest.TestCase):
    """Graph tests."""

    def test_limits(self):
        """Check line graphs."""
        # TODO - Fix GD so that the same min/max is used for all three lines?
        points = 1000
        scale = math.pi * 2.0 / points
        data1 = [math.sin(x * scale) for x in range(points)]
        data2 = [math.cos(x * scale) for x in range(points)]
        data3 = [2 * math.sin(2 * x * scale) for x in range(points)]

        gdd = Diagram(
            "Test Diagram",
            circular=False,
            y=0.01,
            yt=0.01,
            yb=0.01,
            x=0.01,
            xl=0.01,
            xr=0.01,
        )
        gdt_data = gdd.new_track(1, greytrack=False)
        gds_data = gdt_data.new_set("graph")
        for data_values, _name, color in zip(
            [data1, data2, data3], ["sin", "cos", "2sin2"], ["red", "green", "blue"]
        ):
            data = list(zip(range(points), data_values))
            gds_data.new_graph(
                data, "", style="line", color=color, altcolor=color, center=0
            )

        gdd.draw(
            format="linear",
            tracklines=False,
            pagesize=(15 * cm, 15 * cm),
            fragments=1,
            start=0,
            end=points,
        )
        gdd.write(os.path.join("Graphics", "line_graph.pdf"), "pdf")
        # Circular diagram
        gdd.draw(
            tracklines=False,
            pagesize=(15 * cm, 15 * cm),
            circular=True,  # Data designed to be periodic
            start=0,
            end=points,
            circle_core=0.5,
        )
        gdd.write(os.path.join("Graphics", "line_graph_c.pdf"), "pdf")

    def test_slicing(self):
        """Check GraphData slicing."""
        gd = GraphData()
        gd.set_data([(1, 10), (5, 15), (20, 40)])
        gd.add_point((10, 20))

        self.assertEqual(
            gd[4:16],
            [(5, 15), (10, 20)],  # noqa 231
            "Unable to insert and retrieve points correctly",
        )


class LabelTest(unittest.TestCase):
    """Check label positioning."""

    def setUp(self):
        """Start a diagram."""
        self.gdd = Diagram(
            "Test Diagram",
            circular=False,
            y=0.01,
            yt=0.01,
            yb=0.01,
            x=0.01,
            xl=0.01,
            xr=0.01,
        )

    def finish(self, name, circular=True):
        """Draw it..."""
        tracks = len(self.gdd.tracks)
        # Work around the page orientation code being too clever
        # and flipping the h & w round:
        if tracks <= 3:
            orient = "landscape"
        else:
            orient = "portrait"
        self.gdd.draw(
            format="linear",
            orientation=orient,
            tracklines=False,
            pagesize=(15 * cm, 5 * cm * tracks),
            fragments=1,
            start=0,
            end=400,
        )
        self.gdd.write(os.path.join("Graphics", name + ".pdf"), "pdf")
        global renderPM
        if renderPM:
            try:
                # For the tutorial this is useful:
                self.gdd.write(os.path.join("Graphics", name + ".png"), "png")
            except renderPM.RenderPMError:
                # Probably a font problem, e.g.
                # RenderPMError: Can't setFont(Times-Roman) missing the T1 files?
                # Originally <type 'exceptions.TypeError'>: makeT1Font() argument 2
                # must be string, not None
                renderPM = None
            except OSError:
                # Probably a library problem, e.g.
                # OSError: encoder zip not available
                renderPM = None
        if circular:
            # Circular diagram
            self.gdd.draw(
                tracklines=False,
                pagesize=(15 * cm, 15 * cm),
                fragments=1,
                circle_core=0.5,
                start=0,
                end=400,
            )
            self.gdd.write(os.path.join("Graphics", name + "_c.pdf"), "pdf")

    def add_track_with_sigils(self, **kwargs):
        """Add track with sigils."""
        self.gdt_features = self.gdd.new_track(1, greytrack=False)
        self.gds_features = self.gdt_features.new_set()
        for i in range(18):
            start = int((400 * i) / 18.0)
            end = start + 17
            if i % 3 == 0:
                strand = None
                name = "Strandless"
                color = colors.orange
            elif i % 3 == 1:
                strand = +1
                name = "Forward"
                color = colors.red
            else:
                strand = -1
                name = "Reverse"
                color = colors.blue
            feature = SeqFeature(SimpleLocation(start, end, strand=strand))
            self.gds_features.add_feature(
                feature, name=name, color=color, label=True, **kwargs
            )

    def test_label_default(self):
        """Feature labels - default."""
        self.add_track_with_sigils()
        self.finish("labels_default")


class SigilsTest(unittest.TestCase):
    """Check the different feature sigils.

    These figures are intended to be used in the Tutorial...
    """

    def setUp(self):
        """Initialise diagram."""
        self.gdd = Diagram(
            "Test Diagram",
            circular=False,
            y=0.01,
            yt=0.01,
            yb=0.01,
            x=0.01,
            xl=0.01,
            xr=0.01,
        )

    def add_track_with_sigils(self, track_caption="", **kwargs):
        """Add a track of features."""
        self.gdt_features = self.gdd.new_track(
            1, greytrack=(track_caption != ""), name=track_caption, greytrack_labels=1
        )
        # We'll just use one feature set for these features,
        self.gds_features = self.gdt_features.new_set()
        # Add three features to show the strand options,
        feature = SeqFeature(SimpleLocation(25, 125, strand=+1))
        self.gds_features.add_feature(feature, name="Forward", **kwargs)
        feature = SeqFeature(SimpleLocation(150, 250, strand=None))
        self.gds_features.add_feature(feature, name="strandless", **kwargs)
        feature = SeqFeature(SimpleLocation(275, 375, strand=-1))
        self.gds_features.add_feature(feature, name="Reverse", **kwargs)

    def finish(self, name, circular=True):
        """Draw it..."""
        tracks = len(self.gdd.tracks)
        # Work around the page orientation code being too clever
        # and flipping the h & w round:
        if tracks <= 3:
            orient = "landscape"
        else:
            orient = "portrait"
        self.gdd.draw(
            format="linear",
            orientation=orient,
            tracklines=False,
            pagesize=(15 * cm, 5 * cm * tracks),
            fragments=1,
            start=0,
            end=400,
        )
        self.gdd.write(os.path.join("Graphics", name + ".pdf"), "pdf")
        global renderPM
        if renderPM:
            # For the tutorial this might be useful:
            try:
                self.gdd.write(os.path.join("Graphics", name + ".png"), "png")
            except renderPM.RenderPMError:
                # Probably a font problem
                renderPM = None
        if circular:
            # Circular diagram
            self.gdd.draw(
                tracklines=False,
                pagesize=(15 * cm, 15 * cm),
                fragments=1,
                circle_core=0.5,
                start=0,
                end=400,
            )
            self.gdd.write(os.path.join("Graphics", name + "_c.pdf"), "pdf")

    def test_all_sigils(self):
        """All sigils."""
        for glyph in ["BOX", "OCTO", "JAGGY", "ARROW", "BIGARROW"]:
            self.add_track_with_sigils(track_caption=f'  sigil="{glyph}"', sigil=glyph)
        self.finish("GD_sigils")

    def test_labels(self):
        """Feature labels."""
        self.add_track_with_sigils(label=True)
        self.add_track_with_sigils(
            label=True,
            color="green",
            # label_position left as default!
            label_size=25,
            label_angle=0,
        )
        self.add_track_with_sigils(
            label=True,
            color="purple",
            label_position="end",
            label_size=4,
            label_angle=90,
        )
        self.add_track_with_sigils(
            label=True,
            color="blue",
            label_position="middle",
            label_size=6,
            label_angle=-90,
        )
        self.add_track_with_sigils(
            label=True,
            color="cyan",
            label_position="start",
            label_size=6,
            label_angle=-90,
        )
        self.assertEqual(len(self.gdd.tracks), 5)
        self.finish("GD_sigil_labels", circular=True)

    def test_arrow_shafts(self):
        """Feature arrow sigils, varying shafts."""
        self.add_track_with_sigils(sigil="ARROW")
        self.add_track_with_sigils(sigil="ARROW", color="brown", arrowshaft_height=1.0)
        self.add_track_with_sigils(sigil="ARROW", color="teal", arrowshaft_height=0.2)
        self.add_track_with_sigils(
            sigil="ARROW", color="darkgreen", arrowshaft_height=0.1
        )
        self.assertEqual(len(self.gdd.tracks), 4)
        self.finish("GD_sigil_arrow_shafts")

    def test_big_arrow_shafts(self):
        """Feature big-arrow sigils, varying shafts."""
        self.add_track_with_sigils(sigil="BIGARROW")
        self.add_track_with_sigils(
            sigil="BIGARROW", color="orange", arrowshaft_height=1.0
        )
        self.add_track_with_sigils(
            sigil="BIGARROW", color="teal", arrowshaft_height=0.2
        )
        self.add_track_with_sigils(
            sigil="BIGARROW", color="green", arrowshaft_height=0.1
        )
        self.assertEqual(len(self.gdd.tracks), 4)
        self.finish("GD_sigil_bigarrow_shafts")

    def test_arrow_heads(self):
        """Feature arrow sigils, varying heads."""
        self.add_track_with_sigils(sigil="ARROW")
        self.add_track_with_sigils(sigil="ARROW", color="blue", arrowhead_length=0.25)
        self.add_track_with_sigils(sigil="ARROW", color="orange", arrowhead_length=1)
        self.add_track_with_sigils(
            sigil="ARROW", color="red", arrowhead_length=10000
        )  # Triangles
        self.assertEqual(len(self.gdd.tracks), 4)
        self.finish("GD_sigil_arrows")

    def short_sigils(self, glyph):
        """Draw sigils on top of grey box backgrounds."""
        # The blue boxes are only relevant for the BIGARROW
        # Add a track of features, bigger height to emphasise any sigil errors
        self.gdt_features = self.gdd.new_track(1, greytrack=True, height=3)
        # We'll just use one feature set for these features,
        self.gds_features = self.gdt_features.new_set()
        # For the ARROW and BIGARROW sigils:
        # - Green arrows just have small heads (meaning if there is a mitre
        #   it will escape the bounding box).
        # - Red arrows should be small triangles (so short no shaft shown)

        # Forward strand:
        feature = SeqFeature(SimpleLocation(15, 30, strand=-1))
        self.gds_features.add_feature(feature, color="blue")
        feature = SeqFeature(SimpleLocation(15, 30, strand=+1))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Forward", sigil=glyph, arrowhead_length=0.05
        )

        feature = SeqFeature(SimpleLocation(55, 60, strand=-1))
        self.gds_features.add_feature(feature, color="blue")
        feature = SeqFeature(SimpleLocation(55, 60, strand=+1))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Forward", sigil=glyph, arrowhead_length=1000, color="red"
        )

        feature = SeqFeature(SimpleLocation(75, 125, strand=-1))
        self.gds_features.add_feature(feature, color="blue")
        feature = SeqFeature(SimpleLocation(75, 125, strand=+1))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Forward", sigil=glyph, arrowhead_length=0.05
        )

        # Strandless:
        feature = SeqFeature(SimpleLocation(140, 155, strand=None))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Strandless", sigil=glyph, arrowhead_length=0.05
        )

        feature = SeqFeature(SimpleLocation(180, 185, strand=None))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Strandless", sigil=glyph, arrowhead_length=1000, color="red"
        )

        feature = SeqFeature(SimpleLocation(200, 250, strand=None))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Strandless", sigil=glyph, arrowhead_length=0.05
        )

        # Reverse strand:
        feature = SeqFeature(SimpleLocation(265, 280, strand=+1))
        self.gds_features.add_feature(feature, color="blue")
        feature = SeqFeature(SimpleLocation(265, 280, strand=-1))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Reverse", sigil=glyph, arrowhead_length=0.05
        )

        feature = SeqFeature(SimpleLocation(305, 310, strand=+1))
        self.gds_features.add_feature(feature, color="blue")
        feature = SeqFeature(SimpleLocation(305, 310, strand=-1))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Reverse", sigil=glyph, arrowhead_length=1000, color="red"
        )

        feature = SeqFeature(SimpleLocation(325, 375, strand=+1))
        self.gds_features.add_feature(feature, color="blue")
        feature = SeqFeature(SimpleLocation(325, 375, strand=-1))
        self.gds_features.add_feature(feature, color="grey")
        self.gds_features.add_feature(
            feature, name="Reverse", sigil=glyph, arrowhead_length=0.05
        )

        self.finish(f"GD_sigil_short_{glyph}")

    def test_short_arrow(self):
        """Feature arrow sigil heads within bounding box."""
        self.short_sigils("ARROW")

    def test_short_bigarrow(self):
        """Feature big-arrow sigil heads within bounding box."""
        self.short_sigils("BIGARROW")

    def test_short_jaggy(self):
        """Feature arrow sigil heads within bounding box."""
        self.short_sigils("JAGGY")

    def test_short_octo(self):
        """Feature big-arrow sigil heads within bounding box."""
        self.short_sigils("OCTO")

    def long_sigils(self, glyph):
        """Check feature sigils within bounding box."""
        # Add a track of features, bigger height to emphasise any sigil errors
        self.gdt_features = self.gdd.new_track(1, greytrack=True, height=3)
        # We'll just use one feature set for these features if strand specific
        self.gds_features = self.gdt_features.new_set()
        if glyph in ["BIGARROW"]:
            # These straddle the axis, so don't want to draw them on top of each other
            feature = SeqFeature(SimpleLocation(25, 375, strand=None))
            self.gds_features.add_feature(feature, color="lightblue")
            feature = SeqFeature(SimpleLocation(25, 375, strand=+1))
        else:
            feature = SeqFeature(SimpleLocation(25, 375, strand=+1))
            self.gds_features.add_feature(feature, color="lightblue")
        self.gds_features.add_feature(
            feature, name="Forward", sigil=glyph, color="blue", arrowhead_length=2.0
        )

        if glyph in ["BIGARROW"]:
            # These straddle the axis, so don't want to draw them on top of each other
            self.gdt_features = self.gdd.new_track(1, greytrack=True, height=3)
            self.gds_features = self.gdt_features.new_set()
            feature = SeqFeature(SimpleLocation(25, 375, strand=None))
            self.gds_features.add_feature(feature, color="pink")
            feature = SeqFeature(SimpleLocation(25, 375, strand=-1))
        else:
            feature = SeqFeature(SimpleLocation(25, 375, strand=-1))
            self.gds_features.add_feature(feature, color="pink")
        self.gds_features.add_feature(
            feature, name="Reverse", sigil=glyph, color="red", arrowhead_length=2.0
        )
        # Add another track of features, bigger height to emphasise any sigil errors
        self.gdt_features = self.gdd.new_track(1, greytrack=True, height=3)
        # We'll just use one feature set for these features,
        self.gds_features = self.gdt_features.new_set()
        feature = SeqFeature(SimpleLocation(25, 375, strand=None))
        self.gds_features.add_feature(feature, color="lightgreen")
        self.gds_features.add_feature(
            feature, name="Standless", sigil=glyph, color="green", arrowhead_length=2.0
        )
        self.finish(f"GD_sigil_long_{glyph}")

    def test_long_arrow_heads(self):
        """Feature ARROW sigil heads within bounding box."""
        self.long_sigils("ARROW")

    def test_long_bigarrow_heads(self):
        """Feature BIGARROW sigil heads within bounding box."""
        self.long_sigils("BIGARROW")

    def test_long_octo_heads(self):
        """Feature OCTO sigil heads within bounding box."""
        self.long_sigils("OCTO")

    def test_long_jaggy(self):
        """Feature JAGGY sigil heads within bounding box."""
        self.long_sigils("JAGGY")


class DiagramTest(unittest.TestCase):
    """Creating feature sets, graph sets, tracks etc individually for the diagram."""

    def setUp(self):
        """Test setup, just loads a GenBank file as a SeqRecord."""
        with open(os.path.join("GenBank", "NC_005816.gb")) as handle:
            self.record = SeqIO.read(handle, "genbank")

        self.gdd = Diagram("Test Diagram")
        # Add a track of features,
        self.gdd.new_track(
            1, greytrack=True, name="CDS Features", greytrack_labels=0, height=0.5
        )

    def tearDown(self):
        """Release the drawing objects."""
        del self.gdd

    def test_str(self):
        """Test diagram's info as string."""
        expected = (
            "\n<<class 'Bio.Graphics.GenomeDiagram._Diagram.Diagram'>: Test Diagram>"
            "\n1 tracks"
            "\nTrack 1: "
            "\n<<class 'Bio.Graphics.GenomeDiagram._Track.Track'>: CDS Features>"
            "\n0 sets"
            "\n"
        )
        self.assertEqual(expected, str(self.gdd))

    def test_add_track(self):
        """Add track."""
        track = Track(name="Annotated Features")
        self.gdd.add_track(track, 2)
        self.assertEqual(2, len(self.gdd.get_tracks()))

    def test_add_track_to_occupied_level(self):
        """Add track to occupied level."""
        new_track = self.gdd.get_tracks()[0]
        self.gdd.add_track(new_track, 1)
        self.assertEqual(2, len(self.gdd.get_tracks()))

    def test_add_track_error(self):
        """Test adding unspecified track."""
        self.assertRaises(ValueError, self.gdd.add_track, None, 1)

    def test_del_tracks(self):
        """Delete track."""
        self.gdd.del_track(1)
        self.assertEqual(0, len(self.gdd.get_tracks()))

    def test_get_tracks(self):
        """Get track."""
        self.assertEqual(1, len(self.gdd.get_tracks()))

    def test_move_track(self):
        """Move a track."""
        self.gdd.move_track(1, 2)
        expected = (
            "\n<<class 'Bio.Graphics.GenomeDiagram._Diagram.Diagram'>: Test Diagram>"
            "\n1 tracks"
            "\nTrack 2: "
            "\n<<class 'Bio.Graphics.GenomeDiagram._Track.Track'>: CDS Features>"
            "\n0 sets"
            "\n"
        )
        self.assertEqual(expected, str(self.gdd))

    def test_renumber(self):
        """Test renumbering tracks."""
        self.gdd.renumber_tracks(0)
        expected = (
            "\n<<class 'Bio.Graphics.GenomeDiagram._Diagram.Diagram'>: Test Diagram>"
            "\n1 tracks"
            "\nTrack 0: "
            "\n<<class 'Bio.Graphics.GenomeDiagram._Track.Track'>: CDS Features>"
            "\n0 sets"
            "\n"
        )
        self.assertEqual(expected, str(self.gdd))

    def test_write_arguments(self):
        """Check how the write methods respond to output format arguments."""
        gdd = Diagram("Test Diagram")
        gdd.drawing = None  # Hack - need the ReportLab drawing object to be created.
        filename = os.path.join("Graphics", "error.txt")
        # We (now) allow valid formats in any case.
        for output in ["XXX", "xxx", None, 123, 5.9]:
            with self.assertRaises(ValueError):
                gdd.write(filename, output)
            with self.assertRaises(ValueError):
                gdd.write_to_string(output)

    def test_partial_diagram(self):
        """Construct and draw SVG and PDF for just part of a SeqRecord."""
        genbank_entry = self.record
        start = 6500
        end = 8750

        gdd = Diagram(
            "Test Diagram",
            # For the circular diagram we don't want a closed circle:
            circular=False,
        )
        # Add a track of features,
        gdt_features = gdd.new_track(
            1,
            greytrack=True,
            name="CDS Features",
            scale_largetick_interval=1000,
            scale_smalltick_interval=100,
            scale_format="SInt",
            greytrack_labels=False,
            height=0.5,
        )
        # We'll just use one feature set for these features,
        gds_features = gdt_features.new_set()
        for feature in genbank_entry.features:
            if feature.type != "CDS":
                # We're going to ignore these.
                continue
            # These may miss fuzzy locations where the integer sorting is a simplification
            if feature.location.end < start:
                # Out of frame (too far left)
                continue
            if feature.location.start > end:
                # Out of frame (too far right)
                continue

            # This URL should work in SVG output from recent versions
            # of ReportLab.  You need ReportLab 2.4 or later
            try:
                url = (
                    "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id="
                    + str(feature.qualifiers["protein_id"][0])
                )
            except KeyError:
                url = None

            # Note that I am using strings for color names, instead
            # of passing in color objects.  This should also work!
            if len(gds_features) % 2 == 0:
                color = "white"  # for testing the automatic black border!
            else:
                color = "red"
            # Checking it can cope with the old UK spelling colour.
            # Also show the labels perpendicular to the track.
            gds_features.add_feature(
                feature,
                colour=color,
                url=url,
                sigil="ARROW",
                label_position=None,
                label_size=8,
                label_angle=90,
                label=True,
            )

        # And draw it...
        gdd.draw(
            format="linear",
            orientation="landscape",
            tracklines=False,
            pagesize=(10 * cm, 6 * cm),
            fragments=1,
            start=start,
            end=end,
        )
        output_filename = os.path.join("Graphics", "GD_region_linear.pdf")
        gdd.write(output_filename, "PDF")

        # Also check the write_to_string (bytes string) method matches,
        with open(output_filename, "rb") as handle:
            self.assertEqual(handle.read(), gdd.write_to_string("PDF"))

        output_filename = os.path.join("Graphics", "GD_region_linear.svg")
        gdd.write(output_filename, "SVG")

        # Circular with a particular start/end is a bit odd, but by setting
        # circular=False (above) a sweep of 90% is used (a wedge is left out)
        gdd.draw(
            format="circular",
            tracklines=False,
            pagesize=(10 * cm, 10 * cm),
            start=start,
            end=end,
        )
        output_filename = os.path.join("Graphics", "GD_region_circular.pdf")
        gdd.write(output_filename, "PDF")
        output_filename = os.path.join("Graphics", "GD_region_circular.svg")
        gdd.write(output_filename, "SVG")

    def test_diagram_via_methods_pdf(self):
        """Construct and draw PDF using method approach."""
        genbank_entry = self.record
        gdd = Diagram("Test Diagram")

        # Add a track of features,
        gdt_features = gdd.new_track(
            1, greytrack=True, name="CDS Features", greytrack_labels=0, height=0.5
        )
        # We'll just use one feature set for the genes and misc_features,
        gds_features = gdt_features.new_set()
        for feature in genbank_entry.features:
            if feature.type == "gene":
                if len(gds_features) % 2 == 0:
                    color = "blue"
                else:
                    color = "lightblue"
                gds_features.add_feature(
                    feature,
                    color=color,
                    # label_position="middle",
                    # label_position="end",
                    label_position="start",
                    label_size=11,
                    # label_angle=90,
                    sigil="ARROW",
                    label=True,
                )

        # I want to include some strandless features, so for an example
        # will use EcoRI recognition sites etc.
        for site, name, color in [
            ("GAATTC", "EcoRI", "green"),
            ("CCCGGG", "SmaI", "orange"),
            ("AAGCTT", "HindIII", "red"),
            ("GGATCC", "BamHI", "purple"),
        ]:
            index = 0
            while True:
                index = genbank_entry.seq.find(site, start=index)
                if index == -1:
                    break
                feature = SeqFeature(SimpleLocation(index, index + 6, strand=None))

                # This URL should work in SVG output from recent versions
                # of ReportLab.  You need ReportLab 2.4 or later
                try:
                    url = (
                        "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&id="
                        + str(feature.qualifiers["protein_id"][0])
                    )
                except KeyError:
                    url = None

                gds_features.add_feature(
                    feature,
                    color=color,
                    url=url,
                    # label_position="middle",
                    label_size=10,
                    label_color=color,
                    # label_angle=90,
                    name=name,
                    label=True,
                )
                index += len(site)
            del index

        # Now add a graph track...
        gdt_at_gc = gdd.new_track(
            2, greytrack=True, name="AT and GC content", greytrack_labels=True
        )
        gds_at_gc = gdt_at_gc.new_set(type="graph")

        step = len(genbank_entry) // 200
        gds_at_gc.new_graph(
            apply_to_window(genbank_entry.seq, step, calc_gc_content, step),
            "GC content",
            style="line",
            color=colors.lightgreen,
            altcolor=colors.darkseagreen,
        )
        gds_at_gc.new_graph(
            apply_to_window(genbank_entry.seq, step, calc_at_content, step),
            "AT content",
            style="line",
            color=colors.orange,
            altcolor=colors.red,
        )

        # Finally draw it in both formats,
        gdd.draw(
            format="linear",
            orientation="landscape",
            tracklines=0,
            pagesize="A4",
            fragments=3,
        )
        output_filename = os.path.join("Graphics", "GD_by_meth_linear.pdf")
        gdd.write(output_filename, "PDF")

        gdd.draw(
            format="circular",
            tracklines=False,
            circle_core=0.8,
            pagesize=(20 * cm, 20 * cm),
            circular=True,
        )
        output_filename = os.path.join("Graphics", "GD_by_meth_circular.pdf")
        gdd.write(output_filename, "PDF")

    def test_diagram_via_object_pdf(self):
        """Construct and draw PDF using object approach."""
        genbank_entry = self.record
        gdd = Diagram("Test Diagram")

        gdt1 = Track(
            "CDS features",
            greytrack=True,
            scale_largetick_interval=1e4,
            scale_smalltick_interval=1e3,
            greytrack_labels=10,
            greytrack_font_color="red",
            scale_format="SInt",
        )
        gdt2 = Track("gene features", greytrack=1, scale_largetick_interval=1e4)

        # First add some feature sets:
        gdfsA = FeatureSet(name="CDS backgrounds")
        gdfsB = FeatureSet(name="gene background")

        gdfs1 = FeatureSet(name="CDS features")
        gdfs2 = FeatureSet(name="gene features")
        gdfs3 = FeatureSet(name="misc_features")
        gdfs4 = FeatureSet(name="repeat regions")

        prev_gene = None
        cds_count = 0
        for feature in genbank_entry.features:
            if feature.type == "CDS":
                cds_count += 1
                if prev_gene:
                    # Assuming it goes with this CDS!
                    if cds_count % 2 == 0:
                        dark, light = colors.peru, colors.tan
                    else:
                        dark, light = colors.burlywood, colors.bisque
                    # Background for CDS,
                    a = gdfsA.add_feature(
                        SeqFeature(
                            SimpleLocation(
                                feature.location.start, feature.location.end, strand=0
                            )
                        ),
                        color=dark,
                    )
                    # Background for gene,
                    b = gdfsB.add_feature(
                        SeqFeature(
                            SimpleLocation(
                                prev_gene.location.start,
                                prev_gene.location.end,
                                strand=0,
                            )
                        ),
                        color=dark,
                    )
                    # Cross link,
                    gdd.cross_track_links.append(CrossLink(a, b, light, dark))
                    prev_gene = None
            if feature.type == "gene":
                prev_gene = feature

        # Some cross links on the same linear diagram fragment,
        f, c = fill_and_border(colors.red)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(2220, 2230)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(2200, 2210)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c))

        f, c = fill_and_border(colors.blue)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(2150, 2200)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(2220, 2290)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c, flip=True))

        f, c = fill_and_border(colors.green)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(2250, 2560)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(2300, 2860)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c))

        # Some cross links where both parts are saddling the linear diagram fragment boundary,
        f, c = fill_and_border(colors.red)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(3155, 3250)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(3130, 3300)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c))
        # Nestled within that (drawn on top),
        f, c = fill_and_border(colors.blue)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(3160, 3275)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(3180, 3225)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, f, c, flip=True))

        # Some cross links where two features are on either side of the linear diagram fragment boundary,
        f, c = fill_and_border(colors.green)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(6450, 6550)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(6265, 6365)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c))
        f, c = fill_and_border(colors.gold)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(6265, 6365)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(6450, 6550)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c))
        f, c = fill_and_border(colors.red)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(6275, 6375)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(6430, 6530)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c, flip=True))
        f, c = fill_and_border(colors.blue)
        a = gdfsA.add_feature(SeqFeature(SimpleLocation(6430, 6530)), color=f, border=c)
        b = gdfsB.add_feature(SeqFeature(SimpleLocation(6275, 6375)), color=f, border=c)
        gdd.cross_track_links.append(CrossLink(a, b, color=f, border=c, flip=True))

        cds_count = 0
        for feature in genbank_entry.features:
            if feature.type == "CDS":
                cds_count += 1
                if cds_count % 2 == 0:
                    gdfs1.add_feature(feature, color=colors.pink, sigil="ARROW")
                else:
                    gdfs1.add_feature(feature, color=colors.red, sigil="ARROW")

            if feature.type == "gene":
                # Note we set the colour of ALL the genes later on as a test,
                gdfs2.add_feature(feature, sigil="ARROW")

            if feature.type == "misc_feature":
                gdfs3.add_feature(feature, color=colors.orange)

            if feature.type == "repeat_region":
                gdfs4.add_feature(feature, color=colors.purple)

        # gdd.cross_track_links = gdd.cross_track_links[:1]

        gdfs1.set_all_features("label", 1)
        gdfs2.set_all_features("label", 1)
        gdfs3.set_all_features("label", 1)
        gdfs4.set_all_features("label", 1)

        gdfs3.set_all_features("hide", 0)
        gdfs4.set_all_features("hide", 0)

        # gdfs1.set_all_features('color', colors.red)
        gdfs2.set_all_features("color", colors.blue)

        gdt1.add_set(gdfsA)  # Before CDS so under them!
        gdt1.add_set(gdfs1)

        gdt2.add_set(gdfsB)  # Before genes so under them!
        gdt2.add_set(gdfs2)

        gdt3 = Track(
            "misc features and repeats", greytrack=1, scale_largetick_interval=1e4
        )
        gdt3.add_set(gdfs3)
        gdt3.add_set(gdfs4)

        # Now add some graph sets:

        # Use a fairly large step so we can easily tell the difference
        # between the bar and line graphs.
        step = len(genbank_entry) // 200
        gdgs1 = GraphSet("GC skew")

        graphdata1 = apply_to_window(genbank_entry.seq, step, calc_gc_skew, step)
        gdgs1.new_graph(
            graphdata1,
            "GC Skew",
            style="bar",
            color=colors.violet,
            altcolor=colors.purple,
        )

        gdt4 = Track(
            "GC Skew (bar)", height=1.94, greytrack=1, scale_largetick_interval=1e4
        )
        gdt4.add_set(gdgs1)

        gdgs2 = GraphSet("GC and AT Content")
        gdgs2.new_graph(
            apply_to_window(genbank_entry.seq, step, calc_gc_content, step),
            "GC content",
            style="line",
            color=colors.lightgreen,
            altcolor=colors.darkseagreen,
        )

        gdgs2.new_graph(
            apply_to_window(genbank_entry.seq, step, calc_at_content, step),
            "AT content",
            style="line",
            color=colors.orange,
            altcolor=colors.red,
        )

        gdt5 = Track(
            "GC Content(green line), AT Content(red line)",
            height=1.94,
            greytrack=1,
            scale_largetick_interval=1e4,
        )
        gdt5.add_set(gdgs2)

        gdgs3 = GraphSet("Di-nucleotide count")
        step = len(genbank_entry) // 400  # smaller step
        gdgs3.new_graph(
            apply_to_window(genbank_entry.seq, step, calc_dinucleotide_counts, step),
            "Di-nucleotide count",
            style="heat",
            color=colors.red,
            altcolor=colors.orange,
        )
        gdt6 = Track("Di-nucleotide count", height=0.5, greytrack=False, scale=False)
        gdt6.add_set(gdgs3)

        # Add the tracks (from both features and graphs)
        # Leave some white space in the middle/bottom
        gdd.add_track(gdt4, 3)  # GC skew
        gdd.add_track(gdt5, 4)  # GC and AT content
        gdd.add_track(gdt1, 5)  # CDS features
        gdd.add_track(gdt2, 6)  # Gene features
        gdd.add_track(gdt3, 7)  # Misc features and repeat feature
        gdd.add_track(gdt6, 8)  # Feature depth

        # Finally draw it in both formats, and full view and partial
        gdd.draw(
            format="circular", orientation="landscape", tracklines=0, pagesize="A0"
        )
        output_filename = os.path.join("Graphics", "GD_by_obj_circular.pdf")
        gdd.write(output_filename, "PDF")

        gdd.circular = False
        gdd.draw(
            format="circular",
            orientation="landscape",
            tracklines=0,
            pagesize="A0",
            start=3000,
            end=6300,
        )
        output_filename = os.path.join("Graphics", "GD_by_obj_frag_circular.pdf")
        gdd.write(output_filename, "PDF")

        gdd.draw(
            format="linear",
            orientation="landscape",
            tracklines=0,
            pagesize="A0",
            fragments=3,
        )
        output_filename = os.path.join("Graphics", "GD_by_obj_linear.pdf")
        gdd.write(output_filename, "PDF")

        gdd.set_all_tracks("greytrack_labels", 2)
        gdd.draw(
            format="linear",
            orientation="landscape",
            tracklines=0,
            pagesize=(30 * cm, 10 * cm),
            fragments=1,
            start=3000,
            end=6300,
        )
        output_filename = os.path.join("Graphics", "GD_by_obj_frag_linear.pdf")
        gdd.write(output_filename, "PDF")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
