# Copyright (C) 2021 by Robyn Wright (R-Wright-1)
# Copyright (C) 2026 by Martin Mokrejš
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for multi-orientation tree drawing in Bio.Phylo._utils.draw().

Tests the orient_tree, horizontal_direction, vertical_direction,
circular_span, draw_labels, and align_labels parameters added in
PR #3693.
"""

import unittest
from io import StringIO

from Bio import MissingExternalDependencyError
from Bio import Phylo

try:
    import matplotlib
except ImportError:
    raise MissingExternalDependencyError(
        "Install matplotlib if you want to use Bio.Phylo._utils."
    ) from None

matplotlib.use("ps")

try:
    from matplotlib import pyplot
except ImportError:
    raise MissingExternalDependencyError(
        "Install matplotlib if you want to use Bio.Phylo._utils."
    ) from None


# Use a simple Newick tree for consistent testing
SIMPLE_NEWICK = "(((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6):0.7,E:0.8);"
LARGER_NEWICK = (
    "((((A:0.1,B:0.2):0.3,C:0.4):0.5,"
    "(D:0.6,(E:0.7,F:0.8):0.9):1.0):1.1,"
    "(G:1.2,H:1.3):1.4);"
)


def _make_tree(newick=SIMPLE_NEWICK):
    """Create a tree from a Newick string."""
    return Phylo.read(StringIO(newick), "newick")


class TestDrawOrientTreeParameter(unittest.TestCase):
    """Test the orient_tree parameter of Phylo.draw()."""

    def setUp(self):
        """Set up non-interactive matplotlib."""
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        """Close all figures."""
        pyplot.close("all")

    # --- Default behaviour (backward compatibility) ---

    def test_default_draw_unchanged(self):
        """Default draw() should behave exactly like the original."""
        tree = _make_tree()
        # Should not raise, should produce a figure
        Phylo.draw(tree, do_show=False)
        fig = pyplot.gcf()
        ax = fig.axes[0]
        # Default is vertical, right-facing
        self.assertEqual(ax.get_xlabel(), "branch length")
        self.assertEqual(ax.get_ylabel(), "taxa")

    # --- orient_tree validation ---

    def test_orient_tree_invalid_raises(self):
        """Invalid orient_tree value should raise ValueError."""
        tree = _make_tree()
        with self.assertRaises(ValueError) as ctx:
            Phylo.draw(tree, orient_tree="diagonal", do_show=False)
        self.assertIn("orient_tree", str(ctx.exception))

    def test_orient_tree_vertical(self):
        """orient_tree='vertical' should produce a standard tree."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="vertical", do_show=False)
        ax = pyplot.gcf().axes[0]
        self.assertEqual(ax.get_xlabel(), "branch length")
        self.assertEqual(ax.get_ylabel(), "taxa")

    def test_orient_tree_horizontal(self):
        """orient_tree='horizontal' should swap axes labels."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="horizontal", do_show=False)
        ax = pyplot.gcf().axes[0]
        self.assertEqual(ax.get_xlabel(), "taxa")
        self.assertEqual(ax.get_ylabel(), "branch length")

    def test_orient_tree_circular(self):
        """orient_tree='circular' should use polar projection."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="circular", do_show=False)
        ax = pyplot.gcf().axes[0]
        self.assertEqual(ax.name, "polar")


class TestVerticalDirection(unittest.TestCase):
    """Test the vertical_direction parameter."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_vertical_right(self):
        """vertical_direction='right' should have x-axis ascending."""
        tree = _make_tree()
        Phylo.draw(
            tree, orient_tree="vertical", vertical_direction="right", do_show=False
        )
        ax = pyplot.gcf().axes[0]
        xlim = ax.get_xlim()
        # Right-facing: left < right
        self.assertLess(xlim[0], xlim[1])

    def test_vertical_left(self):
        """vertical_direction='left' should have x-axis descending (flipped)."""
        tree = _make_tree()
        Phylo.draw(
            tree, orient_tree="vertical", vertical_direction="left", do_show=False
        )
        ax = pyplot.gcf().axes[0]
        xlim = ax.get_xlim()
        # Left-facing: left > right (flipped)
        self.assertGreater(xlim[0], xlim[1])

    def test_vertical_direction_invalid_raises(self):
        """Invalid vertical_direction should raise ValueError."""
        tree = _make_tree()
        with self.assertRaises(ValueError) as ctx:
            Phylo.draw(
                tree,
                orient_tree="vertical",
                vertical_direction="up",
                do_show=False,
            )
        self.assertIn("vertical_direction", str(ctx.exception))

    def test_vertical_direction_ignored_for_horizontal(self):
        """vertical_direction should be silently ignored for horizontal trees."""
        tree = _make_tree()
        # Should not raise even with a non-standard vertical_direction
        # when orient_tree is horizontal
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            vertical_direction="left",
            do_show=False,
        )


class TestHorizontalDirection(unittest.TestCase):
    """Test the horizontal_direction parameter."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_horizontal_down(self):
        """horizontal_direction='down' should have y-axis descending."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            horizontal_direction="down",
            do_show=False,
        )
        ax = pyplot.gcf().axes[0]
        ylim = ax.get_ylim()
        # Down-facing: top > bottom (descending y)
        self.assertGreater(ylim[0], ylim[1])

    def test_horizontal_up(self):
        """horizontal_direction='up' should have y-axis ascending."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            horizontal_direction="up",
            do_show=False,
        )
        ax = pyplot.gcf().axes[0]
        ylim = ax.get_ylim()
        # Up-facing: bottom < top (ascending y)
        self.assertLess(ylim[0], ylim[1])

    def test_horizontal_direction_invalid_raises(self):
        """Invalid horizontal_direction should raise ValueError."""
        tree = _make_tree()
        with self.assertRaises(ValueError) as ctx:
            Phylo.draw(
                tree,
                orient_tree="horizontal",
                horizontal_direction="left",
                do_show=False,
            )
        self.assertIn("horizontal_direction", str(ctx.exception))


class TestCircularTrees(unittest.TestCase):
    """Test circular tree drawing."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_circular_default_span(self):
        """Circular tree with default span (355°) should work."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="circular", do_show=False)

    def test_circular_custom_span(self):
        """Circular tree with custom span should work."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="circular", circular_span=270, do_show=False)

    def test_circular_full_circle(self):
        """Circular tree with span=360 should work (full circle)."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="circular", circular_span=360, do_show=False)

    def test_circular_half_circle(self):
        """Circular tree with span=180 should work (half circle)."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="circular", circular_span=180, do_show=False)

    def test_circular_needs_polar_axes(self):
        """Pre-created non-polar axes should raise ValueError for circular."""
        tree = _make_tree()
        fig, ax = pyplot.subplots()
        with self.assertRaises(ValueError) as ctx:
            Phylo.draw(tree, orient_tree="circular", axes=ax, do_show=False)
        self.assertIn("polar", str(ctx.exception))

    def test_circular_with_polar_axes(self):
        """Pre-created polar axes should work for circular trees."""
        tree = _make_tree()
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection="polar")
        Phylo.draw(tree, orient_tree="circular", axes=ax, do_show=False)

    def test_circular_with_larger_tree(self):
        """Circular tree should work with a larger tree."""
        tree = _make_tree(LARGER_NEWICK)
        Phylo.draw(tree, orient_tree="circular", do_show=False)


class TestDrawLabels(unittest.TestCase):
    """Test the draw_labels parameter."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_draw_labels_true_default(self):
        """draw_labels=True (default) should add text objects."""
        tree = _make_tree()
        Phylo.draw(tree, draw_labels=True, do_show=False)
        ax = pyplot.gcf().axes[0]
        texts = ax.texts
        # Should have at least the 5 leaf labels (A, B, C, D, E)
        self.assertGreaterEqual(len(texts), 5)

    def test_draw_labels_false(self):
        """draw_labels=False should suppress all text labels."""
        tree = _make_tree()
        Phylo.draw(tree, draw_labels=False, do_show=False)
        ax = pyplot.gcf().axes[0]
        texts = ax.texts
        # No labels should be drawn
        self.assertEqual(len(texts), 0)

    def test_draw_labels_false_horizontal(self):
        """draw_labels=False for horizontal tree should suppress labels."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            draw_labels=False,
            do_show=False,
        )
        ax = pyplot.gcf().axes[0]
        self.assertEqual(len(ax.texts), 0)

    def test_draw_labels_false_circular(self):
        """draw_labels=False for circular tree should suppress labels."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="circular",
            draw_labels=False,
            do_show=False,
        )
        ax = pyplot.gcf().axes[0]
        self.assertEqual(len(ax.texts), 0)


class TestAlignLabels(unittest.TestCase):
    """Test the align_labels parameter."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_align_labels_vertical(self):
        """align_labels=True on vertical tree should not raise."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="vertical", align_labels=True, do_show=False)

    def test_align_labels_horizontal(self):
        """align_labels=True on horizontal tree should not raise."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="horizontal", align_labels=True, do_show=False)

    def test_align_labels_circular(self):
        """align_labels=True on circular tree should not raise."""
        tree = _make_tree()
        Phylo.draw(tree, orient_tree="circular", align_labels=True, do_show=False)

    def test_align_labels_affects_ylim_horizontal_down(self):
        """align_labels should increase y-range for horizontal-down trees."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            horizontal_direction="down",
            align_labels=False,
            do_show=False,
        )
        ax_no_align = pyplot.gcf().axes[0]
        ylim_no_align = ax_no_align.get_ylim()
        pyplot.close("all")

        Phylo.draw(
            tree,
            orient_tree="horizontal",
            horizontal_direction="down",
            align_labels=True,
            do_show=False,
        )
        ax_align = pyplot.gcf().axes[0]
        ylim_align = ax_align.get_ylim()
        # Aligned labels need more space, so ylim range should be wider
        range_no_align = abs(ylim_no_align[0] - ylim_no_align[1])
        range_align = abs(ylim_align[0] - ylim_align[1])
        self.assertGreater(range_align, range_no_align)

    def test_align_labels_affects_ylim_circular(self):
        """align_labels should increase y-range for circular trees."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="circular",
            align_labels=False,
            do_show=False,
        )
        ax_no_align = pyplot.gcf().axes[0]
        ylim_no_align = ax_no_align.get_ylim()
        pyplot.close("all")

        Phylo.draw(
            tree,
            orient_tree="circular",
            align_labels=True,
            do_show=False,
        )
        ax_align = pyplot.gcf().axes[0]
        ylim_align = ax_align.get_ylim()
        self.assertGreater(ylim_align[1], ylim_no_align[1])


class TestCombinations(unittest.TestCase):
    """Test various parameter combinations for completeness."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_vertical_right_labels_aligned(self):
        """Vertical right with aligned labels should not raise."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="vertical",
            vertical_direction="right",
            align_labels=True,
            draw_labels=True,
            do_show=False,
        )

    def test_vertical_left_labels_aligned(self):
        """Vertical left with aligned labels should not raise."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="vertical",
            vertical_direction="left",
            align_labels=True,
            draw_labels=True,
            do_show=False,
        )

    def test_horizontal_up_labels_aligned(self):
        """Horizontal up with aligned labels should not raise."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            horizontal_direction="up",
            align_labels=True,
            draw_labels=True,
            do_show=False,
        )

    def test_horizontal_down_labels_aligned(self):
        """Horizontal down with aligned labels should not raise."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            horizontal_direction="down",
            align_labels=True,
            draw_labels=True,
            do_show=False,
        )

    def test_circular_aligned_no_labels(self):
        """Circular aligned with labels off should not raise."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="circular",
            align_labels=True,
            draw_labels=False,
            do_show=False,
        )

    def test_all_orientations_with_larger_tree(self):
        """All orientations should work with a larger tree."""
        tree = _make_tree(LARGER_NEWICK)
        for orient in ("vertical", "horizontal", "circular"):
            with self.subTest(orient=orient):
                Phylo.draw(tree, orient_tree=orient, do_show=False)
                pyplot.close("all")

    def test_all_orientations_with_confidence(self):
        """All orientations should render confidence values when available."""
        tree = _make_tree()
        for orient in ("vertical", "horizontal"):
            with self.subTest(orient=orient):
                Phylo.draw(
                    tree,
                    orient_tree=orient,
                    show_confidence=True,
                    do_show=False,
                )
                pyplot.close("all")

    def test_label_colors_with_orientation(self):
        """label_colors should work with all orientations."""
        tree = _make_tree()
        colors = {"A": "red", "B": "blue", "C": "green"}
        for orient in ("vertical", "horizontal", "circular"):
            with self.subTest(orient=orient):
                Phylo.draw(
                    tree,
                    orient_tree=orient,
                    label_colors=colors,
                    do_show=False,
                )
                pyplot.close("all")


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and unusual trees."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_single_node_tree(self):
        """Single-node tree should draw without errors."""
        tree = _make_tree("A:1.0;")
        for orient in ("vertical", "horizontal"):
            with self.subTest(orient=orient):
                Phylo.draw(tree, orient_tree=orient, do_show=False)
                pyplot.close("all")

    def test_two_taxa_tree(self):
        """Two-taxa tree should draw in all orientations."""
        tree = _make_tree("(A:0.1,B:0.2);")
        for orient in ("vertical", "horizontal", "circular"):
            with self.subTest(orient=orient):
                Phylo.draw(tree, orient_tree=orient, do_show=False)
                pyplot.close("all")

    def test_no_branch_lengths(self):
        """Tree without branch lengths should draw in all orientations."""
        tree = _make_tree("(((A,B),(C,D)),E);")
        for orient in ("vertical", "horizontal", "circular"):
            with self.subTest(orient=orient):
                Phylo.draw(tree, orient_tree=orient, do_show=False)
                pyplot.close("all")

    def test_custom_label_func(self):
        """Custom label_func should work with new orientations."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            orient_tree="horizontal",
            label_func=lambda c: c.name.upper() if c.name else None,
            do_show=False,
        )

    def test_pre_created_axes_vertical(self):
        """Pre-created standard axes should work for vertical trees."""
        tree = _make_tree()
        fig, ax = pyplot.subplots()
        Phylo.draw(tree, orient_tree="vertical", axes=ax, do_show=False)

    def test_pre_created_axes_horizontal(self):
        """Pre-created standard axes should work for horizontal trees."""
        tree = _make_tree()
        fig, ax = pyplot.subplots()
        Phylo.draw(tree, orient_tree="horizontal", axes=ax, do_show=False)

    def test_phyloxml_tree(self):
        """PhyloXML tree should work with new orientations."""
        import os

        xml_path = os.path.join(os.path.dirname(__file__), "PhyloXML", "apaf.xml")
        tree = Phylo.read(xml_path, "phyloxml")
        for orient in ("vertical", "horizontal", "circular"):
            with self.subTest(orient=orient):
                Phylo.draw(tree, orient_tree=orient, do_show=False)
                pyplot.close("all")


class TestBackwardCompatibility(unittest.TestCase):
    """Ensure the existing draw() API is not broken."""

    def setUp(self):
        pyplot.ioff()
        pyplot.close("all")

    def tearDown(self):
        pyplot.close("all")

    def test_existing_branch_labels_dict(self):
        """Existing branch_labels dict usage should still work."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            branch_labels={tree.root: "Root"},
            do_show=False,
        )

    def test_existing_branch_labels_callable(self):
        """Existing branch_labels callable usage should still work."""
        tree = _make_tree()
        Phylo.draw(
            tree,
            branch_labels=lambda c: c.branch_length,
            do_show=False,
        )

    def test_existing_show_confidence_false(self):
        """show_confidence=False should still work."""
        tree = _make_tree()
        Phylo.draw(tree, show_confidence=False, do_show=False)

    def test_title_from_tree_name(self):
        """Tree with a name should have it as the axes title."""
        tree = _make_tree()
        tree.name = "Test Tree"
        Phylo.draw(tree, do_show=False)
        ax = pyplot.gcf().axes[0]
        self.assertEqual(ax.get_title(), "Test Tree")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
