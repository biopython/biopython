#!/usr/bin/env python
# Copyright (C) 2021 by Robyn Wright (R-Wright-1)
# Copyright (C) 2026 by Martin Mokrejš
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Example: Drawing phylogenetic trees in multiple orientations.

Demonstrates the orient_tree, vertical_direction, horizontal_direction,
circular_span, draw_labels and align_labels parameters of Phylo.draw().

Requires matplotlib, numpy and scipy.

Usage::

    python phylo_draw_orientations.py

This will create 'phylo_draw_orientations.png' in the current directory.
"""

from io import StringIO

from Bio import Phylo

try:
    import matplotlib

    matplotlib.use("Agg")  # non-interactive backend
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("This example requires matplotlib.") from None

# A small tree to illustrate all orientations
NEWICK = "(((Alpha:0.1,Beta:0.25):0.3,(Gamma:0.4,Delta:0.5):0.6):0.7,Epsilon:0.8);"
tree = Phylo.read(StringIO(NEWICK), "newick")

fig = plt.figure(figsize=(20, 24))

# --- 1. Default vertical tree (root on the left) ---
ax1 = fig.add_subplot(3, 2, 1)
Phylo.draw(
    tree, orient_tree="vertical", vertical_direction="right", axes=ax1, do_show=False
)
ax1.set_title("Vertical (default, root left)", fontweight="bold")

# --- 2. Vertical tree facing left (root on the right) ---
ax2 = fig.add_subplot(3, 2, 2)
Phylo.draw(
    tree, orient_tree="vertical", vertical_direction="left", axes=ax2, do_show=False
)
ax2.set_title("Vertical (root right)", fontweight="bold")

# --- 3. Horizontal tree facing down (root at top) ---
ax3 = fig.add_subplot(3, 2, 3)
Phylo.draw(
    tree,
    orient_tree="horizontal",
    horizontal_direction="down",
    align_labels=True,
    axes=ax3,
    do_show=False,
)
ax3.set_title("Horizontal down + aligned labels", fontweight="bold")

# --- 4. Horizontal tree facing up (root at bottom) ---
ax4 = fig.add_subplot(3, 2, 4)
Phylo.draw(
    tree,
    orient_tree="horizontal",
    horizontal_direction="up",
    axes=ax4,
    do_show=False,
)
ax4.set_title("Horizontal up", fontweight="bold")

# --- 5. Circular tree (default 355° span) with aligned labels ---
ax5 = fig.add_subplot(3, 2, 5, projection="polar")
Phylo.draw(
    tree,
    orient_tree="circular",
    align_labels=True,
    axes=ax5,
    do_show=False,
)
ax5.set_title("Circular (355°) + aligned labels", fontweight="bold", pad=20)

# --- 6. Circular tree (270° span) without labels ---
ax6 = fig.add_subplot(3, 2, 6, projection="polar")
Phylo.draw(
    tree,
    orient_tree="circular",
    circular_span=270,
    draw_labels=False,
    axes=ax6,
    do_show=False,
)
ax6.yaxis.grid(False)
ax6.set_xticks([])
ax6.set_yticklabels([])
ax6.set_title("Circular (270°) no labels", fontweight="bold", pad=20)

plt.subplots_adjust(hspace=0.35, wspace=0.3)
outfile = "phylo_draw_orientations.png"
plt.savefig(outfile, bbox_inches="tight", dpi=150)
print(f"Saved {outfile}")
