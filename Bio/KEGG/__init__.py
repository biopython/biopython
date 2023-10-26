# Copyright 2001 by Tarjei Mikkelsen.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Code to work with data from the KEGG database.

References:
Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes.
Nucleic Acids Res. 28, 29-34 (2000).

URL: http://www.genome.ad.jp/kegg/

"""


KEGG_ITEM_LENGTH = 12
KEGG_LINE_LENGTH = 80
KEGG_DATA_LENGTH = KEGG_LINE_LENGTH - KEGG_ITEM_LENGTH


def _default_wrap(indent):
    """Return default wrap rule for _wrap_kegg (PRIVATE).

    A wrap rule is a list with the following elements:
    [indent, connect, (splitstr, connect, splitafter, keep), ...]
    """
    return [indent, "", (" ", "", 1, 0)]


def _struct_wrap(indent):
    """Return wrap rule for KEGG STRUCTURE (PRIVATE)."""
    return [indent, "", ("  ", "", 1, 1)]


def _wrap_kegg(line, max_width=KEGG_DATA_LENGTH, wrap_rule=_default_wrap):
    """Wrap the input line for KEGG output (PRIVATE).

    Arguments:
     - info - String holding the information we want wrapped
       for KEGG output.
     - max_width - Maximum width of a line.
     - wrap_rule - A wrap rule (see above) for deciding how to split
       strings that must be wrapped.

    """
    s = ""
    wrapped_line = ""
    indent = " " * wrap_rule[0]
    connect = wrap_rule[1]
    rules = wrap_rule[2:]
    while True:
        if len(line) <= max_width:
            wrapped_line = wrapped_line + line
            s = s + wrapped_line
            break
        else:
            did_split = 0
            for rule in rules:
                to = max_width
                if not rule[2]:
                    to = to + len(rule[0])
                split_idx = line.rfind(rule[0], 0, to)
                if split_idx > -1:
                    if rule[2] and rule[3]:
                        split_idx = split_idx + len(rule[0])
                    wrapped_line = wrapped_line + line[0:split_idx] + "\n"
                    if not rule[3]:
                        split_idx = split_idx + len(rule[0])
                    line = indent + rule[1] + line[split_idx:]
                    did_split = 1
                    break
            if not did_split:
                wrapped_line = wrapped_line + line[0:max_width] + "\n"
                line = indent + connect + line[max_width:]
    return s


def _write_kegg(item, info, indent=KEGG_ITEM_LENGTH):
    """Write a indented KEGG record item (PRIVATE).

    Arguments:
     - item - The name of the item to be written.
     - info - The (wrapped) information to write.
     - indent - Width of item field.

    """
    s = ""
    for line in info:
        partial_lines = line.splitlines()
        for partial in partial_lines:
            s += item.ljust(indent) + partial + "\n"
            if item:  # ensure item is only written on first line
                item = ""
    return s
