# Copyright 2001 by Tarjei Mikkelsen.
# Copyright 2012 by Kevin Wu.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to work with data from the KEGG database.

References:

Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes.
Nucleic Acids Res. 28, 29-34 (2000).

URL: http://www.genome.ad.jp/kegg/
"""
KEGG_ITEM_LENGTH = 12
KEGG_LINE_LENGTH = 80
KEGG_DATA_LENGTH = KEGG_LINE_LENGTH - KEGG_ITEM_LENGTH

# wrap rule = [indent, connect, (splitstr, connect, splitafter, keep), ...]
_default_wrap = lambda indent: [indent, "", (" ", "", 1, 0)]


def _wrap_kegg(line, max_width=KEGG_DATA_LENGTH, wrap_rule=_default_wrap):
    """Wraps the input line  for KEGG output.

    Arguments:

    o info - String holding the information we want wrapped
    for KEGG output.
    o max_width - Maximum width of a line.
    o wrap_rule - A wrap rule (see above) for deciding how to split
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
    """Write a indented KEGG record item.

    Arguments:

    o item - The name of the item to be written.
    o info - The (wrapped) information to write.
    o indent - Width of item field.
    """
    s = ""
    for line in info:
        partial_lines = line.splitlines()
        for l in partial_lines:
            s = s + item.ljust(indent) + l + "\n"
            if item is not "":  # ensure item is only written on first line
                item = ""
    return s


def _q(op, arg1, arg2=None, arg3=None):

  import urllib2
  URL = "http://rest.kegg.jp/%s"

  if arg2 and arg3:
    args = "%s/%s/%s/%s" % (op, arg1, arg2, arg3)
  elif arg2:
    args = "%s/%s/%s" % (op, arg1, arg2)
  else:
    args = "%s/%s" % (op, arg1)

  req = urllib2.Request(URL % (args))
  resp = urllib2.urlopen(req)

  return resp


# http://www.kegg.jp/kegg/rest/keggapi.html
def info(arg1):
    # http://rest.kegg.jp/info/<database>
    #
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome |<org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme | genomes | genes | ligand | kegg
    # <org> = KEGG organism code or T number
    if arg1:
        resp = _q("info", arg1)
    else:
        raise Exception("Invalid info operation.")

    return resp


def list_(arg1, arg2=None):
    #  http://rest.kegg.jp/list/<database>/<org>
    #
    #  <database> = pathway | module
    #  <org> = KEGG organism code
    if isinstance(arg1, str) and (arg1 == "pathway" or arg1 == "module") and arg2:
        resp = _q("list", arg1, arg2)

    # http://rest.kegg.jp/list/<database>
    #
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome | <org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme | organism
    # <org> = KEGG organism code or T number
    #
    #
    # http://rest.kegg.jp/list/<dbentries>
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome | <org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme
    # <org> = KEGG organism code or T number
    elif arg1 and not arg2:
        if isinstance(arg1, list) and len(arg1) <= 100:
            arg1 = ("+").join(arg1)
        elif isinstance(arg1, list) and len(arg1) > 100:
            raise Exception("Maximuim number of queries is 100")
        resp = _q("list", arg1)

    else:
        raise Exception("Invalid list operation.")

    return resp


def find(arg1, arg2=None, arg3=None):
    # http://rest.kegg.jp/find/<database>/<query>/<option>
    #
    # <database> = compound | drug
    # <option> = formula | exact_mass | mol_weight
    if arg1 in ["compound", "drug"] and arg2 and \
          arg3 in ["formula", "exact_mass", "mol_weight"]:
        resp = _q("find", arg1, arg2, arg3)

    # http://rest.kegg.jp/find/<database>/<query>
    #
    # <database> = pathway | module | disease | drug | environ | ko |
    #              genome | <org> | compound | glycan | reaction | rpair |
    #              rclass | enzyme | genes | ligand
    # <org> = KEGG organism code or T number
    elif arg1 and arg2 and not arg3:

        if isinstance(arg2, list):
            arg2 = "+".join(arg2)
        resp = _q("find", arg1, arg2)

    else:
        raise Exception("Invalid find operation.")

    return resp


def get(arg1, arg2=None):
    if isinstance(arg1, list) and len(arg1) <= 10:
        arg1 = "+".join(arg1)
    elif isinstance(arg1, list) and len(arg1) > 10:
        raise Exception("Maximuim number of queries is 10")

    # http://rest.kegg.jp/get/<dbentries>[/<option>]
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome | <org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme
    # <org> = KEGG organism code or T number
    #
    # <option> = aaseq | ntseq | mol | kcf | image
    if arg1 and arg2 in ["aaseq" , "ntseq", "mol", "kcf", "image"]:
        resp = _q("get", arg1, arg2)
    elif arg1 and not arg2:
        resp = _q("get", arg1)
    else:
        raise Exception("Invalid get operation.")

    return resp


def conv(arg1, arg2=None):
    # http://rest.kegg.jp/conv/<target_db>/<source_db>
    #
    # (<target_db> <source_db>) = (<kegg_db> <outside_db>) |
    #                             (<outside_db> <kegg_db>)
    #
    # For gene identifiers:
    # <kegg_db> = <org>
    # <org> = KEGG organism code or T number
    # <outside_db> = ncbi-gi | ncbi-geneid | uniprot
    #
    # For chemical substance identifiers:
    # <kegg_db> = drug | compound | glycan
    # <outside_db> = pubchem | chebi
    if ((arg1 and arg2 in ["ncbi-gi", "ncbi-geneid", "uniprot"]) or
        (arg1 in ["ncbi-gi", "ncbi-geneid", "uniprot"] and arg2) or
        (arg1 in ["drug", "compound", "glycan"] and
         arg2 in ["pubchem", "glycan"]) or
        (arg1 in ["pubchem", "glycan"] and
         arg2 in ["drug", "compound", "glycan"])) \
       and isinstance(arg2, str):
        resp = _q("conv", arg1, arg2)

    # http://rest.kegg.jp/conv/<target_db>/<dbentries>
    #
    # For gene identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = <org> | ncbi-gi | ncbi-geneid | uniprot
    # <org> = KEGG organism code or T number
    #
    # For chemical substance identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = drug | compound | glycan | pubchem | chebi
    elif arg1 and arg2:
        if isinstance(arg2, list):
            arg2 = "+".join(arg2)
        resp = _q("conv", arg1, arg2)

    else:
        raise Exception("Invalid conv operation.")

    return resp


def link(arg1, arg2):
    # http://rest.kegg.jp/link/<target_db>/<source_db>
    #
    # <target_db> = <database>
    # <source_db> = <database>
    #
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome | <org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme
    #
    # http://rest.kegg.jp/link/<target_db>/<dbentries>
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome | <org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme
    if arg1 and arg2:
        if isinstance(arg2, list):
            arg2 = "+".join(arg2)
        resp = _q("link", arg1, arg2)

    else:
        raise Exception("Invalid link operation.")

    return resp
