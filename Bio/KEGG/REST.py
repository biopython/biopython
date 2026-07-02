# Copyright 2014 by Kevin Wu.
# Revisions copyright 2014 by Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Provides code to access the REST-style KEGG online API.

This module aims to make the KEGG online REST-style API easier to use. See:
https://www.kegg.jp/kegg/rest/keggapi.html

The KEGG REST-style API provides simple access to a range of KEGG databases.
This works using simple URLs (which this module will construct for you),
with any errors indicated via HTTP error levels.

The functionality is somewhat similar to Biopython's Bio.TogoWS and Bio.Entrez
modules.

Currently KEGG does not provide any usage guidelines (unlike the NCBI whose
requirements are reasonably clear). To avoid risking overloading the service,
Biopython will only allow three calls per second.

References:
Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes.
Nucleic Acids Res. 28, 29-34 (2000).

"""

import io
import time
from urllib.request import urlopen

from Bio._utils import function_with_previous


@function_with_previous
def _q(op, arg1, arg2=None, arg3=None):
    delay = 0.333333333  # one third of a second
    current = time.time()
    wait = _q.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _q.previous = current + wait
    else:
        _q.previous = current

    URL = "https://rest.kegg.jp/%s"
    if arg2 and arg3:
        args = f"{op}/{arg1}/{arg2}/{arg3}"
    elif arg2:
        args = f"{op}/{arg1}/{arg2}"
    else:
        args = f"{op}/{arg1}"
    resp = urlopen(URL % (args))

    if "image" == arg2:
        return resp

    handle = io.TextIOWrapper(resp, encoding="UTF-8")
    handle.url = resp.url
    return handle


_q.previous = 0


# https://www.kegg.jp/kegg/rest/keggapi.html
def kegg_info(database):
    """KEGG info - Displays the current statistics of a given database.

    db - database or organism (string)

    The argument db can be a KEGG database name (e.g. 'pathway' or its
    official abbreviation, 'path'), or a KEGG organism code or T number
    (e.g. 'hsa' or 'T01001' for human).

    A valid list of organism codes and their T numbers can be obtained
    via kegg_info('organism') or https://rest.kegg.jp/list/organism

    """
    # TODO - return a string (rather than the handle?)
    # TODO - cache and validate the organism code / T numbers?
    # TODO - can we parse the somewhat formatted output?
    #
    # https://rest.kegg.jp/info/<database>
    #
    # <database> = kegg | pathway | brite | module | ko | genes | <org> | ag |
    #              vg | vp | genome | vtax | vgenome | compound | glycan |
    #              reaction | rclass | rmodule | enzyme | network | ntmap |
    #              variant | disease | drug | dgroup
    # <org> = KEGG organism code or T number
    return _q("info", database)


def kegg_list(database, org=None, option=None):
    """KEGG list - Entry list for database, or specified database entries.

    db - database or organism (string)
    org - optional organism (string), see below.
    option - optional brite or genome (string), see below.

    For the pathway databases the optional organism can be used to
    restrict the results by setting `org`.
    For the brite and genome databases also can be used to restrict
    the results by setting `option`.

    """
    # TODO - split into two functions (dbentries seems separate)?
    #
    #  https://rest.kegg.jp/list/pathway/<org>
    #
    #  <org> = KEGG organism code
    if database == "pathway" and org and not option:
        return _q("list", database, org)

    # https://rest.kegg.jp/list/brite/<option>
    #
    # <option> = br | jp | ko | <org>
    #
    # https://rest.kegg.jp/list/genome/<option>
    #
    # <option> = <group_name> | <rank_id>
    # <group_name> = KEGG organism group name
    # <rank_id> = Taxonomy ID for phylum, class, order, family, genus and species
    if database in ("brite", "genome") and not org and option:
        return _q("list", database, option)

    # https://rest.kegg.jp/list/<database>
    #
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome |compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup
    # <org> = KEGG organism code or T number
    if isinstance(database, str) and not org and not option:
        return _q("list", database)

    # https://rest.kegg.jp/list/<dbentries>
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome |compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup
    # <org> = KEGG organism code or T number
    if isinstance(database, list) and not org and not option:
        if len(database) > 10:
            raise ValueError("Maximum number of databases is 10 for kegg list query")
        database = ("+").join(database)
        return _q("list", database)

    raise ValueError("Invalid database arg for kegg list request.")


def kegg_find(database, query, option=None):
    """KEGG find - Data search.

    Finds entries with matching query keywords or other query data in
    a given database.

    db - database or organism (string)
    query - search terms (string)
    option - search option (string), see below.

    For the compound and drug database, set option to the string 'formula',
    'exact_mass', 'mol_weight' or 'nop' to search on that field only. The
    chemical formula search is a partial match irrespective of the order
    of atoms given. The exact mass (or molecular weight) is checked by
    rounding off to the same decimal place as the query data. A range of
    values may also be specified with the minus(-) sign.

    """
    # TODO - return list of tuples?
    #
    # https://rest.kegg.jp/find/<database>/<query>/<option>
    #
    # <database> = compound | drug
    # <option> = formula | exact_mass | mol_weight | nop
    if database in ["compound", "drug"] and option in [
        "formula",
        "exact_mass",
        "mol_weight",
        "nop",
    ]:
        resp = _q("find", database, query, option)
    elif option:
        raise ValueError("Invalid option arg for kegg find request.")

    # https://rest.kegg.jp/find/<database>/<query>
    #
    #
    # <database> = pathway | brite | module | ko | genes | <org> | ag | vg |
    #              vp | genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | enzyme | network | ntmap | variant | disease |
    #              drug | dgroup
    # <org> = KEGG organism code or T number
    else:
        if isinstance(query, list):
            query = "+".join(query)
        resp = _q("find", database, query)

    return resp


def kegg_get(dbentries, option=None):
    """KEGG get - Data retrieval.

    dbentries - Identifiers (single string, or list of strings), see below.
    option - One of "aaseq", "ntseq", "mol", "kcf", "image", "conf",
    "kgml", "json" (string)

    The input is limited up to 10 entries.
    The input is limited to one pathway entry with the image or kgml option.
    The input is limited to one compound/glycan/drug entry with the image option.

    Returns a handle.
    """
    if isinstance(dbentries, list) and len(dbentries) <= 10:
        dbentries = "+".join(dbentries)
    elif isinstance(dbentries, list) and len(dbentries) > 10:
        raise ValueError("Maximum number of dbentries is 10 for kegg get query")

    # https://rest.kegg.jp/get/<dbentries>[/<option>]
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup | disease_ja | drug_ja | dgroup_ja |
    #              compound_ja
    # <org> = KEGG organism code or T number
    #
    # <option> = aaseq | ntseq | mol | kcf | image | conf | kgml | json
    if option in ["aaseq", "ntseq", "mol", "kcf", "image", "conf", "kgml", "json"]:
        resp = _q("get", dbentries, option)
    elif option:
        raise ValueError("Invalid option arg for kegg get request.")
    else:
        resp = _q("get", dbentries)

    return resp


def kegg_conv(target_db, source_db, option=None):
    """KEGG conv - convert KEGG identifiers to/from outside identifiers.

    Arguments:
     - target_db - Target database
     - source_db_or_dbentries - source database or database entries
     - option - depricated

    """
    # https://rest.kegg.jp/conv/<target_db>/<source_db>
    #
    # (<target_db> <source_db>) = (<kegg_db> <outside_db>) |
    #                             (<outside_db> <kegg_db>)
    #
    # For gene identifiers:
    # <kegg_db> = <org>
    # <org> = KEGG organism code or T number
    # <outside_db> = ncbi-proteinid | ncbi-geneid | uniprot
    #
    # For chemical substance identifiers:
    # <kegg_db> = drug | compound | glycan
    # <outside_db> = pubchem | chebi
    #
    # https://rest.kegg.jp/conv/<target_db>/<dbentries>
    #
    # For gene identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = <org> | ncbi-proteinid | ncbi-geneid | uniprot
    # <org> = KEGG organism code or T number
    #
    # For chemical substance identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = drug | compound | glycan | pubchem | chebi
    _ = option

    if isinstance(source_db, list):
        source_db = "+".join(source_db)

    if (
        target_db in ["ncbi-proteinid", "ncbi-geneid", "uniprot"]
        or source_db in ["ncbi-proteinid", "ncbi-geneid", "uniprot"]
        or (
            target_db in ["drug", "compound", "glycan"]
            and source_db in ["pubchem", "glycan"]
        )
        or (
            target_db in ["pubchem", "glycan"]
            and source_db in ["drug", "compound", "glycan"]
        )
    ):
        return _q("conv", target_db, source_db)
    else:
        raise ValueError("Bad argument target_db or source_db for kegg conv request.")


def kegg_link(target_db, source_db, option=None):
    """KEGG link - find related entries by using database cross-references.

    target_db - Target database
    source_db_or_dbentries - source database
    option - Can be "species", "genus", "family", "order", "class", "phylum".
    If getting data with RDF, can be "turtle" or "n-triple" (string).
    """
    if isinstance(source_db, list):
        source_db = "+".join(source_db)

    # https://rest.kegg.jp/link/<target_db>/<source_db>[/<option>]
    #
    # <target_db> = <database>
    # <source_db> = <database>
    #
    # <database> = drug | atc | jtc
    #
    # <option> = turtle | n-triple
    # https://rest.kegg.jp/link/<target_db>/<dbentries>[/<option>]
    #
    # <dbentries> = KEGG database entries of the following <database>
    # <database> = drug | atc | jtc
    #
    # <option> = turtle | n-triple
    if option in ["turtle", "n-triple"] and target_db in ["drug", "atc", "jtc"]:
        return _q("link", target_db, source_db, option)

    # https://rest.kegg.jp/link/<target_db>/<source_db>[/<option>]
    #
    # <target_db> = <database>
    # <source_db> = <database>
    #
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup | <outside_db>
    # <outside_db> = pubmed | taxonomy | atc | jtc | ndc | yk
    #
    # <option> = species | genus | family | order | class | phylum
    # https://rest.kegg.jp/link/<target_db>/<dbentries>[/<option>]
    # <database> = pathway | brite | module | ko | <org> | ag | vg | vp |
    #              genome | vtax | vgenome | compound | glycan | reaction |
    #              rclass | rmodule | enzyme | network | ntmap | variant |
    #              disease | drug | dgroup | <outside_db>
    # <outside_db> = pubmed | taxonomy | atc | jtc | ndc | yk
    #
    # <option> = species | genus | family | order | class | phylum
    if option in ["species", "genus", "family", "order", "class", "phylum"]:
        return _q("link", target_db, source_db, option)

    if not option:
        return _q("link", target_db, source_db)

    raise ValueError("Invalid option arg for kegg conv request.")
