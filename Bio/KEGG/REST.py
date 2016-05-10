# Copyright 2014 by Kevin Wu.
# Revisions copyright 2014 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access the REST-style KEGG online API.

This module aims to make the KEGG online REST-style API easier to use. See:
http://www.kegg.jp/kegg/rest/keggapi.html

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

from Bio._py3k import urlopen as _urlopen


def _q(op, arg1, arg2=None, arg3=None):
    URL = "http://rest.kegg.jp/%s"
    if arg2 and arg3:
        args = "%s/%s/%s/%s" % (op, arg1, arg2, arg3)
    elif arg2:
        args = "%s/%s/%s" % (op, arg1, arg2)
    else:
        args = "%s/%s" % (op, arg1)
    return _urlopen(URL % (args))


# http://www.kegg.jp/kegg/rest/keggapi.html
def kegg_info(database):
    """KEGG info - Displays the current statistics of a given database.

    db - database or organism (string)

    The argument db can be a KEGG database name (e.g. 'pathway' or its
    official abbreviation, 'path'), or a KEGG organism code or T number
    (e.g. 'hsa' or 'T01001' for human).

    A valid list of organism codes and their T numbers can be obtained
    via kegg_info('organism') or http://rest.kegg.jp/list/organism

    """
    # TODO - return a string (rather than the handle?)
    # TODO - chache and validate the organism code / T numbers?
    # TODO - can we parse the somewhat formatted output?
    #
    # http://rest.kegg.jp/info/<database>
    #
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome |<org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme | genomes | genes | ligand | kegg
    # <org> = KEGG organism code or T number
    return _q("info", database)


def kegg_list(database, org=None):
    """KEGG list - Entry list for database, or specified database entries.

    db - database or organism (string)
    org - optional organism (string), see below.

    For the pathway and module databases the optional organism can be
    used to restrict the results.

    """
    # TODO - split into two functions (dbentries seems separate)?
    #
    #  http://rest.kegg.jp/list/<database>/<org>
    #
    #  <database> = pathway | module
    #  <org> = KEGG organism code
    if isinstance(database, str) and (database in ["pathway", "module"]) and org:
        resp = _q("list", database, org)
    elif isinstance(database, str) and database and org:
        raise Exception("Invalid database arg for kegg list request.")

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
    else:
        if isinstance(database, list) and len(database) <= 100:
            database = ("+").join(database)
        elif isinstance(database, list) and len(database) > 100:
            raise Exception("Maximuim number of databases is 100 for kegg list query")
        resp = _q("list", database)

    return resp


def kegg_find(database, query, option=None):
    """KEGG find - Data search.

    Finds entries with matching query keywords or other query data in
    a given database.

    db - database or organism (string)
    query - search terms (string)
    option - search option (string), see below.

    For the compound and drug database, set option to the string 'formula',
    'exact_mass' or 'mol_weight' to search on that field only. The
    chemical formula search is a partial match irrespective of the order
    of atoms given. The exact mass (or molecular weight) is checked by
    rounding off to the same decimal place as the query data. A range of
    values may also be specified with the minus(-) sign.

    """
    # TODO - return list of tuples?
    #
    # http://rest.kegg.jp/find/<database>/<query>/<option>
    #
    # <database> = compound | drug
    # <option> = formula | exact_mass | mol_weight
    if database in ["compound", "drug"] and \
          option in ["formula", "exact_mass", "mol_weight"]:
        resp = _q("find", database, query, option)
    elif option:
        raise Exception("Invalid option arg for kegg find request.")

    # http://rest.kegg.jp/find/<database>/<query>
    #
    # <database> = pathway | module | disease | drug | environ | ko |
    #              genome | <org> | compound | glycan | reaction | rpair |
    #              rclass | enzyme | genes | ligand
    # <org> = KEGG organism code or T number
    else:
        if isinstance(query, list):
            query = "+".join(query)
        resp = _q("find", database, query)

    return resp


def kegg_get(dbentries, option=None):
    """KEGG get - Data retrieval.

    dbentries - Identifiers (single string, or list of strings), see below.
    option - One of "aaseq", "ntseq", "mol", "kcf", "image", "kgml" (string)

    The input is limited up to 10 entries.
    The input is limited to one pathway entry with the image or kgml option.
    The input is limited to one compound/glycan/drug entry with the image option.

    Returns a handle.
    """
    if isinstance(dbentries, list) and len(dbentries) <= 10:
        dbentries = "+".join(dbentries)
    elif isinstance(dbentries, list) and len(dbentries) > 10:
        raise Exception("Maximum number of dbentries is 10 for kegg get query")

    # http://rest.kegg.jp/get/<dbentries>[/<option>]
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | disease | drug | environ |
    #              ko | genome | <org> | compound | glycan | reaction |
    #              rpair | rclass | enzyme
    # <org> = KEGG organism code or T number
    #
    # <option> = aaseq | ntseq | mol | kcf | image
    if option in ["aaseq", "ntseq", "mol", "kcf", "image", "kgml"]:
        resp = _q("get", dbentries, option)
    elif option:
        raise Exception("Invalid option arg for kegg get request.")
    else:
        resp = _q("get", dbentries)

    return resp


def kegg_conv(target_db, source_db, option=None):
    """KEGG conv - convert KEGG identifiers to/from outside identifiers

    target_db - Target database
    source_db_or_dbentries - source database or database entries
    option - Can be "turtle" or "n-triple" (string).
    """
    # http://rest.kegg.jp/conv/<target_db>/<source_db>[/<option>]
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
    #
    # <option> = turtle | n-triple
    #
    # http://rest.kegg.jp/conv/<target_db>/<dbentries>[/<option>]
    #
    # For gene identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = <org> | ncbi-gi | ncbi-geneid | uniprot
    # <org> = KEGG organism code or T number
    #
    # For chemical substance identifiers:
    # <dbentries> = database entries involving the following <database>
    # <database> = drug | compound | glycan | pubchem | chebi
    #
    # <option> = turtle | n-triple
    if option and option not in ["turtle", "n-triple"]:
        raise Exception("Invalid option arg for kegg conv request.")

    if isinstance(source_db, list):
        source_db = "+".join(source_db)

    if target_db in ["ncbi-gi", "ncbi-geneid", "uniprot"] or \
       source_db in ["ncbi-gi", "ncbi-geneid", "uniprot"] or \
       (target_db in ["drug", "compound", "glycan"] and
           source_db in ["pubchem", "glycan"]) or \
       (target_db in ["pubchem", "glycan"] and
           source_db in ["drug", "compound", "glycan"]):

        if option:
            resp = _q("conv", target_db, source_db, option)
        else:
            resp = _q("conv", target_db, source_db)

        return resp
    else:
        raise Exception("Bad argument target_db or source_db for kegg conv request.")


def kegg_link(target_db, source_db, option=None):
    """KEGG link - find related entries by using database cross-references.

    target_db - Target database
    source_db_or_dbentries - source database
    option - Can be "turtle" or "n-triple" (string).
    """
    # http://rest.kegg.jp/link/<target_db>/<source_db>[/<option>]
    #
    # <target_db> = <database>
    # <source_db> = <database>
    #
    # <database> = pathway | brite | module | ko | genome | <org> | compound |
    #              glycan | reaction | rpair | rclass | enzyme | disease |
    #              drug | dgroup | environ
    #
    # <option> = turtle | n-triple
    # http://rest.kegg.jp/link/<target_db>/<dbentries>[/<option>]
    #
    # <dbentries> = KEGG database entries involving the following <database>
    # <database> = pathway | brite | module | ko | genome | <org> | compound |
    #              glycan | reaction | rpair | rclass | enzyme | disease |
    #              drug | dgroup | environ | genes
    #
    # <option> = turtle | n-triple

    if option and option not in ["turtle", "n-triple"]:
        raise Exception("Invalid option arg for kegg conv request.")

    if isinstance(source_db, list):
        source_db = "+".join(source_db)

    if option:
        resp = _q("link", target_db, source_db, option)
    else:
        resp = _q("link", target_db, source_db)

    return resp
