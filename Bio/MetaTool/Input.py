# Copyright 2001 by Tarjei Mikkelsen. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module provides code to convert Bio.Pathway.System objects into
a text file that can be used as input for the MetaTool program.

For more information on MetaTool, please refer to:

http://www2.bioinf.mdc-berlin.de/metabolic/metatool/

"""

from Bio.Pathway import Reaction
from Bio.Pathway import System


# NOTE: Change enzyme name creation to allow for a function to be passed
#       in that can be applied to r.data to extract a name

def system_to_metatool(system, metext = [], metint = [], generate_names = 1):
    """Converts a Bio.Pathway.System object to a MetaTool input string.

    Note that to be a valid input string, the enzyme names of the reactions
    in the system must conform to the MetaTool requirements.

    Enzyme names are automatically genrated from the catalys attribute of
    each reaction using the following scheme:

    enzyme_name = "_".join([str(x[0]) for x in r.catalysts])

    If an enzyme name has already been used, a positive integer will be
    appended to this name to meet the MetaTool input requirements. If this
    behaviour is undesired, set the optional parameter 'generate_names' to
    false. All enzyme names will the be 'E_x', where x is an unique integer.

    The optional parameters metext and metint can be used to specify the
    external and internal metabolites according to the following rules:

    1. If metext is set, the species in it will be considered external.
       All other species will be considered internal.

    2. Otherwise, if metint is set, the species in it will be considered
       internal. All other species will be considered external.

    3. Otherwise, all species will be considered external.

    If specified, metext and metint must not contains species that are not
    contained in the input system.
    """
    if not isinstance(system, System):
        raise TypeError("Input is not a System object")
    # Build the ENZREV and ENZIRREV strings:
    enzrev    = []
    enzirrev  = []
    enz_names = {}
    enz_map   = {}
    for r in system.reactions():
        # build an unique enzyme name
        enz_name = ""
        if r.catalysts and generate_names:
            enz_name += "_".join([str(x[0]) for x in r.catalysts])
        else:
            enz_name += "E"
        if enz_name in enz_names:
            enz_names[enz_name] += 1
            enz_name += "_" + str(enz_names[enz_name])
        else:
            enz_names[enz_name] = 0
        # keep (name, reaction) pair for later
        enz_map[enz_name] = r
        # add to the corresponding list
        if r.reversible:
            enzrev.append(enz_name)
        else:
            enzirrev.append(enz_name)
    # create the actual strings:
    enzrev_str   = "-ENZREV\n" + " ".join(enzrev) + "\n"
    enzirrev_str = "-ENZIRREV\n" + " ".join(enzirrev) + "\n"       
    # Build the METINT and METEXT strings:
    metabolites = system.species()
    metint_str = "-METINT\n"
    metext_str = "-METEXT\n"
    if metext:
        for m in metext:
            if m in metabolites:
                metabolites.remove(m)
            else:
                raise ValueError("metext contains an unknown metabolite")
        for m in metint:
            if not m in metabolites:
                raise ValueError("metint contains an unknown metabolite")
        metext_str += " ".join([str(m) for m in metext]) + "\n"
        metint_str += " ".join([str(m) for m in metabolites]) + "\n"
    elif metint:
        for m in metint:
            if m in metabolites:
                metabolites.remove(m)
            else:
                raise ValueError("metint contains an unknown metabolite")
        for m in metext:
            if not m in metabolites:
                raise ValueError("metext contains an unknown metabolite")
        metint_str += " ".join([str(m) for m in metint]) + "\n"
        metext_str += " ".join([str(m) for m in metabolites]) + "\n"        
    else:
        metext_str += " ".join([str(m) for m in metabolites]) + "\n"
        metint_str += "\n"
    # Build the CAT string
    cat_str = "-CAT\n"
    for e in enz_map.keys():
        r = enz_map[e]
        cat_str += e + " : "
        # change the reaction string rep to the MetaTool format
        reaction_str = str(r)
        reaction_str = reaction_str.replace("-->","=")
        reaction_str = reaction_str.replace("<=>","=")
        cat_str += reaction_str + " .\n"
    # Return the complete MetaTool input string:
    return enzrev_str   + "\n" + \
           enzirrev_str + "\n" + \
           metint_str   + "\n" + \
           metext_str   + "\n" + \
           cat_str      + "\n"




