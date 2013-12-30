# Copyright (C) 2012, 2013 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import print_function

import os.path
import sys

from Bio._py3k import range


VERSIONS = ["4_1", "4_3", "4_4", "4_4c", "4_5", "4_6", "4_7"]


def codeml(vers=None, verbose=False):
    from Bio.Phylo.PAML import codeml
    if vers is not None:
        versions = [vers]
    else:
        versions = VERSIONS
    tests = [("aa_model0", "aa_alignment.phylip", "species.tree"),
            ("aa_pairwise", "aa_alignment.phylip", "species.tree"),
            ("all_NSsites", "alignment.phylip", "species.tree"),
            ("branchsiteA", "alignment.phylip", "species.tree"),
            ("clademodelC", "alignment.phylip", "species.tree"),
            ("freeratio", "alignment.phylip", "species.tree"),
            ("ngene2_mgene02", "lysinYangSwanson2002.nuc", "lysin.trees"),
            ("ngene2_mgene34", "lysinYangSwanson2002.nuc", "lysin.trees"),
            ("pairwise", "alignment.phylip", "species.tree"),
            ("SE", "alignment.phylip", "species.tree")]

    for test in tests:
        print(test[0])
        cml = codeml.Codeml()
        cml.working_dir = "temp"
        ctl_file = os.path.join("Control_files",
                                "codeml",
                                '.'.join([test[0], "ctl"]))
        alignment = os.path.join("Alignments", test[1])
        tree = os.path.join("Trees", test[2])
        cml.read_ctl_file(ctl_file)
        cml.alignment = alignment
        cml.tree = tree
        for version in versions:
            print("\t{0}".format(version.replace('_', '.')))
            if test[0] in ["ngene2_mgene02", "ngene2_mgene34"] and \
               version == "4_6":
                cml.tree = ".".join([cml.tree, "4.6"])
            out_file = '.'.join(['-'.join([test[0], version]), "out"])
            cml.out_file = os.path.join("Results", "codeml", test[0], out_file)
            bin = ''.join(["codeml", version])
            cml.run(command=bin, verbose=verbose)


def baseml(vers=None, verbose=False):
    from Bio.Phylo.PAML import baseml
    if vers is not None:
        versions = [vers]
    else:
        versions = VERSIONS
    tests = [("model", list(range(0, 9))), ("nhomo", [1, 3, 4]),
            ("nparK", list(range(1, 5))), ("alpha1rho1", None), ("SE", None)]
    alignment = os.path.join("Alignments", "alignment.phylip")
    tree = os.path.join("Trees", "species.tree")
    for test in tests:
        print(test[0])
        bml = baseml.Baseml()
        for version in versions:
            print("\t{0}".format(version.replace('_', '.')))
            if test[1] is not None:
                for n in test[1]:
                    if (version in ["4_3", "4_4", "4_4c", "4_5"] and
                            test[0] == "nparK" and n in [3, 4]):
                        continue
                    print("\t\tn = {0}".format(n))
                    ctl_file = os.path.join("Control_files", "baseml",
                        "{0}{1}.ctl".format(test[0], n))
                    bml.read_ctl_file(ctl_file)
                    bml.alignment = alignment
                    bml.tree = tree
                    out_file = "{0}{1}-{2}.out".format(test[0], n, version)
                    bml.out_file = os.path.join("Results", "baseml", test[0],
                        out_file)
                    bin = "baseml{0}".format(version)
                    bml.run(command=bin, verbose=verbose)
            else:
                if (version in ["4_3", "4_4", "4_4c", "4_5"] and
                        test[0] == "alpha1rho1"):
                    continue
                ctl_file = os.path.join("Control_files", "baseml",
                    "{0}.ctl".format(test[0]))
                bml.read_ctl_file(ctl_file)
                bml.alignment = alignment
                bml.tree = tree
                out_file = "{0}-{1}.out".format(test[0], version)
                bml.out_file = os.path.join("Results", "baseml", test[0],
                    out_file)
                bin = "baseml{0}".format(version)
                bml.run(command=bin, verbose=verbose)


def yn00(vers=None, verbose=False):
    from Bio.Phylo.PAML import yn00
    if vers is not None:
        versions = [vers]
    else:
        versions = VERSIONS
    tests = ["yn00"]
    alignment = os.path.join("Alignments", "alignment.phylip")
    for test in tests:
        print(test[0])
        yn = yn00.Yn00()
        for version in versions:
            print("\t{0}".format(version.replace('_', '.')))
            ctl_file = os.path.join("Control_files", "yn00",
                "{0}.ctl".format(test))
            yn.read_ctl_file(ctl_file)
            yn.alignment = alignment
            out_file = "{0}-{1}.out".format(test, version)
            yn.out_file = os.path.join("Results", "yn00", out_file)
            bin = "yn00{0}".format(version)
            yn.run(command=bin, verbose=verbose)


def print_usage():
    versions = ", ".join(vers.replace("_", ".") for vers in VERSIONS)
    usage = \
'''Usage: gen_results.py [-v] PROGRAM [VERSION]
Generate result files to be used in Bio.Phylo.PAML unit tests.

  -v         Use verbose output
  PROGRAM    codeml, baseml or yn00
  VERSION    %s

To use this, the PAML programs must be in your executable path and
they must be named programX_Y, where X and Y are the version numbers
(i.e. baseml4_5 or codeml4_4c). If VERSION is not specified, test
results will be generated for all versions listed above.
'''%(versions)
    sys.exit(usage)

if __name__ == "__main__":
    programs = ["codeml", "baseml", "yn00"]
    prog = None
    verbose = False
    vers = None
    if len(sys.argv) < 2:
        print_usage()
    for arg in sys.argv[1:]:
        if arg == "-v":
            verbose = True
        elif arg in programs:
            if prog is not None:
                print("Only one program at a time, please.")
                print_usage()
            prog = arg
        elif arg.replace(".", "_") in VERSIONS:
            if vers is not None:
                print("Only one version at a time, sorry.")
            vers = arg.replace(".", "_")
        else:
            print("Unrecognized argument")
            print_usage()
    if prog is None:
        print("No program specified")
        print_usage()
    if prog == "codeml":
        codeml(vers, verbose)
    elif prog == "baseml":
        baseml(vers, verbose)
    elif prog == "yn00":
        yn00(vers, verbose)
