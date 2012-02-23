import os.path
import sys


VERSIONS = ["4_1", "4_3", "4_4", "4_4c", "4_5"]
VERBOSE = True
    
def codeml():
    from Bio.Phylo.PAML import codeml
    tests = [("aa_model0", "aa_alignment.phylip", "species.tree"),
            ("aa_pairwise", "aa_alignment.phylip", "species.tree"),
            ("all_NSsites", "alignment.phylip", "species.tree"),
            ("branchsiteA", "alignment.phylip", "species.tree"),
            ("clademodelC", "alignment.phylip", "species.tree"),
            ("freeratio", "alignment.phylip", "species.tree"),
            ("ngene2_mgene012", "lysinYangSwanson2002.nuc", "lysin.trees"),
            ("ngene2_mgene34", "lysinYangSwanson2002.nuc", "lysin.trees"),
            ("pairwise", "alignment.phylip", "species.tree"),
            ("SE", "alignment.phylip", "species.tree")]

    for test in tests:
        print test[0]
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
        for version in VERSIONS:
            print "\t{0}".format(version.replace('_', '.'))
            out_file = '.'.join(['-'.join([test[0], version]), "out"])
            cml.out_file = os.path.join("Results", "codeml", test[0], out_file)
            bin = ''.join(["codeml", version])
            cml.run(command=bin, verbose=VERBOSE)
         

def baseml():
    from Bio.Phylo.PAML import baseml
    
    tests = [("model", range(0, 9)), ("nhomo", [1, 3, 4]),
            ("nparK", range(1, 5)), ("alpha1rho1", None), ("SE", None)]
    alignment = os.path.join("Alignments", "alignment.phylip")
    tree = os.path.join("Trees", "species.tree")
    for test in tests:
        print test[0]
        bml = baseml.Baseml()
        for version in VERSIONS:
            print "\t{0}".format(version.replace('_', '.'))
            if test[1] is not None:
                for n in test[1]:
                    if (version in ["4_3", "4_4", "4_4c", "4_5"] and 
                            test[0] == "nparK" and n in [3, 4]):
                        continue
                    print "\t\tn = {0}".format(n)
                    ctl_file = os.path.join("Control_files", "baseml", 
                        "{0}{1}.ctl".format(test[0], n))                   
                    bml.read_ctl_file(ctl_file)
                    bml.alignment = alignment
                    bml.tree = tree                
                    out_file = "{0}{1}-{2}.out".format(test[0], n, version)  
                    bml.out_file = os.path.join("Results", "baseml", test[0],
                        out_file)
                    bin = "baseml{0}".format(version)
                    bml.run(command=bin, verbose=VERBOSE)               
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
                bml.run(command=bin, verbose=VERBOSE)

def yn00():
    from Bio.Phylo.PAML import yn00
    
    tests = ["yn00"]
    alignment = os.path.join("Alignments", "alignment.phylip")
    for test in tests:
        print test[0]
        yn = yn00.Yn00()
        for version in VERSIONS:
            print "\t{0}".format(version.replace('_', '.'))
            ctl_file = os.path.join("Control_files", "yn00",
                "{0}.ctl".format(test))
            yn.read_ctl_file(ctl_file)
            yn.alignment = alignment
            out_file = "{0}-{1}.out".format(test, version)
            yn.out_file = os.path.join("Results", "yn00", out_file)
            bin = "yn00{0}".format(version)
            yn.run(command=bin, verbose=VERBOSE)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("specify a paml program")
    if sys.argv[1] == "codeml":
        codeml()
    elif sys.argv[1] == "baseml":
        baseml()
    elif sys.argv[1] == "yn00":
        yn00()
    else:
        sys.exit("specify a paml program")
