# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import os
import os.path
from _paml import Paml, PamlError, _relpath
import _parse_codeml

#TODO - Restore use of with statement for closing handles automatically
#after dropping Python 2.4

class CodemlError(EnvironmentError):
    """CODEML has failed. Run with verbose = True to view CODEML's error
message"""

class Codeml(Paml):
    """This class implements an interface to CODEML, part of the PAML package."""

    def __init__(self, alignment = None, tree = None, working_dir = None,
                out_file = None):
        """Initialize the codeml instance. 
        
        The user may optionally pass in strings specifying the locations
        of the input alignment and tree files, the working directory and
        the final output file. Other options found in the CODEML control
        have typical settings by default to run site class models 0, 1 and
        2 on a nucleotide alignment.
        """
        Paml.__init__(self, alignment, working_dir, out_file)
        if tree is not None:
            if not os.path.exists(tree):
                raise IOError, "The specified tree file does not exist."
        self.tree = tree
        self.ctl_file = "codeml.ctl"
        self._options = {"noisy": None, 
                        "verbose": None, 
                        "runmode": None,
                        "seqtype": None, 
                        "CodonFreq": None, 
                        "ndata": None,
                        "clock": None, 
                        "aaDist": None,
                        "aaRatefile": None, 
                        "model": None,
                        "NSsites": None, 
                        "icode": None, 
                        "Mgene": None,
                        "fix_kappa": None, 
                        "kappa": None, 
                        "fix_omega": None,
                        "omega": None, 
                        "fix_alpha": None, 
                        "alpha": None,
                        "Malpha": None, 
                        "ncatG": None, 
                        "getSE": None,
                        "RateAncestor": None, 
                        "Small_Diff": None,
                        "cleandata": None, 
                        "fix_blength": None, 
                        "method": None,
                        "rho": None,
                        "fix_rho": None}
        
    def write_ctl_file(self):
        """Dynamically build a CODEML control file from the options.
        
        The control file is written to the location specified by the 
        ctl_file property of the codeml class.
        """
        # Make sure all paths are relative to the working directory
        self._set_rel_paths()
        if True: #Dummy statement to preserve indentation for diff
            ctl_handle = open(self.ctl_file, 'w')
            ctl_handle.write("seqfile = %s\n" % self._rel_alignment)
            ctl_handle.write("outfile = %s\n" % self._rel_out_file)
            ctl_handle.write("treefile = %s\n" % self._rel_tree)
            for option in self._options.items():
                if option[1] == None:
                    # If an option has a value of None, there's no need
                    # to write it in the control file; it's normally just
                    # commented out.
                    continue
                if option[0] == "NSsites":
                    # NSsites is stored in Python as a list but in the 
                    # control file it is specified as a series of numbers
                    # separated by spaces.
                    NSsites = " ".join([str(site) for site in option[1]])
                    ctl_handle.write("%s = %s\n" % (option[0], NSsites))
                else:
                    ctl_handle.write("%s = %s\n" % (option[0], option[1]))
            ctl_handle.close()
    
    def read_ctl_file(self, ctl_file):
        """Parse a control file and load the options into the Codeml instance.
        """
        temp_options = {}
        if not os.path.isfile(ctl_file):
            raise IOError("File not found: %r" % ctl_file)
        else:
            ctl_handle = open(ctl_file)
            for line in ctl_handle:
                line = line.strip()
                uncommented = line.split("*",1)[0]
                if uncommented != "":
                    if "=" not in uncommented:
                        ctl_handle.close()
                        raise AttributeError, \
                            "Malformed line in control file:\n%r" % line
                    (option, value) = uncommented.split("=")
                    option = option.strip()
                    value = value.strip()
                    if option == "seqfile":
                        self.alignment = value
                    elif option == "treefile":
                        self.tree = value
                    elif option == "outfile":
                        self.out_file = value
                    elif option == "NSsites":
                        site_classes = value.split(" ")
                        for n in range(len(site_classes)):
                            try:
                                site_classes[n] = int(site_classes[n])
                            except:
                                ctl_handle.close()
                                raise TypeError, \
                                    "Invalid site class: %s" % site_classes[n]
                        temp_options["NSsites"] = site_classes
                    elif option not in self._options:
                        ctl_handle.close()
                        raise KeyError, "Invalid option: %s" % option
                    else:
                        if "." in value:
                            try:
                                converted_value = float(value)
                            except:
                                converted_value = value
                        else:
                            try:
                                converted_value = int(value)
                            except:
                                converted_value = value
                        temp_options[option] = converted_value
            ctl_handle.close()
        for option in self._options.keys():
            if option in temp_options.keys():
                self._options[option] = temp_options[option]
            else:
                self._options[option] = None
                            
    def print_options(self):
        """Print out all of the options and their current settings."""
        for option in self._options.items():
            if option[0] == "NSsites" and option[1] is not None:
                # NSsites is stored in Python as a list but in the 
                # control file it is specified as a series of numbers
                # separated by spaces.
                NSsites = " ".join([str(site) for site in option[1]])
                print "%s = %s" % (option[0], NSsites)
            else:
                print "%s = %s" % (option[0], option[1])
        
    def _set_rel_paths(self):
        """Convert all file/directory locations to paths relative to the current working directory.
        
        CODEML requires that all paths specified in the control file be
        relative to the directory from which it is called rather than 
        absolute paths.
        """
        Paml._set_rel_paths(self)
        if self.tree is not None:
            self._rel_tree = _relpath(self.tree, self.working_dir)
        
    def run(self, ctl_file = None, verbose = False, command = "codeml",
                parse = True):
        """Run codeml using the current configuration and then parse the results. 
        
        Return a process signal so the user can determine if
        the execution was successful (return code 0 is successful, -N
        indicates a failure). The arguments may be passed as either 
        absolute or relative paths, despite the fact that CODEML 
        requires relative paths.
        """
        if self.tree is None:
            raise ValueError, "Tree file not specified."
        if not os.path.exists(self.tree):
            raise IOError, "The specified tree file does not exist."
        Paml.run(self, ctl_file, verbose, command)
        if parse:
            results = read(self.out_file)
        else:
            results = None
        return results        

def read(results_file):
    """Parse a CODEML results file."""
    results = {}
    if not os.path.exists(results_file):
        raise IOError, "Results file does not exist."
    handle = open(results_file)
    lines = handle.readlines()
    handle.close()
    (results, multi_models) = _parse_codeml.parse_basics(lines, results)
    results = _parse_codeml.parse_nssites(lines, results, multi_models)
    results = _parse_codeml.parse_pairwise(lines, results)
    results = _parse_codeml.parse_distances(lines, results)
    if len(results) == 0:
        raise ValueError, "Invalid results file"
    return results
