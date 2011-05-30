# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement
import os
import os.path
from _paml import Paml, PamlError
import _parse_baseml

class BasemlError(EnvironmentError):
    """BASEML has failed. Run with verbose = True to view BASEML's error
message"""

class Baseml(Paml):
    """This class implements an interface to BASEML, part of the PAML package."""

    def __init__(self, alignment = None, tree = None, working_dir = None,
                out_file = None):
        """Initialize the Baseml instance. 
        
        The user may optionally pass in strings specifying the locations
        of the input alignment and tree files, the working directory and
        the final output file. 
        """
        Paml.__init__(self, alignment, working_dir, out_file)
        if tree is not None:
            if not os.path.exists(tree):
                raise IOError, "The specified tree file does not exist."
        self.tree = tree
        self.ctl_file = "baseml.ctl"
        self._options = {"noisy": 2,
                        "verbose": 0,
                        "runmode": 0,
                        "model": 7,
                        "model_options": None,
                        "Mgene": 0,
                        "ndata": None,
                        "clock": 0,
                        "fix_kappa": 0,
                        "kappa": 5,
                        "fix_alpha": 0,
                        "alpha": 0.5,
                        "Malpha": 0,
                        "ncatG": 5,
                        "fix_rho": 1,
                        "rho": 0,
                        "nparK": 0,
                        "nhomo": 0,
                        "getSE": 0,
                        "RateAncestor": 0,
                        "Small_Diff": 7e-6,
                        "cleandata": 1,
                        "icode": None,
                        "fix_blength": None,
                        "method": 0}

    def write_ctl_file(self):
        """Dynamically build a BASEML control file from the options.

        The control file is written to the location specified by the 
        ctl_file property of the baseml class.
        """
        # Make sure all paths are relative to the working directory
        try:
            self._set_rel_paths()
        except (AttributeError, ValueError) as error:
            raise error
        with open(self.ctl_file, 'w') as ctl_handle:
            ctl_handle.write("seqfile = {0}\n".format(self._rel_alignment))
            ctl_handle.write("outfile = {0}\n".format(self._rel_out_file))
            ctl_handle.write("treefile = {0}\n".format(self._rel_tree))
            for option in self._options.items():
                if option[1] == None:
                    # If an option has a value of None, there's no need
                    # to write it in the control file; it's normally just
                    # commented out.
                    continue
                if option[0] == "model_options":
                    continue
                # If "model" is 9 or 10, it may be followed in the cotnrol
                # file by further options such as
                # [1 (TC CT AG GA)]
                # or
                # [5 (AC CA) (AG GA) (AT TA) (CG GC) (CT TC)]
                # which are to be stored in "model_options" as a string.
                if option[0] == "model" and option[1] in [9, 10]:
                    if self._options["model_options"] is not None:
                        ctl_handle.write("model = {0}  {1}".format(option[1],
                            self._options["model_options"]))
                        continue
                ctl_handle.write("{0} = {1}\n".format(option[0],
                    option[1]))

    def read_ctl_file(self, ctl_file):
        """Parse a control file and load the options into the Baseml instance.
        """
        temp_options = {}
        with open(ctl_file) as ctl_handle:
            for line in ctl_handle:
                line = line.strip()
                uncommented = line.partition("*")[0]
                if uncommented != "":
                    if "=" not in uncommented:
                        raise AttributeError, \
                            "Malformed line in control file:\n{0}".format(line)
                    (option, value) = uncommented.split("=")
                    option = option.strip()
                    value = value.strip()
                    if option == "seqfile":
                        self.alignment = value
                    elif option == "treefile":
                        self.tree = value
                    elif option == "outfile":
                        self.out_file = value
                    elif option not in self._options:
                        raise KeyError, "Invalid option: {0}".format(option)
                    elif option == "model":
                        if len(value) <= 2 and value.isdigit():
                            temp_options["model"] = int(value)
                            temp_options["model_options"] = None
                        else:
                            model_num = value.partition(" ")[0]
                            model_opt = value.partition(" ")[2].strip()
                            temp_options["model"] = int(model_num)
                            temp_options["model_options"] = model_opt
                    else:
                        if "." in value or "e-" in value:
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
        for option in self._options.keys():
            if option in temp_options.keys():
                self._options[option] = temp_options[option]
            else:
                self._options[option] = None

    def _set_rel_paths(self):
        """Convert all file/directory locations to paths relative to the current working directory.

        BASEML requires that all paths specified in the control file be
        relative to the directory from which it is called rather than 
        absolute paths.
        """
        Paml._set_rel_paths(self)
        if self.tree is not None:
            self._rel_tree = os.path.relpath(self.tree, 
                self.working_dir)

    def run(self, ctl_file = None, verbose = False, command = "baseml",
                parse = True):
        """Run baseml using the current configuration and then parse the results. 

        Return a process signal so the user can determine if
        the execution was successful (return code 0 is successful, -N
        indicates a failure). The arguments may be passed as either 
        absolute or relative paths, despite the fact that BASEML 
        requires relative paths.
        """
        if self.tree is None:
            raise ValueError, "Tree file not specified."
        if not os.path.exists(self.tree):
            raise IOError, "The specified tree file does not exist."
        try:
            Paml.run(self, ctl_file, verbose, command)
        except PamlError as strerror:
            raise PamlError, strerror
        if parse:
            try:
                results = read(self._rel_out_file)
            except KeyError as (errorno, strerror):
                raise KeyError, strerror
            except ValueError as (errorno, strerror):
                raise ValueError, strerror
        else:
            results = None
        return results

def read(results_file):
    results = {}
    """Parse a BASEML results file."""
    if not os.path.exists(results_file):
        raise IOError, "Results file does not exist."
    with open(results_file) as results_handle:
        lines = results_handle.readlines()
    num_params = _parse_baseml.parse_basics(lines, results)
    _parse_baseml.parse_parameters(lines, results, num_params)
    if results.get("version") is None:
        raise ValueError, "Invalid results file"
    return results
