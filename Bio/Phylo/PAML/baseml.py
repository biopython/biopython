# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement
import os
import os.path
import subprocess
import re


class BasemlError(EnvironmentError):
    """BASEML has failed. Run with verbose = True to view BASEML's error
message"""

class Baseml(object):
    """This class implements an interface to BASEML, part of the PAML package."""

    def __init__(self, alignment = None, tree = None, working_dir = None,
                out_file = None):
        """Initialize the Baseml instance. 
        
        The user may optionally pass in strings specifying the locations
        of the input alignment and tree files, the working directory and
        the final output file. 
        """
        self.ctl_file = "baseml.ctl"
        if working_dir is None:
            self.working_dir = os.getcwd()
        else:
            self.working_dir = working_dir
        if alignment is not None:
            if not os.path.exists(alignment):
                raise IOError, "The specified alignment file does not exist."
        self.alignment = alignment
        if tree is not None:
            if not os.path.exists(tree):
                raise IOError, "The specified tree file does not exist."
        self.tree = tree
        self.out_file = out_file
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
                            temp_options["model"] = int(value)
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
                            
    def print_options(self):
        """Print out all of the options and their current settings."""
        for option in self._options.items():
            print "{0} = {1}".format(option[0], option[1])
    
    def set_option(self, option, value):
        """Set the value of an option. 
        
        This function abstracts the options dict to prevent the user from 
        adding options that do not exist or mispelling options.
        """
        if not self._options.has_key(option):
            raise KeyError, "Invalid option: " + option
        else:
            self._options[option] = value
    
    def get_option(self, option):
        """Return the value of an option."""
        if not self._options.has_key(option):
            raise KeyError, "Invalid option: " + option
        else:
            return self._options.get(option)
        
    def get_all_options(self):
        """Return the values of all the options."""        
        return self._options.items()
        
    def _set_rel_paths(self):
        """Convert all file/directory locations to paths relative to the current working directory.
        
        BASEML requires that all paths specified in the control file be
        relative to the directory from which it is called rather than 
        absolute paths.
        """
        if self.working_dir is not None:
            self._rel_working_dir = os.path.relpath(self.working_dir)
        if self.alignment is not None:
            self._rel_alignment = os.path.relpath(self.alignment, 
                self.working_dir)
        if self.tree is not None:
            self._rel_tree = os.path.relpath(self.tree, 
                self.working_dir)
        if self.out_file is not None:
            self._rel_out_file = os.path.relpath(self.out_file, 
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
        if self.alignment is None:
            raise ValueError, "Alignment file not specified."
        if not os.path.exists(self.alignment):
            raise IOError, "The specified alignment file does not exist."
        if self.tree is None:
            raise ValueError, "Tree file not specified."
        if not os.path.exists(self.tree):
            raise IOError, "The specified tree file does not exist."
        if self.out_file is None:
            raise ValueError, "Output file not specified."
        if self.working_dir is None:
            raise ValueError, "Working directory not specified."
        # Store the current working directory
        cwd = os.getcwd()
        # Move to the desired working directory
        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)
        os.chdir(self.working_dir)
        # If no external control file was specified...
        if ctl_file is None:
            # Dynamically build a control file
            self.write_ctl_file()
            if verbose:
                result_code = subprocess.call([command, self.ctl_file])
            else:
                result_code = subprocess.call([command, self.ctl_file],
                    stdout=open(os.devnull, 'w'))
        else:
            if not os.path.exists(ctl_file):
                raise IOError, "The specified control file does not exist."
            if verbose:
                result_code = subprocess.call([command, ctl_file])
            else:
                result_code = subprocess.call([command, ctl_file],
                    stdout=open(os.devnull, 'w'))
        if result_code > 0:
            # If baseml fails for any reason
            os.chdir(cwd)
            raise BasemlError, "BASEML has failed. Run with verbose = True to view BASEML's error message"
        if result_code < 0:
            # If the baseml process is killed by a signal somehow
            os.chdir(cwd)
            raise EnvironmentError, "The BASEML process was killed."
        if parse:
            try:
                results = read(self._rel_out_file)
            except KeyError as (errorno, strerror):
                os.chdir(cwd)
                raise KeyError, strerror
            except ValueError as (errorno, strerror):
                os.chdir(cwd)
                raise ValueError, strerror
            except:
                os.chdir(cwd)
        else:
            results = None
        os.chdir(cwd)
        return results

def read(results_file):
    """Parse a BASEML results file."""
    results = {"parameters":{}}
    count = 0
    Q_mat_found = False
    kappa_found = False
    SEs_flag = False
    branch_freqs_found = False
    num_params = -1
    with open(results_file) as results_handle:
        line = None
        while line != "":
            line = results_handle.readline()
            # Find the version number
            # Example match: 
            # "BASEML (in paml version 4.3, August 2009)  alignment.phylip"
            version_re = re.match("BASEML \(in paml version (\d+\.\d+[a-z]*).*",
                line)
            if version_re is not None:
                results["version"] = version_re.group(1)
                break
        if line == "":    
            raise ValueError, "Invalid results file"            
        for line in results_handle:
            # Find all floating point numbers in this line
            line_floats_res = re.findall("-*\d+\.\d+", line)
            line_floats = [float(val) for val in line_floats_res]
            # Find the number of branches.
            # Branches listed for parameters as pairs of numbers, ie 6..7
            branch_res = re.findall("\d+\.\.\d+", line)
            num_branches = len(branch_res)
            # Find max lnL
            # Example match:
            # ln Lmax (unconstrained) = -316.049385
            if "ln Lmax" in line and len(line_floats) == 1:
                results["lnL max"] = line_floats[0]
            # Find lnL values.
            # Example match (lnL = -2021.348300):
            # "lnL(ntime: 19  np: 22):  -2021.348300      +0.000000"
            elif "lnL(ntime:" in line and len(line_floats) > 0:
                results["lnL"] = line_floats[0]
                np_res = re.match("lnL\(ntime:\s+\d+\s+np:\s+(\d+)\)",line)
                if np_res is not None:
                    num_params = int(np_res.group(1))
            # Get parameter list. This can be useful for specifying starting
            # parameters in another run by copying the list of parameters
            # to a file called in.baseml. Since the parameters must be in
            # a fixed order and format, copying and pasting to the file is
            # best. For this reason, they are grabbed here just as a long
            # string and not as individual numbers.
            elif len(line_floats) == num_params and not SEs_flag:
                results["parameters"]["parameter list"] = line.strip()
            # Find SEs. The same format as parameters above is maintained
            # since there is a correspondance between the SE format and
            # the parameter format.
            # Example match:
            # "SEs for parameters:
            # -1.00000 -1.00000 -1.00000 801727.63247 730462.67590 -1.00000 
            elif "SEs for parameters:" in line:
                SEs_flag = True
            elif SEs_flag and len(line_floats) > 0:
                results["parameters"]["SEs"] = line.strip()
                SEs_flag = False
            # Find tree lengths.
            # Example match: "tree length =   1.71931"
            elif "tree length" in line and len(line_floats) == 1:
                results["tree length"] = line_floats[0]
            # Find the estimated tree, only taking the tree if it has
            # branch lengths
            elif re.match("\(+", line) is not None:
                if ":" in line:
                    results["tree"] = line.strip()
            # Find kappa parameter (F84, HKY85, T92 model)
            # Example match:
            # "Parameters (kappa) in the rate matrix (F84) (Yang 1994 J Mol Evol 39:105-111):
            #    3.00749"
            elif "Parameters (kappa)" in line:
                print line_floats
                kappa_found = True
            elif kappa_found and len(line_floats) > 0:
                branch_res = re.match("\s(\d+\.\.\d+)", line)
                if branch_res is None:
                    if len(line_floats) == 1:
                        results["parameters"]["kappa"] = line_floats[0]
                    else:
                        results["parameters"]["kappa"] = line_floats
                    kappa_found = False
                else:
                    if results["parameters"].get("branches") is None:
                        results["parameters"]["branches"] = {}
                    branch = branch_res.group(1)
                    if len(line_floats) > 0:
                        results["parameters"]["branches"][branch] = \
                            {"t":line_floats[0], "kappa":line_floats[1],
                            "TS":line_floats[2], "TV":line_floats[3]}
            # Find rate parameters
            # Example match: 
            # "Rate parameters:   999.00000 145.59775  0.00001  0.00001  0.00001"
            elif "Rate parameters:" in line and len(line_floats) > 0:
                results["parameters"]["rate parameters"] = line_floats
            # Find kappa under REV
            # Example match:
            # kappa under REV: 999.00000 145.76453  0.00001  0.00001  0.00001
            elif "kappa under" in line and len(line_floats) > 0:
                if len(line_floats) == 1:
                    results["parameters"]["kappa"] = line_floats[0]
                else:
                    results["parameters"]["kappa"] = line_floats
            # Find base frequencies
            # Example match:
            # "Base frequencies:   0.20090  0.16306  0.37027  0.26577"  
            elif "Base frequencies" in line and len(line_floats) > 0:
                results["parameters"]["base frequencies"] = {}
                results["parameters"]["base frequencies"]["T"] = line_floats[0]
                results["parameters"]["base frequencies"]["C"] = line_floats[1]
                results["parameters"]["base frequencies"]["A"] = line_floats[2]
                results["parameters"]["base frequencies"]["G"] = line_floats[3]
            # Find Rate matrix Q & average kappa (REV model)
            # Example match:
            # Rate matrix Q, Average Ts/Tv =   3.0308
            #  -2.483179    1.865730    0.617449    0.000000
            #   2.298662   -2.298662    0.000000    0.000000
            #   0.335015    0.000000   -0.338059    0.003044
            #   0.000000    0.000000    0.004241   -0.004241
            elif "matrix Q" in line:
                results["parameters"]["Q matrix"] = {"matrix":[]}
                if len(line_floats) > 0:
                    results["parameters"]["Q matrix"]["average Ts/Tv"] = \
                        line_floats[0]
                Q_mat_found = True
            elif Q_mat_found and len(line_floats) > 0:
                results["parameters"]["Q matrix"]["matrix"].append(line_floats)
                if len(results["parameters"]["Q matrix"]["matrix"]) == 4:
                    Q_mat_found = False
            # Find alpha
            # Example match: "alpha (gamma, K=5) = 192.47918"
            elif "alpha" in line and len(line_floats) > 0:
                results["parameters"]["alpha"] = line_floats[0]
            # Find rates
            # Example match: 
            # "rate:   0.90121  0.96051  0.99831  1.03711  1.10287"
            elif "rate: " in line and len(line_floats) > 0:
                results["parameters"]["rates"] = line_floats
            # Find frequencies
            # Example match: 
            # "freq:   0.90121  0.96051  0.99831  1.03711  1.10287"
            elif "freq: " in line and len(line_floats) > 0:
                results["parameters"]["rate frequencies"] = line_floats
            # Find branch-specific frequency parameters
            # Example match (note: I think it's possible to have 4 more
            # values per line, enclosed in brackets, so I'll account for 
            # this):
            # (frequency parameters for branches)  [frequencies at nodes] (see Yang & Roberts 1995 fig 1)

            # Node #1  ( 0.25824  0.24176  0.25824  0.24176 )
            # Node #2  ( 0.00000  0.50000  0.00000  0.50000 )
            elif "(frequency parameters for branches)" in line:
                results["parameters"]["nodes"] = {}
                branch_freqs_found = True
            elif branch_freqs_found is True:
                if len(line_floats) > 0:
                    node_res = re.match("Node \#(\d+)", line)
                    node_num = int(node_res.group(1))
                    results["parameters"]["nodes"][node_num] = {"root":False}
                    results["parameters"]["nodes"][node_num]\
                        ["frequency parameters"] = line_floats[:4]
                    if len(line_floats) > 4:
                        results["parameters"]["nodes"][node_num]\
                            ["base frequencies"] = {"T":line_floats[4],
                                                    "C":line_floats[5],
                                                    "A":line_floats[6],
                                                    "G":line_floats[7]}
                else:
                    root_res = re.match("Note: node (\d+) is root.", line)
                    if root_res is not None:
                        root_node = int(root_res.group(1))
                        results["parameters"]["nodes"][root_node]["root"] =\
                            True
                        branch_freqs_found = False
    if len(results["parameters"]) == 0:
        raise ValueError, "Invalid results file"
    return results
