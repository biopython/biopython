# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement
import os
import os.path
import re
from _paml import Paml, PamlError
import _parse_codeml

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
        self._options = {"noisy": 9, 
                        "verbose": 1, 
                        "runmode": 0,
                        "seqtype": 2, 
                        "CodonFreq": 2, 
                        "ndata": None,
                        "clock": 0, 
                        "aaDist": 0,
                        "aaRatefile": "dat/jones.dat", 
                        "model": 2,
                        "NSsites": [0], 
                        "icode": 0, 
                        "Mgene": 0,
                        "fix_kappa": 0, 
                        "kappa": 2, 
                        "fix_omega": 0,
                        "omega": .4, 
                        "fix_alpha": 1, 
                        "alpha": 0,
                        "Malpha": 0, 
                        "ncatG": 8, 
                        "getSE": 0,
                        "RateAncestor": 1, 
                        "Small_Diff": .5e-6,
                        "cleandata": 1, 
                        "fix_blength": None, 
                        "method": 0}
                        
    def write_ctl_file(self):
        """Dynamically build a CODEML control file from the options.
        
        The control file is written to the location specified by the 
        ctl_file property of the codeml class.
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
                if option[0] == "NSsites":
                    # NSsites is stored in Python as a list but in the 
                    # control file it is specified as a series of numbers
                    # separated by spaces.
                    NSsites = " ".join(["{0}".format(site) for site in option[1]])
                    ctl_handle.write("{0} = {1}\n".format(option[0], NSsites))
                else:
                    ctl_handle.write("{0} = {1}\n".format(option[0], 
                        option[1]))
    
    def read_ctl_file(self, ctl_file):
        """Parse a control file and load the options into the Codeml instance.
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
                    elif option == "NSsites":
                        site_classes = value.split(" ")
                        for n in range(len(site_classes)):
                            try:
                                site_classes[n] = int(site_classes[n])
                            except:
                                raise TypeError, \
                                    "Invalid site class: {0}".format(site_classes[n])
                        temp_options["NSsites"] = site_classes
                    elif option not in self._options:
                        raise KeyError, "Invalid option: {0}".format(option)
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
        for option in self._options.keys():
            if option in temp_options.keys():
                self._options[option] = temp_options[option]
            else:
                self._options[option] = None
                            
    def print_options(self):
        """Print out all of the options and their current settings."""
        for option in self._options.items():
            if option[0] == "NSsites":
                # NSsites is stored in Python as a list but in the 
                # control file it is specified as a series of numbers
                # separated by spaces.
                NSsites = " ".join(["{0}".format(site) for site in option[1]])
                print "{0} = {1}".format(option[0], NSsites)
            else:
                print "{0} = {1}".format(option[0], option[1])
        
    def _set_rel_paths(self):
        """Convert all file/directory locations to paths relative to the current working directory.
        
        CODEML requires that all paths specified in the control file be
        relative to the directory from which it is called rather than 
        absolute paths.
        """
        Paml._set_rel_paths(self)
        if self.tree is not None:
            self._rel_tree = os.path.relpath(self.tree, 
                self.working_dir)
        
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
    """Parse a CODEML results file."""
    results = {}
    current_model = 0
    dN_tree_flag = False
    dS_tree_flag = False
    w_tree_flag = False
    multi_models = False
    siteclass_model = None
    SEs_flag = False
    num_params = -1
    current_pair = None
    raw_aa_distances_flag = False
    ml_aa_distances_flag = False
    sequences = []
    floats_re = re.compile("-*\d+\.\d+")
    if not os.path.exists(results_file):
        raise IOError, "Results file does not exist."
    with open(results_file) as results_handle:
        lines = results_handle.readlines()
    _parse_codeml.parse_basics(lines, results)
    _parse_codeml.parse_nssites(lines, results)
    _parse_codeml.parse_pairwise(lines, results)
    _parse_codeml.parse_distances(lines, results)


           # Find site class distributions.
            # Example match 1 (normal model, 2 site classes):
            # "p:   0.77615  0.22385"
            # Example match 2 (branch site A model, 4 site classes):
            # "proportion       0.00000  0.00000  0.73921  0.26079"
            elif line[0:2] == "p:" or line[0:10] == "proportion":
                if len(line_floats) > 0: 
                    results["NSsites"][current_model]["parameters"]\
                        ["site classes"] = {}
                    for n in range(len(line_floats)):
                        results["NSsites"][current_model]["parameters"]\
                            ["site classes"][n] = {}
                        results["NSsites"][current_model]["parameters"]\
                            ["site classes"][n]["proportion"] = \
                                line_floats[n]
            # Find the omega value corresponding to each site class
            # Example match (2 site classes): "w:   0.10224  1.00000"
            elif line[0:2] == "w:":
                if len(line_floats) > 0:
                    for n in range(len(line_floats)):
                        if results["NSsites"][current_model]\
                                ["parameters"].get("site classes") \
                                is not None:
                            results["NSsites"][current_model]\
                                ["parameters"]["site classes"][n]\
                                ["omega"] = line_floats[n]
            # Find the omega values corresponding to a branch type from  
            # the clade model C for each site class
            # Example match:
            # "branch type 0:    0.31022   1.00000   0.00000"
            elif "branch type " in line:
                branch_type = re.match("branch type (\d)", line)
                if branch_type is not None:
                    branch_type_no = int(branch_type.group(1))
                    if len(line_floats) > 0:
                        for n in range(len(line_floats)):
                            if results["NSsites"][current_model]\
                                    ["parameters"].get("site classes")\
                                    is not None:
                                if results["NSsites"][current_model]\
                                        ["parameters"]["site classes"]\
                                        [n].get("branch types") is None:
                                    results["NSsites"][current_model]\
                                        ["parameters"]["site classes"]\
                                        [n]["branch types"] = {}
                                results["NSsites"][current_model]\
                                    ["parameters"]["site classes"][n]\
                                    ["branch types"][branch_type_no] =\
                                    line_floats[n]
            # Find the omega values of the foreground branch for each site
            # class in the branch site A model
            # Example match:
            # "foreground w     0.07992  1.00000 134.54218 134.54218"
            elif line[0:12] == "foreground w":
                if len(line_floats) > 0:
                    for n in range(len(line_floats)):
                        if results["NSsites"][current_model]\
                                ["parameters"].get("site classes") \
                                is not None:
                            if results["NSsites"][current_model]\
                                    ["parameters"]["site classes"]\
                                    [n].get("branch types") \
                                    is None:
                                results["NSsites"][current_model]\
                                    ["parameters"]["site classes"][n]\
                                    ["branch types"] = {}
                            results["NSsites"][current_model]\
                                ["parameters"]["site classes"][n]\
                                ["branch types"]["foreground omega"] =\
                                line_floats[n]
            # Find the omega values of the background for each site
            # class in the branch site A model
            # Example match:
            # "background w     0.07992  1.00000  0.07992  1.00000"
            elif line[0:12] == "background w":
                if len(line_floats) > 0:
                    for n in range(len(line_floats)):
                        if results["NSsites"][current_model]\
                                ["parameters"].get("site classes")\
                                is not None:
                            if results["NSsites"][current_model]\
                                    ["parameters"]["site classes"]\
                                    [n].get("branch types") is None:
                                results["NSsites"][current_model]\
                                    ["parameters"]["site classes"][n]\
                                    ["branch types"] = {}
                            results["NSsites"][current_model]\
                                ["parameters"]["site classes"][n]\
                                ["branch types"]["background omega"] =\
                                line_floats[n]
            # Get the dS, dN and omega trees for free-branch models
            elif dS_tree_flag:
                results["NSsites"][current_model]["dS tree"] = line.strip()
                dS_tree_flag = False
            elif dN_tree_flag :
                results["NSsites"][current_model]["dN tree"] = line.strip()
                dN_tree_flag = False
            elif w_tree_flag:
                line_edit = line.replace(" '#", ": ")
                line_edit = line_edit.replace("'", "")
                line_edit = line_edit.replace(" ,", ",")
                line_edit = line_edit.replace(" )", ")")
                results["NSsites"][current_model]["omega tree"] = \
                    line_edit.strip()
                w_tree_flag = False
            elif "dS tree:" in line:
                dS_tree_flag = True
            elif "dN tree:" in line:
                dN_tree_flag = True
            elif "w ratios as labels for TreeView:" in line:
                w_tree_flag = True
            elif "AA distances" in line:
                raw_aa_distances_flag = True
                # In current versions, the raw distances always come
                # first but I don't trust this to always be true
                ml_aa_distances_flag = False
            elif "ML distances of aa seqs." in line:
                ml_aa_distances_flag = True
                raw_aa_distances_flag = False
            # Find pairwise comparisons
            # Example:
            # 2 (Pan_troglo) ... 1 (Homo_sapie)
            # lnL = -291.465693
            #  0.01262 999.00000  0.00100
            #
            # t= 0.0126  S=    81.4  N=   140.6  dN/dS= 0.0010  dN= 0.0000  dS= 0.0115
            pair_res = re.match("\d+ \((.+)\) ... \d+ \((.+)\)", line)
            if pair_res is not None:
                seq1 = pair_res.group(1)
                seq2 = pair_res.group(2)
                current_pair = (seq1, seq2)
                if results.get("pairwise") is None:
                    results["pairwise"] = {}
                if results["pairwise"].get(seq1) is None:
                    results["pairwise"][seq1] = {}
                if results["pairwise"].get(seq2) is None:
                    results["pairwise"][seq2] = {}
                results["pairwise"][seq1][seq2] = {}
                results["pairwise"][seq2][seq1] = {}
            if current_pair is not None and len(line_floats) == 1:
                results["pairwise"][seq1][seq2]["lnL"] = line_floats[0]
                results["pairwise"][seq2][seq1]["lnL"] = line_floats[0]
            elif current_pair is not None and len(line_floats) > 3:
                results["pairwise"][seq1][seq2]["t"] = line_floats[0]
                results["pairwise"][seq1][seq2]["S"] = line_floats[1]
                results["pairwise"][seq1][seq2]["N"] = line_floats[2]
                results["pairwise"][seq1][seq2]["omega"] = line_floats[3]
                results["pairwise"][seq1][seq2]["dN"] = line_floats[4]
                results["pairwise"][seq1][seq2]["dS"] = line_floats[5]
                results["pairwise"][seq2][seq1]["t"] = line_floats[0]
                results["pairwise"][seq2][seq1]["S"] = line_floats[1]
                results["pairwise"][seq2][seq1]["N"] = line_floats[2]
                results["pairwise"][seq2][seq1]["omega"] = line_floats[3]
                results["pairwise"][seq2][seq1]["dN"] = line_floats[4]
                results["pairwise"][seq2][seq1]["dS"] = line_floats[5]
            # Find dN & dS for each branch, which is organized in a table
            # The existence of NaN forces me to not use the line_floats
            # method.
            # Example row (some spaces removed to make it smaller...).
            # " 6..7   0.000  167.7  54.3  0.0000  0.0000  0.0000  0.0  0.0"
            branch_res = re.match("\s+(\d+\.\.\d+)[\s+\d+\.\d+]+", line)
            if branch_res is not None and len(line_floats) > 0:
                branch = branch_res.group(1)
                if results["NSsites"][current_model]\
                        ["parameters"].get("branches") is None:
                    results["NSsites"][current_model]["parameters"]\
                        ["branches"] = {}
                params = line.strip().split()[1:]                
                results["NSsites"][current_model]["parameters"]\
                        ["branches"][branch]=\
                    {"t":float(params[0].strip()), 
                    "N":float(params[1].strip()), 
                    "S":float(params[2].strip()),  
                    "omega":float(params[3].strip()), 
                    "dN":float(params[4].strip()), 
                    "dS":float(params[5].strip()), 
                    "N*dN":float(params[6].strip()),  
                    "S*dS":float(params[7].strip())}
            # Find model parameters, which can be spread across multiple
            # lines.
            # Example matches:
            # "  p0=  0.99043  p=  0.36657 q=  1.04445
            #"  (p1=  0.00957) w=  3.25530"
            params = re.findall("([a-z]\d?)=\s+(\d+\.\d+)", line)
            if len(params) > 0:
                float_params = []
                for param in params:
                    float_params.append((param[0], float(param[1])))
                other_params = results["NSsites"][current_model].get("parameters")
                if other_params == None:
                    results["NSsites"][current_model]["parameters"] = \
                        dict(float_params)
                else:
                    results["NSsites"][current_model]["parameters"] = \
                        dict(other_params.items() + float_params)
            # Parse AA distances (raw or ML), in a lower diagonal matrix
            matrix_row_res = re.match("(.+)\s{5,15}",line)
            if matrix_row_res is not None and (raw_aa_distances_flag or \
                    ml_aa_distances_flag):
                seq_name = matrix_row_res.group(1).strip()
                if seq_name not in sequences:
                    sequences.append(seq_name)
                if results.get("distances") is None:
                    results["distances"] = {}
                if raw_aa_distances_flag:
                    if results["distances"].get("raw") is None:
                        results["distances"]["raw"] = {}
                    results["distances"]["raw"][seq_name] = {}
                    for i in range(0, len(line_floats)):
                        results["distances"]["raw"][seq_name]\
                            [sequences[i]] = line_floats[i]
                        results["distances"]["raw"][sequences[i]]\
                            [seq_name] = line_floats[i]
                else:
                    if results["distances"].get("ml") is None:
                        results["distances"]["ml"] = {}                    
                    results["distances"]["ml"][seq_name] = {}
                    for i in range(0, len(line_floats)):
                        results["distances"]["ml"][seq_name]\
                            [sequences[i]] = line_floats[i]
                        results["distances"]["ml"][sequences[i]]\
                            [seq_name] = line_floats[i]
    if len(results) == 1:
        raise ValueError, "Invalid results file"
    # The following is a rather ugly solution. If only NSsites model 0
    # is used, there is no indication in the results file that this is
    # the model being used. However, other cases, such as pairwise AA
    # results, also don't indicate any NSsites models (of course; they
    # don't use them). So, if after parsing the file the NSsites branch
    # of the results is still empty (other than the description), delete
    # it.
    if len(results["NSsites"]) == 1:
        if results["NSsites"].get(0) is not None and \
            len(results["NSsites"][0]) == 1:
                del results["NSsites"]
    return results
