# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement
import os
import os.path
import subprocess
import re

class Yn00Error(EnvironmentError):
    """yn00 has failed. Run with verbose = True to view yn00's error
message"""

class Yn00(object):
    """This class implements an interface to yn00, part of the PAML package."""

    def __init__(self, alignment = None, working_dir = None,
                out_file = None):
        """Initialize the Yn00 instance. 
        
        The user may optionally pass in strings specifying the locations
        of the input alignment, the working directory and
        the final output file. 
        """
        self.ctl_file = "yn00.ctl"
        if working_dir is None:
            self.working_dir = os.getcwd()
        else:
            self.working_dir = working_dir
        if alignment is not None:
            if not os.path.exists(alignment):
                raise IOError, "The specified alignment file does not exist."
        self.alignment = alignment
        self.out_file = out_file
        self._options = {"verbose": 0,
                        "icode": 0,
                        "weighting": 0,
                        "commonf3x4": 0,
                        "ndata": None}

    def write_ctl_file(self):
        """Dynamically build a yn00 control file from the options.
        
        The control file is written to the location specified by the 
        ctl_file property of the yn00 class.
        """
        # Make sure all paths are relative to the working directory
        try:
            self._set_rel_paths()
        except (AttributeError, ValueError) as error:
            raise error
        with open(self.ctl_file, 'w') as ctl_handle:
            ctl_handle.write("seqfile = {0}\n".format(self._rel_alignment))
            ctl_handle.write("outfile = {0}\n".format(self._rel_out_file))
            for option in self._options.items():
                if option[1] == None:
                    # If an option has a value of None, there's no need
                    # to write it in the control file; it's normally just
                    # commented out.
                    continue
                ctl_handle.write("{0} = {1}\n".format(option[0], 
                    option[1]))
    
    def read_ctl_file(self, ctl_file):
        """Parse a control file and load the options into the yn00 instance.
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
                    elif option == "outfile":
                        self.out_file = value
                    elif option not in self._options:
                        raise KeyError, "Invalid option: {0}".format(option)
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
        
        yn00 requires that all paths specified in the control file be
        relative to the directory from which it is called rather than 
        absolute paths.
        """
        if self.working_dir is not None:
            self._rel_working_dir = os.path.relpath(self.working_dir)
        if self.alignment is not None:
            self._rel_alignment = os.path.relpath(self.alignment, 
                self.working_dir)
        if self.out_file is not None:
            self._rel_out_file = os.path.relpath(self.out_file, 
                self.working_dir)
        
    def run(self, ctl_file = None, verbose = False, command = "yn00",
                parse = True):
        """Run cyn00 using the current configuration and then parse the results. 
        
        Return a process signal so the user can determine if
        the execution was successful (return code 0 is successful, -N
        indicates a failure). The arguments may be passed as either 
        absolute or relative paths, despite the fact that yn00 
        requires relative paths.
        """
        if self.alignment is None:
            raise ValueError, "Alignment file not specified."
        if not os.path.exists(self.alignment):
            raise IOError, "The specified alignment file does not exist."
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
            # If yn00 fails for any reason
            os.chdir(cwd)
            raise Yn00Error, "yn00 has failed. Run with verbose = True to view yn00's error message"
        if result_code < 0:
            # If the yn00 process is killed by a signal somehow
            os.chdir(cwd)
            raise EnvironmentError, "The yn00 process was killed."
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
    """Parse a yn00 results file."""
    results = {}
    sequences = []
    current_method = 0
    with open(results_file) as results_handle:          
        for line in results_handle:
            # Find all floating point numbers in this line
            line_floats_res = re.findall("-*\d+\.\d+", line)
            line_floats = [float(val) for val in line_floats_res]
            # The results file is organized by method
            if "(A) Nei-Gojobori (1986) method" in line:
                current_method = 1
                continue
            elif "(B) Yang & Nielsen (2000) method" in line:
                current_method = 2
                continue
            elif "(C) LWL85, LPB93 & LWLm methods" in line:
                current_method = 3
                continue
            if current_method == 1:
                # Nei_Gojobori results are organized in a lower 
                # triangular mattrix, with the sequence names labeling
                # the rows and statistics in the format:
                # w (dN dS) per column
                # Example row (2 columns):
                # 0.0000 (0.0000 0.0207) 0.0000 (0.0000 0.0421)
                matrix_row_res = re.match("(.+)\s{5,15}",line)
                if matrix_row_res is not None:
                    seq_name = matrix_row_res.group(1).strip()
                    sequences.append(seq_name)
                    results[seq_name] = {}
                    for i in range(0, len(line_floats), 3):
                        NG86 = {}
                        NG86["omega"] = line_floats[i]
                        NG86["dN"] = line_floats[i+1]
                        NG86["dS"] = line_floats[i+2]
                        results[seq_name][sequences[i/3]] = {"NG86":NG86}
                        results[sequences[i/3]][seq_name] = {"NG86":NG86}
            elif current_method == 2:
                # Yang & Nielsen results are organized in a table with
                # each row comprising one pairwise species comparison.
                # Rows are labeled by spequence number rather than by
                # sequence name.
                # Example (header row and first table row):
                # seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE
                # 2    1    67.3   154.7   0.0136  3.6564  0.0000 -0.0000 +- 0.0000  0.0150 +- 0.0151
                row_res = re.match("\s+(\d+)\s+(\d+)", line)
                if row_res is not None:
                    seq1 = int(row_res.group(1))
                    seq2 = int(row_res.group(2))
                    seq_name1 = sequences[seq1-1]
                    seq_name2 = sequences[seq2-1]
                    YN00 = {}
                    YN00["S"] = line_floats[0]
                    YN00["N"] = line_floats[1]
                    YN00["t"] = line_floats[2]
                    YN00["kappa"] = line_floats[3]
                    YN00["omega"] = line_floats[4]
                    YN00["dN"] = line_floats[5]
                    YN00["dN SE"] = line_floats[6]
                    YN00["dS"] = line_floats[7]
                    YN00["dS SE"] = line_floats[8]
                    results[seq_name1][seq_name2]["YN00"] = YN00
                    results[seq_name2][seq_name1]["YN00"] = YN00
                seq_name1 = None
                seq_name2 = None
            elif current_method == 3:
                # The remaining methods are grouped together. Statistics
                # for all three are listed for each of the pairwise 
                # species comparisons, with each method's results on its
                # own line.
                # The stats in this section must be handled differently
                # due to the possible presence of NaN values, which won't
                # get caught by my typical "line_floats" method used above.
                # Example:
                # 2 (Pan_troglo) vs. 1 (Homo_sapie)

                # L(i):      143.0      51.0      28.0  sum=    222.0
                # Ns(i):    0.0000    1.0000    0.0000  sum=   1.0000
                # Nv(i):    0.0000    0.0000    0.0000  sum=   0.0000
                # A(i):     0.0000    0.0200    0.0000
                # B(i):    -0.0000   -0.0000   -0.0000
                # LWL85:  dS =  0.0227 dN =  0.0000 w = 0.0000 S =   45.0 N =  177.0
                # LWL85m: dS =    -nan dN =    -nan w =   -nan S =   -nan N =   -nan (rho = -nan)
                # LPB93:  dS =  0.0129 dN =  0.0000 w = 0.0000
                comp_res = re.match("\d+ \((.+)\) vs. \d+ \((.+)\)", line)
                if comp_res is not None:
                    seq_name1 = comp_res.group(1)
                    seq_name2 = comp_res.group(2)
                elif seq_name1 is not None and seq_name2 is not None:
                    if "dS =" in line:
                        stats = {}
                        line_stats = line.split(":")[1].strip()
                        stats_split = line_stats.split()
                        for i in range(0, len(stats_split), 3):
                            stat = stats_split[i].strip("()")
                            value = stats_split[i+2].strip("()")
                            try:
                                stats[stat] = float(value)
                            except:
                                stats[stat] = None
                        if "LWL85:" in line:
                            results[seq_name1][seq_name2]["LWL85"] = stats
                            results[seq_name2][seq_name1]["LWL85"] = stats
                        elif "LWL85m" in line:
                            results[seq_name1][seq_name2]["LWL85m"] = stats
                            results[seq_name2][seq_name1]["LWL85m"] = stats
                        elif "LPB93" in line:
                            results[seq_name1][seq_name2]["LPB93"] = stats
                            results[seq_name2][seq_name1]["LPB93"] = stats
    if len(results) == 0:
        raise ValueError, "Invalid results file"
    return results
