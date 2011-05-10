# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement
import os
import re
from _paml import Paml, PamlError

class Yn00Error(EnvironmentError):
    """yn00 has failed. Run with verbose = True to view yn00's error
message"""

class Yn00(Paml):
    """This class implements an interface to yn00, part of the PAML package."""

    def __init__(self, alignment = None, working_dir = None,
                out_file = None):
        """Initialize the Yn00 instance. 
        
        The user may optionally pass in strings specifying the locations
        of the input alignment, the working directory and
        the final output file. 
        """
        Paml.__init__(self, alignment, working_dir, out_file)
        self.ctl_file = "yn00.ctl"
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
                
    def run(self, ctl_file = None, verbose = False, command = "yn00",
                parse = True):
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
