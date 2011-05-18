# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import with_statement
import os.path
from _paml import Paml, PamlError
import _parse_yn00

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
    if not os.path.exists(results_file):
        raise IOError, "Results file does not exist."
    results = _parse_yn00.parse(results_file)
    return results
