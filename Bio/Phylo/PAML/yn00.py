# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import os.path
from ._paml import Paml
from . import _parse_yn00

# TODO - Restore use of with statement for closing handles automatically
# after dropping Python 2.4


class Yn00Error(EnvironmentError):
    """yn00 has failed. Run with verbose = True to view yn00's error
message"""


class Yn00(Paml):
    """This class implements an interface to yn00, part of the PAML package."""

    def __init__(self, alignment=None, working_dir=None,
                out_file=None):
        """Initialize the Yn00 instance.

        The user may optionally pass in strings specifying the locations
        of the input alignment, the working directory and
        the final output file.
        """
        Paml.__init__(self, alignment, working_dir, out_file)
        self.ctl_file = "yn00.ctl"
        self._options = {"verbose": None,
                        "icode": None,
                        "weighting": None,
                        "commonf3x4": None,
                        "ndata": None}

    def write_ctl_file(self):
        """Dynamically build a yn00 control file from the options.

        The control file is written to the location specified by the
        ctl_file property of the yn00 class.
        """
        # Make sure all paths are relative to the working directory
        self._set_rel_paths()
        with open(self.ctl_file, 'w') as ctl_handle:
            ctl_handle.write("seqfile = %s\n" % self._rel_alignment)
            ctl_handle.write("outfile = %s\n" % self._rel_out_file)
            for option in self._options.items():
                if option[1] is None:
                    # If an option has a value of None, there's no need
                    # to write it in the control file; it's normally just
                    # commented out.
                    continue
                ctl_handle.write("%s = %s\n" % (option[0], option[1]))

    def read_ctl_file(self, ctl_file):
        """Parse a control file and load the options into the yn00 instance.
        """
        temp_options = {}
        if not os.path.isfile(ctl_file):
            raise IOError("File not found: %r" % ctl_file)
        else:
            with open(ctl_file) as ctl_handle:
                for line in ctl_handle:
                    line = line.strip()
                    uncommented = line.split("*", 1)[0]
                    if uncommented != "":
                        if "=" not in uncommented:
                            raise AttributeError(
                                "Malformed line in control file:\n%r" % line)
                        (option, value) = uncommented.split("=")
                        option = option.strip()
                        value = value.strip()
                        if option == "seqfile":
                            self.alignment = value
                        elif option == "outfile":
                            self.out_file = value
                        elif option not in self._options:
                            raise KeyError("Invalid option: %s" % option)
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
        for option in self._options:
            if option in temp_options:
                self._options[option] = temp_options[option]
            else:
                self._options[option] = None

    def run(self, ctl_file=None, verbose=False, command="yn00",
                parse=True):
        Paml.run(self, ctl_file, verbose, command)
        if parse:
            results = read(self.out_file)
        else:
            results = None
        return results


def read(results_file):
    """Parse a yn00 results file."""
    results = {}
    if not os.path.exists(results_file):
        raise IOError("Results file does not exist.")
    with open(results_file) as handle:
        lines = handle.readlines()
    for line_num in range(len(lines)):
        line = lines[line_num]
        if "(A) Nei-Gojobori (1986) method" in line:
            ng86_start = line_num + 1
        elif "(B) Yang & Nielsen (2000) method" in line:
            (results, sequences) = _parse_yn00.parse_ng86(lines[ng86_start:line_num],
                    results)
            yn00_start = line_num + 1
        elif "(C) LWL85, LPB93 & LWLm methods" in line:
            results = _parse_yn00.parse_yn00(lines[yn00_start:line_num], results,
                    sequences)
            results = _parse_yn00.parse_others(lines[line_num+1:], results,
                    sequences)
    if len(results) == 0:
        raise ValueError("Invalid results file.")
    return results
