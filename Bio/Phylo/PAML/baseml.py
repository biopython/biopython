# Copyright (C) 2011, 2018 by Brandon Invergo (b.invergo@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes for the support of baseml.

Maximum likelihood analysis of nucleotide sequences.
"""

import os
import os.path
from ._paml import Paml
from . import _parse_baseml


class BasemlError(EnvironmentError):
    """BASEML failed. Run with verbose=True to view BASEML's error message."""


class Baseml(Paml):
    """An interface to BASEML, part of the PAML package."""

    def __init__(self, alignment=None, tree=None, working_dir=None, out_file=None):
        """Initialize the Baseml instance.

        The user may optionally pass in strings specifying the locations
        of the input alignment and tree files, the working directory and
        the final output file.
        """
        Paml.__init__(self, alignment, working_dir, out_file)
        if tree is not None:
            if not os.path.exists(tree):
                raise FileNotFoundError("The specified tree file does not exist.")
        self.tree = tree
        self.ctl_file = "baseml.ctl"
        self._options = {
            "noisy": None,
            "verbose": None,
            "runmode": None,
            "model": None,
            "model_options": None,
            "Mgene": None,
            "ndata": None,
            "clock": None,
            "fix_kappa": None,
            "kappa": None,
            "fix_alpha": None,
            "alpha": None,
            "Malpha": None,
            "ncatG": None,
            "fix_rho": None,
            "rho": None,
            "nparK": None,
            "nhomo": None,
            "getSE": None,
            "RateAncestor": None,
            "Small_Diff": None,
            "cleandata": None,
            "icode": None,
            "fix_blength": None,
            "method": None,
        }

    def write_ctl_file(self):
        """Dynamically build a BASEML control file from the options.

        The control file is written to the location specified by the
        ctl_file property of the baseml class.
        """
        # Make sure all paths are relative to the working directory
        self._set_rel_paths()
        with open(self.ctl_file, "w") as ctl_handle:
            ctl_handle.write(f"seqfile = {self._rel_alignment}\n")
            ctl_handle.write(f"outfile = {self._rel_out_file}\n")
            ctl_handle.write(f"treefile = {self._rel_tree}\n")
            for option in self._options.items():
                if option[1] is None:
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
                        ctl_handle.write(
                            f"model = {option[1]}  {self._options['model_options']}"
                        )
                        continue
                ctl_handle.write(f"{option[0]} = {option[1]}\n")

    def read_ctl_file(self, ctl_file):
        """Parse a control file and load the options into the Baseml instance."""
        temp_options = {}
        if not os.path.isfile(ctl_file):
            raise FileNotFoundError(f"File not found: {ctl_file!r}")
        else:
            with open(ctl_file) as ctl_handle:
                for line in ctl_handle:
                    line = line.strip()
                    uncommented = line.split("*", 1)[0]
                    if uncommented != "":
                        if "=" not in uncommented:
                            raise AttributeError(
                                f"Malformed line in control file:\n{line!r}"
                            )
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
                            raise KeyError(f"Invalid option: {option}")
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
                                except ValueError:
                                    converted_value = value
                            else:
                                try:
                                    converted_value = int(value)
                                except ValueError:
                                    converted_value = value
                            temp_options[option] = converted_value
        for option in self._options:
            if option in temp_options:
                self._options[option] = temp_options[option]
            else:
                self._options[option] = None

    def _set_rel_paths(self):
        """Make file/directory paths relative to the PWD (PRIVATE).

        BASEML requires that all paths specified in the control file be
        relative to the directory from which it is called rather than
        absolute paths.
        """
        Paml._set_rel_paths(self)
        if self.tree is not None:
            self._rel_tree = os.path.relpath(self.tree, self.working_dir)

    def run(self, ctl_file=None, verbose=False, command="baseml", parse=True):
        """Run baseml using the current configuration.

        Check that the tree attribute is specified and exists, and then
        run baseml. If parse is True then read and return the result,
        otherwise return none.

        The arguments may be passed as either absolute or relative paths,
        despite the fact that BASEML requires relative paths.
        """
        if self.tree is None:
            raise ValueError("Tree file not specified.")
        if not os.path.exists(self.tree):
            raise FileNotFoundError("The specified tree file does not exist.")
        Paml.run(self, ctl_file, verbose, command)
        if parse:
            return read(self.out_file)
        return None


def read(results_file):
    """Parse a BASEML results file."""
    results = {}
    if not os.path.exists(results_file):
        raise FileNotFoundError("Results file does not exist.")
    with open(results_file) as handle:
        lines = handle.readlines()
    if not lines:
        raise ValueError(
            "Empty results file.  Did BASEML exit successfully?  "
            "Run 'Baseml.run()' with 'verbose=True'."
        )
    (results, num_params) = _parse_baseml.parse_basics(lines, results)
    results = _parse_baseml.parse_parameters(lines, results, num_params)
    if results.get("version") is None:
        raise ValueError("Invalid results file")
    return results
