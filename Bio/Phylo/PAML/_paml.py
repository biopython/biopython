# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Base class for the support of PAML, Phylogenetic Analysis by Maximum Likelihood."""


import os
import subprocess


class PamlError(EnvironmentError):
    """paml has failed.

    Run with verbose=True to view the error message.
    """


class Paml:
    """Base class for wrapping PAML commands."""

    def __init__(self, alignment=None, working_dir=None, out_file=None):
        """Initialize the class."""
        if working_dir is None:
            self.working_dir = os.getcwd()
        else:
            self.working_dir = working_dir
        if alignment is not None:
            if not os.path.exists(alignment):
                raise FileNotFoundError("The specified alignment file does not exist.")
        self.alignment = alignment
        self.out_file = out_file
        self._options = {}  # will be set in subclasses

    def write_ctl_file(self):
        """Write control file."""
        pass

    def read_ctl_file(self):
        """Read control file."""
        pass

    def print_options(self):
        """Print out all of the options and their current settings."""
        for option in self._options.items():
            print("%s = %s" % (option[0], option[1]))

    def set_options(self, **kwargs):
        """Set the value of an option.

        This function abstracts the options dict to prevent the user from
        adding options that do not exist or misspelling options.
        """
        for option, value in kwargs.items():
            if option not in self._options:
                raise KeyError("Invalid option: " + option)
            else:
                self._options[option] = value

    def get_option(self, option):
        """Return the value of an option."""
        if option not in self._options:
            raise KeyError("Invalid option: " + option)
        else:
            return self._options.get(option)

    def get_all_options(self):
        """Return the values of all the options."""
        return list(self._options.items())

    def _set_rel_paths(self):
        """Convert all file/directory locations to paths relative to the current working directory (PRIVATE).

        paml requires that all paths specified in the control file be
        relative to the directory from which it is called rather than
        absolute paths.
        """
        if self.working_dir is not None:
            self._rel_working_dir = os.path.relpath(self.working_dir)
        if self.alignment is not None:
            self._rel_alignment = os.path.relpath(self.alignment, self.working_dir)
        if self.out_file is not None:
            self._rel_out_file = os.path.relpath(self.out_file, self.working_dir)

    def run(self, ctl_file, verbose, command):
        """Run a paml program using the current configuration.

        Check that the class attributes exist and raise an error
        if not. Then run the command and check if it succeeds with
        a return code of 0, otherwise raise an error.

        The arguments may be passed as either absolute or relative
        paths, despite the fact that paml requires relative paths.
        """
        if self.alignment is None:
            raise ValueError("Alignment file not specified.")
        if not os.path.exists(self.alignment):
            raise FileNotFoundError("The specified alignment file does not exist.")
        if self.out_file is None:
            raise ValueError("Output file not specified.")
        if self.working_dir is None:
            raise ValueError("Working directory not specified.")
        # Get the current working directory
        cwd = os.getcwd()
        # Move to the desired working directory
        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)
        os.chdir(self.working_dir)
        # If no external control file was specified...
        if ctl_file is None:
            # Dynamically build a control file
            self.write_ctl_file()
            ctl_file = self.ctl_file
        else:
            if not os.path.exists(ctl_file):
                raise FileNotFoundError("The specified control file does not exist.")
        if verbose:
            result_code = subprocess.call([command, ctl_file])
        else:
            with open(os.devnull) as dn:
                result_code = subprocess.call([command, ctl_file], stdout=dn, stderr=dn)
        os.chdir(cwd)
        if result_code > 0:
            # If the program fails for any reason
            raise PamlError(
                "%s has failed (return code %i). Run with verbose = True to view error message"
                % (command, result_code)
            )
        if result_code < 0:
            # If the paml process is killed by a signal somehow
            raise OSError(
                "The %s process was killed (return code %i)." % (command, result_code)
            )
