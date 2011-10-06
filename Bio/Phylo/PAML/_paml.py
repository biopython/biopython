# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import os
import subprocess

try:
    from os.path import relpath as _relpath
except ImportError:
    #New in Python 2.6
    def _relpath(path, start=None):
        """Return a relative version of a path.

        Implementation by James Gardner in his BareNecessities
        package, under MIT licence.

        With a fix for Windows where posixpath.sep (and functions like
        join) use the Unix slash not the Windows slash.
        """
        import posixpath
        if start is None:
            start = posixpath.curdir
        else:
            start = start.replace(os.path.sep, posixpath.sep)
        if not path:
            raise ValueError("no path specified")
        else:
            path = path.replace(os.path.sep, posixpath.sep)
        start_list = posixpath.abspath(start).split(posixpath.sep)
        path_list = posixpath.abspath(path).split(posixpath.sep)
        # Work out how much of the filepath is shared by start and path.
        i = len(posixpath.commonprefix([start_list, path_list]))
        rel_list = [posixpath.pardir] * (len(start_list)-i) + path_list[i:]
        if not rel_list:
            return posixpath.curdir.replace(posixpath.sep, os.path.sep)
        return posixpath.join(*rel_list).replace(posixpath.sep, os.path.sep)

class PamlError(EnvironmentError):
    """paml has failed. Run with verbose = True to view the error
message"""

class Paml(object):
    
    def __init__(self, alignment = None, working_dir = None,
                out_file = None):
        if working_dir is None:
            self.working_dir = os.getcwd()
        else:
            self.working_dir = working_dir
        if alignment is not None:
            if not os.path.exists(alignment):
                raise IOError, "The specified alignment file does not exist."
        self.alignment = alignment
        self.out_file = out_file
        
    def write_ctl_file(self):
        pass
        
    def read_ctl_file(self):
        pass
        
    def print_options(self):
        """Print out all of the options and their current settings."""
        for option in self._options.items():
            print "%s = %s" % (option[0], option[1])
 
    def set_options(self, **kwargs):
        """Set the value of an option. 
        
        This function abstracts the options dict to prevent the user from 
        adding options that do not exist or mispelling options.
        """
        for option, value in kwargs.items():
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
        
        paml requires that all paths specified in the control file be
        relative to the directory from which it is called rather than 
        absolute paths.
        """
        if self.working_dir is not None:
            self._rel_working_dir = _relpath(self.working_dir)
        if self.alignment is not None:
            self._rel_alignment = _relpath(self.alignment, 
                self.working_dir)
        if self.out_file is not None:
            self._rel_out_file = _relpath(self.out_file, self.working_dir)
        
    def run(self, ctl_file, verbose, command):
        """Run a paml program using the current configuration and then parse the results. 
        
        Return a process signal so the user can determine if
        the execution was successful (return code 0 is successful, -N
        indicates a failure). The arguments may be passed as either 
        absolute or relative paths, despite the fact that paml 
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
            if verbose:
                result_code = subprocess.call([command, self.ctl_file])
            else:
                # To suppress output, redirect it to a pipe to nowhere
                result_code = subprocess.call([command, self.ctl_file],
                    stdout=subprocess.PIPE)
        else:
            if not os.path.exists(ctl_file):
                raise IOError, "The specified control file does not exist."
            if verbose:
                result_code = subprocess.call([command, ctl_file])
            else:
                result_code = subprocess.call([command, ctl_file],
                    stdout=subprocess.PIPE)
        os.chdir(cwd)
        if result_code > 0:
            # If the program fails for any reason
            raise PamlError, \
            "%s has failed (return code %i). Run with verbose = True to view error message" \
            % (command, result_code)
        if result_code < 0:
            # If the paml process is killed by a signal somehow
            raise EnvironmentError, "The %s process was killed (return code %i)." \
                  % (command, result_code)
