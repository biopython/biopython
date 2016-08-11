# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.
# Revisions copyright 2014 by Melissa Gymrek <mgymrek@mit.edu>.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This module allows you to control Simcoal2 and FastSimcoal (DEPRECATED)."""

import os
import sys
from Bio.Application import AbstractCommandline, _Option, _Switch


class SimCoalController(object):
    def __init__(self, simcoal_dir):
        """Initializes the controller. (DEPRECATED)

        simcoal_dir is the directory where simcoal is.

        The initializer checks for existence and executability of binaries.
        """
        self.simcoal_dir = simcoal_dir
        self.os_name = os.name  # remove this?
        dir_contents = os.listdir(self.simcoal_dir)
        # We expect the tool to be installed as simcoal2(.exe)
        # without any trailing version number.
        self.bin_name = "simcoal2"
        if self.bin_name not in dir_contents:
            # Try case insensitive,
            dir_contents = [x.lower() for x in dir_contents]
        if self.bin_name not in dir_contents:
            # Try with .exe
            self.bin_name += '.exe'
        if self.bin_name not in dir_contents:
            raise IOError("SimCoal not available")
        if not os.access(os.path.join(self.simcoal_dir, self.bin_name),
                         os.X_OK):
            raise IOError("SimCoal not executable")

    def run_simcoal(self, par_file, num_sims, ploydi='1', par_dir='.'):
        """Executes SimCoal.
        """
        if par_dir is None:
            par_dir = os.sep.join([".", 'SimCoal', 'runs'])
        curr_dir = os.getcwd()
        # TODO - Make sure we change drive on Windows as well?
        os.chdir(par_dir)
        exe = os.path.join(self.simcoal_dir, self.bin_name)
        if " " in exe:
            exe = '"' + exe + '"'
        cmd = exe + ' ' + par_file + ' ' + str(num_sims) + ' ' + ploydi
        # TODO - Better way to spot if on Jython on Windows?
        if sys.platform == "win32" or self.bin_name.endswith(".exe"):
            # There is no /dev/nul on Windows
            cmd += ' > nul 2>nul'
        else:
            cmd += ' >/dev/null 2>&1'
        os.system(cmd)
        os.chdir(curr_dir)


class _FastSimCoalCommandLine(AbstractCommandline):
    """ Command Line Wrapper for Fastsimcoal
    """
    def __init__(self, fastsimcoal_dir=None, cmd='fastsimcoal', **kwargs):
        self.parameters = [
            _Option(["-i", "--ifile", "parfile"], "Name of the parameter file",
                    filename=True, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-n", "--numsims", "numsims"], "Number of simulations to perform",
                    filename=False, equate=False, is_required=True,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-t", "--tfile", "tfile"], "Name of template parameter file",
                    filename=True, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-f", "--dfile", "dfile"], "Name of parameter definition file",
                    filename=True, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-F", "--dFile", "dFile"],
                    """Same as -f but only uses simple parameters defined
                    in the template file. Complex params are recomputed""",
                    filename=True, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-e", "--efile", "efile"],
                    """Parameter prior definition file.
                    Parameters drawn from specified distributions are
                    substituted into template file.""",
                    filename=True, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Option(["-E", "--numest", "numest"],
                    """Number of estimations from parameter priors.
                    Listed parameter values are substituted in template file.""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-g", "--genotypic", "genotypic"], "Generates Arlequin projects with genotypic data"),
            _Switch(["-p", "--phased", "phased"], "Specifies that phase is known in Arlequin output"),
            _Option(["-s", "--dnatosnp", "dnatosnp"],
                    """"Output DNA as SNP data (0: ancestral, 1: derived
                    and specify maximum no. SNPs to output.""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-S", "--allsites", "allsites"],
                    """Output the whole DNA sequence, including monomorphic sites"""),
            _Switch(["-I", "--inf", "inf"],
                    """Generates DNA mutations according to an
                    infinite sites (IS) mutation model."""),
            _Switch(["-d", "--dsfs", "dsfs"], "Computes derived site frequency spectrum"),
            _Switch(["-m", "--msfs", "msfs"], "Computes minor site frequency spectrum"),
            _Option(["-o", "--oname", "oname"], "Generic name for observed SFS files",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, str)),
            _Switch(["-H", "--header", "header"], "Generates header in site frequency spectrum files."),
            _Switch(["-q", "--quiet", "quiet"], "Minimal messages output to console"),
            _Switch(["-T", "--tree", "tree"], "Output coalescent tree in nexus format."),
            _Option(["-k", "--keep", "keep"],
                    """Number of simulated polymorphic sites kept in memory.
                    If the simulated no. is larger, then temporary files are created.""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["--seed", "seed"], "Seed for the random number generator (positive int <=1E6)",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-x", "--noarloutput", "noarloutput"], "Does not generate Arlequin output"),
            _Switch(["-D", "--dadioutput", "dadioutput"], "Output SFS in dadi format"),
            _Option(["-M", "--maxlhood", "maxlhood"],
                    """Perform parameter estimation by max lhood from SFS, and
                    define stop criterion as min., rel., diff. in parameter
                    values between iterations""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, float)),
            _Option(["-N", "--maxnumsims", "maxnumsims"],
                    """Maximum number of simulations to perform during
                    likelihood maximization.""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-l", "--minnumloops", "minnumloops"],
                    """Minimum number of iteration loops to perform during
                    likelihood maximization.""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-L", "--maxnumloops", "maxnumloops"],
                    """Maximum number of iterations to perform during
                    likelihood maximization""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-C", "--minSFSCount", "minSFSCount"],
                    """Minimum observed SFS entry count taken into account
                    in likelihood computation""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-0", "--removeZeroSFS", "removeZeroSFS"],
                    """Do not take into account monomorphic sites for
                    SFS likelihood computation."""),
            _Option(["-a", "--ascDeme", "ascDeme"],
                    """This is the deme id where ascertainment is performed
                    when simulating SNPs.""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Option(["-A", "--ascSize", "ascSize"],
                    """Number of ascertained chromosomes used to define SNPs in
                    a given deme.""",
                    filename=False, equate=False, is_required=False,
                    checker_function=lambda x: isinstance(x, int)),
            _Switch(["-u", "--multiSFS", "multiSFS"],
                    "Generate or use multidimensional SFS")]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class FastSimCoalController(object):
    def __init__(self, fastsimcoal_dir=None, bin_name="fsc252"):
        """Initializes the controller.

        fastsimcoal_dir is the directory where fastsimcoal is.
        By default the binary should be called fsc252.
        bin_name specifies a different name for the binary.

        The initializer checks for existence and executability of binaries
        and sets up the command line controller.

        Fastsimcoal2 is available here: http://cmpg.unibe.ch/software/fastsimcoal2/.
        This wrapper was written and tested for fastsimcoal version 2.51.
        """
        self.bin_name = bin_name
        self.fastsimcoal_dir = fastsimcoal_dir
        if fastsimcoal_dir is None:
            for path in os.environ["PATH"].split(os.pathsep):
                if os.path.isfile(os.path.join(path, self.bin_name)):
                    self.fastsimcoal_dir = path
            if self.fastsimcoal_dir is None:
                raise IOError("Fastsimcoal not available")
        else:
            dir_contents = os.listdir(fastsimcoal_dir)
            if self.bin_name not in dir_contents:
                raise IOError("Fastsimcoal not available")
        if not os.access(os.path.join(self.fastsimcoal_dir, self.bin_name), os.X_OK):
            raise IOError("Fastsimcoal not executable")

    def run_fastsimcoal(self, par_file, num_sims, par_dir='.', opts=None):
        """Executes Fastsimcoal.

        par_file is the input parameter file (--ifile) for fastsimcoal.
        num_sims is the number of simulations to perform.
        par_dir is the directory where par_file is and where output will be written.
        opts is a dictionary of additional options to fastsimcoal.
        """
        if opts is None:
            opts = {}
        if par_dir is None:
            par_dir = os.sep.join([".", "Fastsimcoal", "runs"])
            if not os.path.exists(par_dir):
                os.mkdir(par_dir)
        curr_dir = os.getcwd()
        os.chdir(par_dir)
        if par_file is None:  # Must use .tpl for -t instead if no par_file
            controller = _FastSimCoalCommandLine(cmd=os.path.join(self.fastsimcoal_dir, self.bin_name),
                                                 numsims=num_sims, **opts)
        else:
            controller = _FastSimCoalCommandLine(cmd=os.path.join(self.fastsimcoal_dir, self.bin_name),
                                                 parfile=par_file, numsims=num_sims, **opts)
        controller()
        os.chdir(curr_dir)
