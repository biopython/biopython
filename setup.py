"""Distutils based setup script for Biopython.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install

For more in-depth instructions, see the installation section of the
biopython manual, linked to from:

http://biopython.org/wiki/html/BioPython/BiopythonCode.html

Or for more details about the options available from distutils, look at
the 'Installing Python Modules' distutils documentation, available from:

http://python.org/sigs/distutils-sig/doc/

Or, if all else fails, feel free to write to the biopython list
at biopython@biopython.org and ask for help.
"""

import sys
import os
try:
    from distutils.core import setup
    from distutils.command.install import install
    from distutils.core import Command
except ImportError:
    print "Biopython installation requires distutils, available with python 2.0"
    print "or better, or from:"
    print "http://python.org/sigs/distutils-sig/download.html"
    sys.exit(0)

# check if the distutils has the new extension class stuff
try:
    from distutils.extension import Extension
except ImportError:
    print "Your version of distutils is really old. You need to upgrade"
    print "to a newer version. The latest releases of distutils are available"
    print "from http://python.org/sigs/distutils-sig/download.html"
    sys.exit(0)

# --- check for installed programs needed by Biopython

def check_install(name, check_library, location, other_messages = None):
    """Check if a program is installed and print a warning message if not.

    This helps users at least know they are missing some installed stuff
    and where to get it when they install biopython.

    Arguments:
    
    o check_library -- a function to check whether or not the specified
    program and version is present, returns 1 if it is, 0 otherwise.

    o name -- the name of the library we are looking for

    o location -- a URL where the library can be downloaded

    o other_messages -- other random messages to print if the library
    is not present (ie. version information, etc...)
    """
    if not(check_library()):
        print "\nWARNING -- %s is not installed." % name
        print "You should install this from:"
        print location
        print "to get full functionality from Biopython."
        if other_messages:
            print other_messages

# -- functions to check for specific libraries and versions.

def check_Martel():
    """Check for Martel version 0.5 or better.
    """
    try:
        import Martel
        if str(Martel.__version__) >= "0.5":
            return 1
        else:
            return 0
    except ImportError:
        return 0

def check_mxTextTools():
    try:
        import TextTools
        return 1
    except ImportError:
        try:
            from mx import TextTools
            return 1
        except ImportError:
            return 0

def check_Numpy():
    try:
        import Numeric
        return 1
    except ImportError:
        return 0

def check_reportlab():
    try:
        import reportlab
        return 1
    except ImportError:
        return 0
                  
class my_install(install):
    """Override the standard install to check for dependencies.

    This will just run the normal install, and then print warning messages
    if packages are missing.
    """
    def run(self):
        # run the normal install and everthing
        install.run(self)

        # now print warning messages if we are missing stuff
        check_install("Martel", check_Martel,
                      "http://www.biopython.org/~dalke/Martel/",
                      "Version 0.5 or better required")

        check_install("mxTextTools", check_mxTextTools,
                      "http://www.lemburg.com/files/python/mxExtensions.html")

        check_install("Numerical Python", check_Numpy,
                      "http://numpy.sourceforge.net/")

        check_install("Reportlab", check_reportlab,
                      "http://www.reportlab.com/download.html")

class run_tests(Command):
    """Run all of the tests for the package.

    This is a automatic test run class to make distutils kind of act like
    perl. With this you can do:

    python setup.py build
    python setup.py install
    python setup.py test
    """
    description = "Automatically run the test suite for the package."

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("Tests")
        import run_tests
        run_tests.main([])

        # change back to the current directory
        os.chdir(this_dir)

# --- set up the packages we are going to install
# standard biopython packages
biopython_packages = ['Bio',
                      'Bio.Align',
                      'Bio.Alphabet',
                      'Bio.Application',
                      'Bio.Blast',
                      'Bio.Clustalw',
                      'Bio.Data',
                      'Bio.Encodings',
                      'Bio.Emboss',
                      'Bio.Enzyme',
                      'Bio.FSSP',
                      'Bio.Fasta',
                      'Bio.GA',
                      'Bio.GA.Crossover',
                      'Bio.GA.Mutation',
                      'Bio.GA.Repair',
                      'Bio.GA.Selection',
                      'Bio.GenBank',
                      'Bio.Gobase',
                      'Bio.Graphics',
                      'Bio.HMM',
                      'Bio.InterPro',
                      'Bio.IntelliGenetics',
                      'Bio.Kabat',
                      'Bio.KEGG',
                      'Bio.KEGG.Compound',
                      'Bio.KEGG.Enzyme',
                      'Bio.KEGG.Map',
                      'Bio.Medline',
                      'Bio.MetaTool',
                      'Bio.NBRF',
                      'Bio.NeuralNetwork',
                      'Bio.NeuralNetwork.BackPropagation',
                      'Bio.NeuralNetwork.Gene',
                      'Bio.Pathway',
                      'Bio.Pathway.Rep',
                      'Bio.PDB',
                      'Bio.Prosite',
                      'Bio.Rebase',
                      'Bio.SCOP',
                      'Bio.SeqIO',
                      'Bio.SubsMat',
                      'Bio.SwissProt',
                      'Bio.Tools',
                      'Bio.Tools.Classification',
                      'Bio.Tools.Clustering',
                      'Bio.Tools.MultiProc',
                      'Bio.Tools.Parsers',
                      'Bio.UniGene',
                      'Bio.WWW'
                      ]

# Martel -- for building parsers
martel_packages = ['Martel']
# a flag to determine if we should install Martel
INSTALL_MARTEL = 1

# -- setup the packages list
all_packages = biopython_packages
# only install Martel if we have the flag set 
if INSTALL_MARTEL:
    all_packages.extend(martel_packages)

setup(name='biopython', 
      version='1.00a4',
      author='The Biopython Consortium',
      author_email='biopython@biopython.org',
      url='http://www.biopython.org/',

      cmdclass = {"install" : my_install,
                  "test" : run_tests},
      
      packages = all_packages,
      
      ext_modules = [Extension('Bio.Tools.Classification.cSVM',
                               ['Bio/Tools/Classification/cSVMmodule.c']
                               ),
                     Extension('Bio.Tools.clistfns',
                               ['Bio/Tools/clistfnsmodule.c']
                               ),
                     Extension('Bio.Tools.cmathfns',
                               ['Bio/Tools/cmathfnsmodule.c']
                               ),
                     Extension('Bio.Tools.cstringfns',
                               ['Bio/Tools/cstringfnsmodule.c']
                               ),
                     Extension('Bio.Align.csupport',
                               ['Bio/Align/csupportmodule.c']
                               ),
                     Extension('Bio.Align.cfastpairwise',
                               ['Bio/Align/cfastpairwisemodule.c']
                               )
                     ]
      )

