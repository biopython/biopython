"""Distutils based setup script for Biopython.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install

For more in-depth instructions, see the installation section of the
biopython manual, linked to from:

http://www.biopython.org/documentation/

Or for more details about the options available from distutils, look at
the 'Installing Python Modules' distutils documentation, available from:

http://python.org/sigs/distutils-sig/doc/

Or, if all else fails, feel free to write to the biopython list at
biopython@biopython.org and ask for help.
"""
import sys
import os
import shutil

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 2):
    print "Biopython requires Python 2.2.  Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)

from distutils.core import setup
from distutils.core import Command
from distutils.command.install import install
from distutils.command.build_py import build_py
from distutils.extension import Extension

def get_yes_or_no(question, default):
    if default:
        option_str = "(Y/n)"
        default_str = 'y'
    else:
        option_str = "(y/N)"
        default_str = 'n'

    while 1:
        print "%s %s " % (question, option_str),
        response = raw_input().lower()
        if not response:
            response = default_str
        if response[0] in ['y', 'n']:
            break
        print "Please answer y or n."
    return response[0] == 'y'

class install_biopython(install):
    """Override the standard install to check for dependencies.

    This will just run the normal install, and then print warning messages
    if packages are missing.

    This also has a hack to install DTDs needed for EUtils into
    Bio.EUtils.DTDs. This is not a pretty thing since we should 
    really only have pure python modules installed.
    """
    def check_dependencies(self):
        """S.check_dependencies() -> boolean

        Return whether the installation should continue.

        """
        # There should be some way for the user to tell specify not to
        # check dependencies.  For example, it probably should not if
        # the user specified "-q".  However, I'm not sure where
        # distutils stores that information.  Also, install has a
        # --force option that gets saved in self.user_options.  It
        # means overwrite previous installations.  If the user has
        # forced an installation, should we also ignore dependencies?
        dependencies = [
            ("mxTextTools", is_mxTextTools_installed, 1,
             "http://www.lemburg.com/files/python/mxExtensions.html"),
            ("Martel", is_Martel_installed, 1,
             "http://www.biopython.org/~dalke/Martel/"),
            ("Numerical Python", is_Numpy_installed, 0,
             "http://numpy.sourceforge.net/"),
            ("Reportlab", is_reportlab_installed, 0,
             "http://www.reportlab.com/download.html"),
            ]

        for name, is_installed_fn, is_required, url in dependencies:
            if is_installed_fn():
                continue
            
            print "*** %s *** is either not installed or out of date." % name
            if is_required:
                
                print """
This package is required for many Biopython features.  Please install
it before you install Biopython."""
                default = 0
            else:
                print """
This package is optional, which means it is only used in a few
specialized modules in Biopython.  You probably don't need this is you
are unsure.  You can ignore this requirement, and install it later if
you see ImportErrors."""
                default = 1
            print "You can find %s at %s." % (name, url)
            print

            if not get_yes_or_no(
                "Do you want to continue this installation?", default):
                return 0
        return 1

    def install_eutils_dtds(self):
        """This is a hack to install DTDs needed for EUtils into Bio.
        """
        dtds_dir = os.path.join(os.getcwd(), "Bio", "EUtils", "DTDs")
        install_dir = os.path.join(self.install_purelib, "Bio", "EUtils",
                                   "DTDs")
        potential_dtds = os.listdir(dtds_dir)
        for potential_dtd in potential_dtds:
            name, ext = os.path.splitext(potential_dtd)
            if ext == ".dtd":
                shutil.copy(os.path.join(dtds_dir, potential_dtd),
                            install_dir)
         
        
    def run(self):
        if self.check_dependencies():
            # Run the normal install.
            install.run(self)
            self.install_eutils_dtds()

class build_py_biopython(build_py):
    def run(self):
        # Check to see if Martel is installed.  If not, then install
        # it automatically.
        if not is_Martel_installed():
            self.packages.append("Martel")
        build_py.run(self)

class test_biopython(Command):
    """Run all of the tests for the package.

    This is a automatic test run class to make distutils kind of act like
    perl. With this you can do:

    python setup.py build
    python setup.py install
    python setup.py test
    
    """
    description = "Automatically run the test suite for Biopython."
    user_options = []  # distutils complains if this is not here.
    def initialize_options(self):  # distutils wants this
        pass
    def finalize_options(self):    # this too
        pass
    def run(self):
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("Tests")
        import run_tests
        run_tests.main([])

        # change back to the current directory
        os.chdir(this_dir)

def can_import(module_name):
    """can_import(module_name) -> module or None"""
    try:
        return __import__(module_name)
    except ImportError:
        return None
    raise AssertionError, "how did I get here?"

def is_Martel_installed():
    old_path = sys.path[:]

    # First, check the version of the Martel that's bundled with
    # Biopython.
    sys.path.insert(0, '')   # Make sure I'm importing the current one.
    m = can_import("Martel")
    sys.path = old_path
    if m:
        bundled_martel_version = m.__version__
    else:
        bundled_martel_version = None
    del sys.modules["Martel"]   # Delete the old version of Martel.

    # Now try and import a Martel that's not bundled with Biopython.
    # To do that, I need to delete all the references to the current
    # path from sys.path.
    i = 0
    while i < len(sys.path):
        if sys.path[i] in ['', '.', os.getcwd()]:
            del sys.path[i]
        else:
            i += 1
    m = can_import("Martel")
    sys.path = old_path
    if m:
        old_martel_version = m.__version__
    else:
        old_martel_version = None

    installed = 0
    # If the bundled one is the same or older, then ignore it
    if old_martel_version and bundled_martel_version and \
           bundled_martel_version <= old_martel_version:
        installed = 1
    return installed

def is_mxTextTools_installed():
    if can_import("TextTools"):
        return 1
    return can_import("mx.TextTools")

def is_Numpy_installed():
    return can_import("Numeric")

def is_reportlab_installed():
    return can_import("reportlab")
                  
# --- set up the packages we are going to install
# standard biopython packages
PACKAGES = [
    'Bio',
    'Bio.Ais',
    'Bio.Align',
    'Bio.Alphabet',
    'Bio.Application',
    'Bio.Blast',
    'Bio.builders',
    'Bio.builders.Search',
    'Bio.builders.SeqRecord',
    'Bio.CDD',
    'Bio.Clustalw',
    'Bio.Cluster',
    'Bio.config',
    'Bio.Crystal',
    'Bio.Data',
    'Bio.dbdefs',
    'Bio.ECell',
    'Bio.Emboss',
    'Bio.Encodings',
    'Bio.Enzyme',
    'Bio.expressions',
    'Bio.expressions.blast',
    'Bio.expressions.embl',
    'Bio.expressions.swissprot',
    'Bio.EUtils',
    'Bio.EUtils.DTDs',
    'Bio.Fasta',
    'Bio.formatdefs',
    'Bio.FSSP',
    'Bio.GA',
    'Bio.GA.Crossover',
    'Bio.GA.Mutation',
    'Bio.GA.Repair',
    'Bio.GA.Selection',
    'Bio.GenBank',
    'Bio.Geo',
    'Bio.GFF',
    'Bio.Gobase',
    'Bio.Graphics',
    'Bio.HMM',
    'Bio.IntelliGenetics',
    'Bio.InterPro',
    'Bio.Kabat',
    'Bio.KDTree',
    'Bio.KEGG',
    'Bio.KEGG.Compound',
    'Bio.KEGG.Enzyme',
    'Bio.KEGG.Map',
    'Bio.LocusLink',
    'Bio.Medline',
    'Bio.MetaTool',
    'Bio.Mindy',
    'Bio.MultiProc',
    'Bio.NBRF',
    'Bio.Ndb',
    'Bio.NeuralNetwork',
    'Bio.NeuralNetwork.BackPropagation',
    'Bio.NeuralNetwork.Gene',
    'Bio.Parsers',
    'Bio.Pathway',
    'Bio.Pathway.Rep',
    'Bio.PDB',
    'Bio.Prosite',
    'Bio.Rebase',
    'Bio.Saf',
    'Bio.SCOP',
    'Bio.SCOP.tests',
    'Bio.SeqIO',
    'Bio.SeqUtils',
    'Bio.SubsMat',
    'Bio.SVDSuperimposer',
    'Bio.SwissProt',
    'Bio.UniGene',
    'Bio.writers',
    'Bio.writers.SeqRecord',
    'Bio.WWW',
    ]

EXTENSIONS = [
    Extension('Bio.cSVM',
              ['Bio/cSVMmodule.c',
               'Bio/csupport.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.ckMeans',
              ['Bio/ckMeansmodule.c',
               'Bio/csupport.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.clistfns',
              ['Bio/clistfnsmodule.c']
              ),
    Extension('Bio.cmathfns',
              ['Bio/cmathfnsmodule.c',
               'Bio/csupport.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.cstringfns',
              ['Bio/cstringfnsmodule.c']
              ),
    Extension('Bio.cdistance',
              ['Bio/cdistancemodule.c',
               'Bio/csupport.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.cpairwise2',
              ['Bio/cpairwise2module.c',
               'Bio/csupport.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.trie',
              ['Bio/triemodule.c',
               'Bio/trie.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.cMarkovModel',
              ['Bio/cMarkovModelmodule.c',
               'Bio/csupport.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.Cluster.cluster',
              ['Bio/Cluster/clustermodule.c',
               'Bio/Cluster/cluster.c',
               'Bio/Cluster/ranlib.c',
               'Bio/Cluster/com.c',
               'Bio/Cluster/linpack.c'],
              include_dirs=["Bio/Cluster"]
              ),
    #Extension('Bio.KDTree._KDTreecmodule',
    #          ["Bio/KDTree/_KDTree.C", 
    #           "Bio/KDTree/_KDTree.swig.C"],
    #          libraries=["stdc++"]
    #          ),
    ]


# Install BioSQL.
PACKAGES.append("BioSQL")

setup(
    name='biopython',
    version='1.20',
    author='The Biopython Consortium',
    author_email='biopython@biopython.org',
    url='http://www.biopython.org/',
    cmdclass={
        "install" : install_biopython,
        "build_py" : build_py_biopython,
        "test" : test_biopython,
        },
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    )
