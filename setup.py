"""Distutils based setup script for Biopython.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install

For more in-depth instructions, see the installation section of the
Biopython manual, linked to from:

http://biopython.org/wiki/Documentation

Or for more details about the options available from distutils, look at
the 'Installing Python Modules' distutils documentation, available from:

http://python.org/sigs/distutils-sig/doc/

Or, if all else fails, feel free to write to the sign up to the Biopython
mailing list and ask for help.  See:

http://biopython.org/wiki/Mailing_lists
"""
import sys
import os

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 4):
    print "Biopython requires Python 2.4 or better.  Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)

from distutils.core import setup
from distutils.core import Command
from distutils.command.install import install
from distutils.command.build_py import build_py
from distutils.command.build_ext import build_ext
from distutils.extension import Extension

def get_yes_or_no(question, default):
    if default:
        option_str = "(Y/n)"
        default_str = 'y'
    else:
        option_str = "(y/N)"
        default_str = 'n'

    while True:
        print "%s %s " % (question, option_str),
        response = raw_input().lower()
        if not response:
            response = default_str
        if response[0] in ['y', 'n']:
            break
        print "Please answer y or n."
    return response[0] == 'y'

_CHECKED = None
def check_dependencies_once():
    # Call check_dependencies, but cache the result for subsequent
    # calls.
    global _CHECKED
    if _CHECKED is None:
        _CHECKED = check_dependencies()
    return _CHECKED

def check_dependencies():
    """Return whether the installation should continue."""
    # There should be some way for the user to tell specify not to
    # check dependencies.  For example, it probably should not if
    # the user specified "-q".  However, I'm not sure where
    # distutils stores that information.  Also, install has a
    # --force option that gets saved in self.user_options.  It
    # means overwrite previous installations.  If the user has
    # forced an installation, should we also ignore dependencies?

    # We only check for NumPy, as this is a compile time dependency
    if is_Numpy_installed() : return True
    print """
Numerical Python (NumPy) is not installed.

This package is required for many Biopython features.  Please install
it before you install Biopython. You can install Biopython anyway, but
anything dependent on NumPy will not work. If you do this, and later
install NumPy, you should then re-install Biopython.

You can find NumPy at http://numpy.scipy.org
"""
    # exit automatically if running as part of some script
    # (e.g. PyPM, ActiveState's Python Package Manager)
    if not sys.stdout.isatty() :
        sys.exit(-1)
    # We can ask the user
    return get_yes_or_no("Do you want to continue this installation?", False)

class install_biopython(install):
    """Override the standard install to check for dependencies.

    This will just run the normal install, and then print warning messages
    if packages are missing.

    """
    def run(self):
        if check_dependencies_once():
            # Run the normal install.
            install.run(self)

class build_py_biopython(build_py):
    def run(self):
        if not check_dependencies_once():
            return
        # Add software that requires Numpy to be installed.
        if is_Numpy_installed():
            self.packages.extend(NUMPY_PACKAGES)
        build_py.run(self)


class build_ext_biopython(build_ext):
    def run(self):
        if not check_dependencies_once():
            return
        # add software that requires NumPy to install
        if is_Numpy_installed():
            import numpy
            numpy_include_dir = numpy.get_include()
            self.extensions.append(
                Extension('Bio.Cluster.cluster',
                          ['Bio/Cluster/clustermodule.c',
                           'Bio/Cluster/cluster.c'],
                          include_dirs=[numpy_include_dir],
                          ))
            self.extensions.append(
                Extension('Bio.KDTree._CKDTree',
                          ["Bio/KDTree/KDTree.c",
                           "Bio/KDTree/KDTreemodule.c"],
                          include_dirs=[numpy_include_dir],
                          ))
            self.extensions.append(
                Extension('Bio.Motif._pwm',
                          ["Bio/Motif/_pwm.c"],
                          include_dirs=[numpy_include_dir],
                          ))
        build_ext.run(self)


class test_biopython(Command):
    """Run all of the tests for the package.

    This is a automatic test run class to make distutils kind of act like
    perl. With this you can do:

    python setup.py build
    python setup.py install
    python setup.py test
    
    """
    description = "Automatically run the test suite for Biopython."
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("Tests")
        sys.path.insert(0, '')
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

def is_Numpy_installed():
    return bool(can_import("numpy"))

# --- set up the packages we are going to install
# standard biopython packages
PACKAGES = [
    'Bio',
    'Bio.Align',
    'Bio.Align.Applications',
    'Bio.AlignIO',
    'Bio.AlignAce',
    'Bio.Alphabet',
    'Bio.Application',
    'Bio.Blast',
    'Bio.CAPS',
    'Bio.Compass',
    'Bio.Clustalw',
    'Bio.Crystal',
    'Bio.Data',
    'Bio.Emboss',
    'Bio.Encodings',
    'Bio.Entrez',
    'Bio.Enzyme',
    'Bio.ExPASy',
    'Bio.Fasta',
    'Bio.FSSP',
    'Bio.GA',
    'Bio.GA.Crossover',
    'Bio.GA.Mutation',
    'Bio.GA.Repair',
    'Bio.GA.Selection',
    'Bio.GenBank',
    'Bio.Geo',
    'Bio.GFF',
    'Bio.Graphics',
    'Bio.Graphics.GenomeDiagram',
    'Bio.HMM',
    'Bio.InterPro',
    'Bio.KEGG',
    'Bio.KEGG.Compound',
    'Bio.KEGG.Enzyme',
    'Bio.KEGG.Map',
    'Bio.Medline',
    'Bio.MEME',
    'Bio.Motif',
    'Bio.Motif.Parsers',
    'Bio.Motif.Applications',
    'Bio.Ndb',
    'Bio.NeuralNetwork',
    'Bio.NeuralNetwork.BackPropagation',
    'Bio.NeuralNetwork.Gene',
    'Bio.Nexus',
    'Bio.NMR',
    'Bio.Parsers',
    'Bio.Pathway',
    'Bio.Pathway.Rep',
    'Bio.PDB',
    'Bio.PDB.mmCIF',
    'Bio.PopGen',
    'Bio.PopGen.Async',
    'Bio.PopGen.FDist',
    'Bio.PopGen.GenePop',
    'Bio.PopGen.SimCoal',
    'Bio.Prosite',
    'Bio.Restriction',
    'Bio.Restriction._Update',
    'Bio.SCOP',
    'Bio.SeqIO',
    'Bio.SeqUtils',
    'Bio.Sequencing',
    'Bio.Statistics',
    'Bio.SubsMat',
    'Bio.SVDSuperimposer',
    'Bio.SwissProt',
    'Bio.UniGene',
    'Bio.Wise',
    'Bio.WWW',
    #Other top level packages,
    'BioSQL',
    ]

# packages that require Numeric Python
NUMPY_PACKAGES = [
    'Bio.Affy',
    'Bio.Cluster',
    'Bio.KDTree',
]

if os.name == 'java' :
    # Jython doesn't support C extensions
    EXTENSIONS = []
else :
    EXTENSIONS = [
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
#Commented out due to the build dependency on flex, see Bug 2619
#   Extension('Bio.PDB.mmCIF.MMCIFlex',
#              ['Bio/PDB/mmCIF/lex.yy.c',
#               'Bio/PDB/mmCIF/MMCIFlexmodule.c'],
#              include_dirs=["Bio"],
#              libraries=["fl"]
#              ),
    Extension('Bio.Nexus.cnexus',
              ['Bio/Nexus/cnexus.c']
              ),
    Extension('Bio.Restriction.DNAUtils',
              ['Bio/Restriction/DNAUtils.c']
              ),
    ]


#We now define the Biopython version number in Bio/__init__.py
#Here we can't use "import Bio" then "Bio.__version__" as that would
#tell us the version of Biopython already installed (if any).
__version__ = "Undefined"
for line in open('Bio/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

setup(
    name='biopython',
    version=__version__,
    author='The Biopython Consortium',
    author_email='biopython@biopython.org',
    url='http://www.biopython.org/',
    description='Freely available tools for computational molecular biology.',
    download_url='http://biopython.org/DIST/',
    cmdclass={
        "install" : install_biopython,
        "build_py" : build_py_biopython,
        "build_ext" : build_ext_biopython,
        "test" : test_biopython,
        },
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    package_data = {'Bio.Entrez': ['DTDs/*.dtd'],
                    'Bio.PopGen': ['SimCoal/data/*.par'],
                   },
    #install_requires = ['numpy>=1.0'],
    #extras_require = {
    #    'PDF' : ['reportlab>=2.0']
    #    }
    )
