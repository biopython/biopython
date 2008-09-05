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

Or, if all else fails, feel free to write to the sign up to the Biopython
mailing list and ask for help.  See http://biopython.org/wiki/Mailing_lists
"""
import sys
import os

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 3):
    print "Biopython requires Python 2.3 or better.  Python %d.%d detected" % \
          sys.version_info[:2]
    sys.exit(-1)

from distutils.core import setup
from distutils.core import Command
from distutils.command.install import install
from distutils.command.install_data import install_data
from distutils.command.build_py import build_py
from distutils.command.build_ext import build_ext
from distutils.extension import Extension
from distutils import sysconfig

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
    
    #This is a list of tuples, containing:
    # - package name, string
    # - is packaged installed, boolean
    # - is packaged required, boolean
    # - package website, string
    dependencies = [
        ("mxTextTools", is_mxTextTools_installed, 0,
         "http://www.egenix.com/files/python/eGenix-mx-Extensions.html"),
        ("Numerical Python", is_Numpy_installed, 0,
         "http://numpy.sourceforge.net/"),
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
specialized modules in Biopython.  You probably don't need this if you
are unsure.  You can ignore this requirement, and install it later if
you see ImportErrors."""
            default = 1
        print "You can find %s at %s." % (name, url)
        print
        # exit automatically if required packages not installed
        if not(default):
            sys.exit(-1)

        if not get_yes_or_no(
            "Do you want to continue this installation?", default):
            return 0
        
    
    # Compile KDTree ? Not compiled by default
    print "\n*** Bio.KDTree *** NOT built by default "
    kdtree_msg = """
The Bio.PDB.NeighborSearch module depends on the Bio.KDTree module,
which in turn, depends on C++ code that does not compile cleanly on
all platforms. Hence, Bio.KDTree is not built by default.

Would you like to build Bio.KDTree ?"""

    if get_yes_or_no (kdtree_msg, 0):
        NUMPY_PACKAGES.append("Bio.KDTree")
        NUMPY_EXTENSIONS.append(
            CplusplusExtension('Bio.KDTree._CKDTree', 
                               ["Bio/KDTree/KDTree.cpp",
                                "Bio/KDTree/KDTree.swig.cpp"],
                               libraries=["stdc++"],
                               language="c++"))
    
    
    return 1

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
        # Check to see if Martel is installed.  If not, then install
        # it automatically.
        if not is_Martel_installed():
            self.packages.append("Martel")
        # Add software that requires Numpy to be installed.
        if is_Numpy_installed():
            self.packages.extend(NUMPY_PACKAGES)
        build_py.run(self)

        # In addition to installing the data files, we also need to make
        # sure that they are copied to the build directory. Otherwise,
        # the unit tests will fail because they cannot find the data files
        # in the build directory.
        # This is taken care of automatically in Python 2.4 or higher by
        # using package_data.

        import glob
        data_files = self.distribution.data_files
        for entry in data_files:
            if type(entry) is not type(""):
                raise ValueError, "data_files must be strings"
            # Unix- to platform-convention conversion
            entry = os.sep.join(entry.split("/"))
            filenames = glob.glob(entry)
            for filename in filenames:
                dst = os.path.join(self.build_lib, filename)
                dstdir = os.path.split(dst)[0]
                self.mkpath(dstdir)
                self.copy_file(filename, dst)


class CplusplusExtension(Extension):
    """Hack-y wrapper around Extension to support C++ and Python2.2.

    Python2.3 defines an extension attribute, which can be used in
    'build_extension' to work around problems Python has with always
    using the C++ compiler to compile C++ code.
    
    This should be able to be removed once we move to requiring Python 2.3 or
    better.
    """
    def __init__(self, *args, **kw):
        # fix the language -- 2.2 doesn't have languages
        if sys.version_info[1] < 3:
            try:
                self.language = kw['language']
                del kw['language']
            except KeyError:
                pass
        Extension.__init__(self, *args, **kw)

class build_ext_biopython(build_ext):
    def run(self):
        if not check_dependencies_once():
            return
        # add software that requires NumPy to install
        if is_Numpy_installed():
            self.extensions.extend(NUMPY_EXTENSIONS)
        build_ext.run(self)

    def build_extensions(self):
        # Unix C compiler plus others
        if hasattr(self.compiler, "compiler_so"):
            self._original_compiler_so = self.compiler.compiler_so
        # MSVC -- others?
        else:
            self._original_compiler_so = self.compiler.cc

        build_ext.build_extensions(self)

    def build_extension(self, ext):
        """Work around distutils bug which uses the C compiler for C++ code.
        """
        # build this extension by default
        build = 1
        if hasattr(ext, "language") and ext.language == "c++":
            # C++ didn't build in the past with msvc
            # mingw32 seems fine now
            if self.compiler.compiler_type=="msvc":
                build = 0
            # fix for distutils where C++ is not handled well. This includes
            # Python 2.2.x -- need to find the C++ compiler
            cxx = None
            if (sys.version_info[1] < 3) and build: # Python 2.2
                cxx = sysconfig.get_config_vars("CXX")
                if os.environ.has_key("CXX"):
                    cxx = os.environ["CXX"]
            # set the C++ compiler if it doesn't exist in distutils
            if cxx:
                self.compiler.set_executable("compiler", cxx)
                self.compiler.set_executable("compiler_so", cxx)
                self.compiler.set_executable("linker_so",
                        cxx + ["-shared"])
        else:
            self.compiler.compiler_so = self._original_compiler_so

        # C++ extensions just plain won't build on some platforms
        if build:
            build_ext.build_extension(self, ext)

class test_biopython(Command):
    """Run all of the tests for the package.

    This is a automatic test run class to make distutils kind of act like
    perl. With this you can do:

    python setup.py build
    python setup.py install
    python setup.py test
    
    """
    description = "Automatically run the test suite for Biopython."

    user_options = [
        # provide the option to run tests in no-gui mode
        ('no-gui', None, "Do not run in GUI mode")
    ]

    def initialize_options(self):
        self.no_gui = None

    def finalize_options(self):
        pass

    def run(self):
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("Tests")
        sys.path.insert(0, '')
        import run_tests
        if self.no_gui:
            run_tests.main(['--no-gui'])
        else:
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
        del sys.modules["Martel"]   # Unload the bundled version of Martel.
    else:
        #We won't be able to import the bundled version of Martel if an
        #external dependency like mxTextTools is missing.
        #In this case, we can't compare versions to any pre-installed Martel
        #(even if that could be imported).
        bundled_martel_version = None

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
        #Either there is no pre-installed copy of Martel, or if there is
        #it cannot be imported (e.g. missing a dependency like mxTextTools).
        old_martel_version = None

    installed = 0
    # If the bundled one is the older, then ignore it
    if old_martel_version and bundled_martel_version and \
           bundled_martel_version < old_martel_version:
        installed = 1
    return installed

def is_mxTextTools_installed():
    if can_import("TextTools"):
        return 1
    return can_import("mx.TextTools")

def is_Numpy_installed():
    return can_import("Numeric")

# --- set up the packages we are going to install
# standard biopython packages
PACKAGES = [
    'Bio',
    'Bio.Ais',
    'Bio.Align',
    'Bio.AlignIO',
    'Bio.AlignAce',
    'Bio.Alphabet',
    'Bio.Application',
    'Bio.Blast',
    'Bio.builders',
    'Bio.builders.Search',
    'Bio.builders.SeqRecord',
    'Bio.CAPS',
    'Bio.CDD',
    'Bio.Compass',
    'Bio.Clustalw',
    'Bio.config',
    'Bio.Crystal',
    'Bio.Data',
    'Bio.dbdefs',
    'Bio.ECell',
    'Bio.Emboss',
    'Bio.Encodings',
    'Bio.Entrez',
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
    'Bio.KEGG',
    'Bio.KEGG.Compound',
    'Bio.KEGG.Enzyme',
    'Bio.KEGG.Map',
    'Bio.LocusLink',
    'Bio.Medline',
    'Bio.MEME',
    'Bio.MetaTool',
    'Bio.Mindy',
    'Bio.NBRF',
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
    'Bio.Rebase',
    'Bio.Restriction',
    'Bio.Restriction._Update',
    'Bio.Saf',
    'Bio.SCOP',
    'Bio.SeqIO',
    'Bio.SeqUtils',
    'Bio.Sequencing',
    'Bio.Statistics',
    'Bio.SubsMat',
    'Bio.SVDSuperimposer',
    'Bio.SwissProt',
    'Bio.UniGene',
    'Bio.writers',
    'Bio.writers.SeqRecord',
    'Bio.Wise',
    'Bio.WWW',
    ]

# packages that require Numeric Python
NUMPY_PACKAGES = [
    'Bio.Affy',
    'Bio.Cluster',
]

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

# extensions that require numeric python
NUMPY_EXTENSIONS = [
    Extension('Bio.Cluster.cluster',
              ['Bio/Cluster/clustermodule.c',
               'Bio/Cluster/cluster.c'],
              include_dirs=["Bio/Cluster"]
              ),
#   CplusplusExtension('Bio.Affy._cel',  # The file parser in celmodule.cc was
#            ['Bio/Affy/celmodule.cc'],  # replaced by a scanner/consumer in
#            language="c++"              # CelFile.py, using Biopython's
#            ),                          # parser framework
]

DATA_FILES=[
    "Bio/Entrez/DTDs/*.dtd",
    "Bio/EUtils/DTDs/*.dtd",
    "Bio/PopGen/SimCoal/data/*.par"
    ]

# EUtils contains dtd files that need to be installed in the same
# directory as the python modules.  Distutils doesn't have a simple
# way of handling this, and we need to subclass install_data.  This
# code is adapted from the mx.TextTools distribution.

# We can use package_data instead once we require Python 2.4 or higher.

class install_data_biopython(install_data):
    def finalize_options(self):
        if self.install_dir is None:
            installobj = self.distribution.get_command_obj('install')
            self.install_dir = installobj.install_platlib
        install_data.finalize_options(self)

    def run (self):
        import glob
        if not self.dry_run:
            self.mkpath(self.install_dir)
        data_files = self.get_inputs()
        for entry in data_files:
            if type(entry) is not type(""):
                raise ValueError, "data_files must be strings"
            # Unix- to platform-convention conversion
            entry = os.sep.join(entry.split("/"))
            filenames = glob.glob(entry)
            for filename in filenames:
                dst = os.path.join(self.install_dir, filename)
                dstdir = os.path.split(dst)[0]
                if not self.dry_run:
                    self.mkpath(dstdir)
                    outfile = self.copy_file(filename, dst)[0]
                else:
                    outfile = dst
                self.outfiles.append(outfile)


# Install BioSQL.
PACKAGES.append("BioSQL")

setup(
    name='biopython',
    version='1.48',
    author='The Biopython Consortium',
    author_email='biopython@biopython.org',
    url='http://www.biopython.org/',
    cmdclass={
        "install" : install_biopython,
        "build_py" : build_py_biopython,
        "build_ext" : build_ext_biopython,
        "install_data" : install_data_biopython,
        "test" : test_biopython,
        },
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    data_files=DATA_FILES,
    # package_data = {'Bio.Entrez': ['DTDs/*.dtd']}
    ## Use this once we require Python version >= 2.4.
    )
