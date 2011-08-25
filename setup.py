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
import shutil

def is_pypy():
    import platform
    try:
        if platform.python_implementation()=='PyPy':
            return True
    except AttributeError:
        #New in Python 2.6, not in Jython yet either
        pass
    return False

def get_yes_or_no(question, default):
    if default:
        option_str = "(Y/n)"
        default_str = 'y'
    else:
        option_str = "(y/N)"
        default_str = 'n'

    while True:
        print ("%s %s:" % (question, option_str))
        if sys.version_info[0] == 3:
            response = input().lower()
        else:
            response = raw_input().lower()
        if not response:
            response = default_str
        if response[0] in ['y', 'n']:
            break
        print ("Please answer y or n.")
    return response[0] == 'y'

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 5):
    print ("Biopython requires Python 2.5 or better (but not Python 3 " \
          + "yet).  Python %d.%d detected" % sys.version_info[:2])
    sys.exit(-1)
elif sys.version_info[0] == 3:
    print("WARNING - Biopython does not yet officially support Python 3")
    import do2to3
    python3_source = "build/py%i.%i" % sys.version_info[:2]
    if "clean" in sys.argv:
        if os.path.isdir(python3_source):
            shutil.rmtree(python3_source)
        del python3_source #so we don't try to change to it below
    else:
        if not os.path.isdir("build"):
            os.mkdir("build")
        do2to3.main(".", python3_source)

# use setuptools, falling back on core modules if not found
try:
    from setuptools import setup, Command
    from setuptools.command.install import install
    from setuptools.command.build_py import build_py
    from setuptools.command.build_ext import build_ext
    from setuptools.extension import Extension
    _SETUPTOOLS = True
except ImportError:
    from distutils.core import setup
    from distutils.core import Command
    from distutils.command.install import install
    from distutils.command.build_py import build_py
    from distutils.command.build_ext import build_ext
    from distutils.extension import Extension
    _SETUPTOOLS = False

_CHECKED = None
def check_dependencies_once():
    # Call check_dependencies, but cache the result for subsequent
    # calls.
    global _CHECKED
    if _CHECKED is None:
        _CHECKED = check_dependencies()
    return _CHECKED

def get_install_requires():
    install_requires = []
    # skip this with distutils (otherwise get a warning)
    if not _SETUPTOOLS:
        return []
    # skip this with jython and pypy
    if os.name=="java" or is_pypy():
        return []
    # check for easy_install and pip
    is_automated = False
    # easy_install: --dist-dir option passed
    try:
        dist_dir_i = sys.argv.index("--dist-dir")
    except ValueError:
        dist_dir_i = None
    if dist_dir_i is not None:
        dist_dir = sys.argv[dist_dir_i+1]
        if "egg-dist-tmp" in dist_dir:
            is_automated = True
    # pip -- calls from python directly with "-c"
    if sys.argv in [["-c", "develop", "--no-deps"],
                    ["--no-deps", "-c", "develop"],
                    ["-c", "egg_info"]]:
        is_automated = True
    if is_automated:
        global _CHECKED
        if _CHECKED is None: _CHECKED = True
        install_requires.append("numpy >= 1.5.1")
    return install_requires

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

    if os.name=='java':
        return True #NumPy is not avaliable for Jython (for now)
    if is_pypy():
        return True #Full NumPy not available for PyPy (for now)

    print ("""
Numerical Python (NumPy) is not installed.

This package is required for many Biopython features.  Please install
it before you install Biopython. You can install Biopython anyway, but
anything dependent on NumPy will not work. If you do this, and later
install NumPy, you should then re-install Biopython.

You can find NumPy at http://numpy.scipy.org
""")
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
    # Adds support for the single-version-externally-managed flag
    # which is present in setuptools but not distutils. pip requires it.
    # In setuptools this forces installation the "old way" which we
    # only support here, so we just make it a no-op.
    user_options = install.user_options + [
        ('single-version-externally-managed', None,
            "used by system package builders to create 'flat' eggs"),
    ]
    boolean_options = install.boolean_options + [
        'single-version-externally-managed',
    ]
    def initialize_options(self):
        install.initialize_options(self)
        self.single_version_externally_managed = None

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
    if is_pypy():
        return False
    return bool(can_import("numpy"))

# --- set up the packages we are going to install
# standard biopython packages
PACKAGES = [
    'Bio',
    'Bio.Align',
    'Bio.Align.Applications',
    'Bio.AlignIO',
    'Bio.Alphabet',
    'Bio.Application',
    'Bio.Blast',
    'Bio.CAPS',
    'Bio.Compass',
    'Bio.Crystal',
    'Bio.Data',
    'Bio.Emboss',
    'Bio.Entrez',
    'Bio.ExPASy',
    'Bio.FSSP',
    'Bio.GA',
    'Bio.GA.Crossover',
    'Bio.GA.Mutation',
    'Bio.GA.Repair',
    'Bio.GA.Selection',
    'Bio.GenBank',
    'Bio.Geo',
    'Bio.Graphics',
    'Bio.Graphics.GenomeDiagram',
    'Bio.HMM',
    'Bio.KEGG',
    'Bio.KEGG.Compound',
    'Bio.KEGG.Enzyme',
    'Bio.KEGG.Map',
    'Bio.Medline',
    'Bio.Motif',
    'Bio.Motif.Parsers',
    'Bio.Motif.Applications',
    'Bio.NeuralNetwork',
    'Bio.NeuralNetwork.BackPropagation',
    'Bio.NeuralNetwork.Gene',
    'Bio.Nexus',
    'Bio.NMR',
    'Bio.Pathway',
    'Bio.Pathway.Rep',
    'Bio.PDB',
    'Bio.PDB.mmCIF',
    'Bio.PopGen',
    'Bio.PopGen.Async',
    'Bio.PopGen.FDist',
    'Bio.PopGen.GenePop',
    'Bio.PopGen.SimCoal',
    'Bio.Restriction',
    'Bio.Restriction._Update',
    'Bio.SCOP',
    'Bio.SeqIO',
    'Bio.SeqUtils',
    'Bio.Sequencing',
    'Bio.Sequencing.Applications',
    'Bio.Statistics',
    'Bio.SubsMat',
    'Bio.SVDSuperimposer',
    'Bio.SwissProt',
    'Bio.TogoWS',
    'Bio.Phylo',
    'Bio.Phylo.Applications',
    'Bio.Phylo.PAML',
    'Bio.UniGene',
    'Bio.Wise',
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
elif is_pypy():
    # Skip C extensions for now
    EXTENSIONS = []
elif sys.version_info[0] == 3:
    # TODO - Must update our C extensions for Python 3
    EXTENSIONS = [
    Extension('Bio.cpairwise2',
              ['Bio/cpairwise2module.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.Nexus.cnexus',
              ['Bio/Nexus/cnexus.c']
              ),
    ]
else :
    EXTENSIONS = [
    Extension('Bio.cpairwise2',
              ['Bio/cpairwise2module.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.trie',
              ['Bio/triemodule.c',
               'Bio/trie.c'],
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
    ]

#Add extensions that requires NumPy to build
if is_Numpy_installed():
    import numpy
    numpy_include_dir = numpy.get_include()
    EXTENSIONS.append(
        Extension('Bio.Cluster.cluster',
                  ['Bio/Cluster/clustermodule.c',
                   'Bio/Cluster/cluster.c'],
                  include_dirs=[numpy_include_dir],
                  ))
    EXTENSIONS.append(
        Extension('Bio.KDTree._CKDTree',
                  ["Bio/KDTree/KDTree.c",
                   "Bio/KDTree/KDTreemodule.c"],
                  include_dirs=[numpy_include_dir],
                  ))
    EXTENSIONS.append(
        Extension('Bio.Motif._pwm',
                  ["Bio/Motif/_pwm.c"],
                  include_dirs=[numpy_include_dir],
                  ))


#We now define the Biopython version number in Bio/__init__.py
#Here we can't use "import Bio" then "Bio.__version__" as that would
#tell us the version of Biopython already installed (if any).
__version__ = "Undefined"
for line in open('Bio/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

#Simple trick to use the 2to3 converted source under Python 3,
#change the current directory before/after running setup.
#Note as a side effect there will be a build folder underneath
#the python3_source folder.
old_path = os.getcwd()
try:
    src_path = python3_source
except NameError:
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(src_path)
sys.path.insert(0, src_path)

setup_args = {
    "name" : 'biopython',
    "version" : __version__,
    "author" : 'The Biopython Consortium',
    "author_email" : 'biopython@biopython.org',
    "url" : 'http://www.biopython.org/',
    "description" : 'Freely available tools for computational molecular biology.',
    "download_url" : 'http://biopython.org/DIST/',
    "cmdclass" : {
        "install" : install_biopython,
        "build_py" : build_py_biopython,
        "build_ext" : build_ext_biopython,
        "test" : test_biopython,
        },
    "packages" : PACKAGES,
    "ext_modules" : EXTENSIONS,
    "package_data" : {
        'Bio.Entrez': ['DTDs/*.dtd', 'DTDs/*.ent', 'DTDs/*.mod'],
        'Bio.PopGen': ['SimCoal/data/*.par'],
         },
   }

if _SETUPTOOLS:
    setup_args["install_requires"] = get_install_requires()

try:
    setup(**setup_args)
finally:
    del sys.path[0]
    os.chdir(old_path)
