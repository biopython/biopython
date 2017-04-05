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
from __future__ import print_function

import sys
import os
import shutil

if "bdist_wheel" in sys.argv:
    try:
        import setuptools
        import wheel
    except ImportError:
        sys.exit("We need both setuptools AND wheel packages installed for bdist_wheel to work")
    # Import specific bits of setuptools ...
    from setuptools import setup
    from setuptools import Command
    from setuptools.command.install import install
    from setuptools.command.build_py import build_py
    from setuptools.command.build_ext import build_ext
    from setuptools import Extension
else:
    # Except for wheels, stick with standard library's distutils
    from distutils.core import setup
    from distutils.core import Command
    from distutils.command.install import install
    from distutils.command.build_py import build_py
    from distutils.command.build_ext import build_ext
    from distutils.extension import Extension


_CHECKED = None


def osx_clang_fix():
    """Add clang switch to ignore unused arguments to avoid OS X compile error.

    This is a hack to cope with Apple shipping a version of Python compiled
    with the -mno-fused-madd argument which clang from XCode 5.1 does not
    support::

        $ cc -v
        Apple LLVM version 5.1 (clang-503.0.40) (based on LLVM 3.4svn)
        Target: x86_64-apple-darwin13.2.0
        Thread model: posix

        $ which python-config
        /Library/Frameworks/Python.framework/Versions/Current/bin/python-config

        $ python-config --cflags
        -I/Library/Frameworks/Python.framework/Versions/2.5/include/python2.5
        -I/Library/Frameworks/Python.framework/Versions/2.5/include/python2.5
        -arch ppc -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk
        -fno-strict-aliasing -Wno-long-double -no-cpp-precomp -mno-fused-madd
        -fno-common -dynamic -DNDEBUG -g -O3

    We can avoid the clang compilation error with -Qunused-arguments which is
    (currently) harmless if gcc is being used instead (e.g. compiling Biopython
    against a locally compiled Python rather than the Apple provided Python).
    """
    # see http://lists.open-bio.org/pipermail/biopython-dev/2014-April/011240.html
    if sys.platform != "darwin":
        return
    # see also Bio/_py3k/__init__.py (which we can't use in setup.py)
    if sys.version_info[0] >= 3:
        from subprocess import getoutput
    else:
        from commands import getoutput
    cc = getoutput("cc -v")
    if "gcc" in cc or "clang" not in cc:
        return
    for flag in ["CFLAGS", "CPPFLAGS"]:
        if flag not in os.environ:
            os.environ[flag] = "-Qunused-arguments"
        elif "-Qunused-arguments" not in os.environ[flag]:
            os.environ[flag] += " -Qunused-arguments"

osx_clang_fix()


def is_pypy():
    import platform
    try:
        if platform.python_implementation() == 'PyPy':
            return True
    except AttributeError:
        # New in Python 2.6, not in Jython yet either
        pass
    return False


def is_ironpython():
    return sys.platform == "cli"
    # TODO - Use platform as in Pypy test?


def get_yes_or_no(question, default):
    if default:
        option_str = "(Y/n)"
        default_str = 'y'
    else:
        option_str = "(y/N)"
        default_str = 'n'

    while True:
        print("%s %s:" % (question, option_str))
        if sys.version_info[0] == 3:
            response = input().lower()
        else:
            response = raw_input().lower()
        if not response:
            response = default_str
        if response[0] in ['y', 'n']:
            break
        print("Please answer y or n.")
    return response[0] == 'y'


# Make sure we have the right Python version.
if sys.version_info[:2] < (2, 7):
    sys.stderr.write("Biopython requires Python 2.7, or Python 3.3 or later. "
                     "Python %d.%d detected.\n" % sys.version_info[:2])
    sys.exit(1)
elif sys.version_info[0] == 3 and sys.version_info[:2] < (3, 3):
    sys.stderr.write("Biopython requires Python 3.3 or later (or Python 2.7). "
                     "Python %d.%d detected.\n" % sys.version_info[:2])
    sys.exit(1)
elif sys.version_info[:2] == (3, 3):
    sys.stderr.write("WARNING: Biopython support for Python 3.3 is now deprecated.\n")


def check_dependencies_once():
    # Call check_dependencies, but cache the result for subsequent
    # calls.
    global _CHECKED
    if _CHECKED is None:
        _CHECKED = check_dependencies()
    return _CHECKED


def is_automated():
    """Check for installation with easy_install or pip.
    """
    is_automated = False
    # easy_install: --dist-dir option passed
    try:
        dist_dir_i = sys.argv.index("--dist-dir")
    except ValueError:
        dist_dir_i = None
    if dist_dir_i is not None:
        dist_dir = sys.argv[dist_dir_i + 1]
        if "egg-dist-tmp" in dist_dir:
            is_automated = True
    # pip -- calls from python directly with "-c"
    if sys.argv in [["-c", "develop", "--no-deps"],
                    ["--no-deps", "-c", "develop"],
                    ["-c", "egg_info"]] \
                    or "pip-egg-info" in sys.argv \
                    or sys.argv[:3] == ["-c", "install", "--record"] \
                    or sys.argv[:4] == ['-c', 'install', '--single-version-externally-managed',
                                        '--record']:
        is_automated = True
    return is_automated


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
    if is_Numpy_installed():
        return True
    if is_automated():
        return True  # For automated builds go ahead with installed packages
    if os.name == 'java':
        return True  # NumPy is not avaliable for Jython (for now)
    if is_ironpython():
        return True  # We're ignoring NumPy under IronPython (for now)

    print("""
Numerical Python (NumPy) is not installed.

This package is required for many Biopython features.  Please install
it before you install Biopython. You can install Biopython anyway, but
anything dependent on NumPy will not work. If you do this, and later
install NumPy, you should then re-install Biopython.

You can find NumPy at http://www.numpy.org
""")
    # exit automatically if running as part of some script
    # (e.g. PyPM, ActiveState's Python Package Manager)
    if not sys.stdout.isatty():
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
        if os.name == "java" and "Bio.Restriction" in self.packages:
            # Evil hack to work on Jython 2.7
            # This is to avoid java.lang.RuntimeException: Method code too large!
            # from Bio/Restriction/Restriction_Dictionary.py
            self.packages.remove("Bio.Restriction")
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
    'Bio.codonalign',
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
    'Bio.PDB.mmtf',
    'Bio.KEGG.KGML',
    'Bio.Medline',
    'Bio.motifs',
    'Bio.motifs.applications',
    'Bio.motifs.jaspar',
    'Bio.NeuralNetwork',
    'Bio.NeuralNetwork.BackPropagation',
    'Bio.NeuralNetwork.Gene',
    'Bio.Nexus',
    'Bio.NMR',
    'Bio.Pathway',
    'Bio.Pathway.Rep',
    'Bio.PDB',
    'Bio.PopGen',
    'Bio.PopGen.Async',
    'Bio.PopGen.FDist',
    'Bio.PopGen.GenePop',
    'Bio.PopGen.SimCoal',
    'Bio.Restriction',
    'Bio.SCOP',
    'Bio.SearchIO',
    'Bio.SearchIO._model',
    'Bio.SearchIO.BlastIO',
    'Bio.SearchIO.HmmerIO',
    'Bio.SearchIO.ExonerateIO',
    'Bio.SeqIO',
    'Bio.SeqUtils',
    'Bio.Sequencing',
    'Bio.Sequencing.Applications',
    'Bio.Statistics',
    'Bio.SubsMat',
    'Bio.SVDSuperimposer',
    'Bio.PDB.QCPSuperimposer',
    'Bio.SwissProt',
    'Bio.TogoWS',
    'Bio.Phylo',
    'Bio.Phylo.Applications',
    'Bio.Phylo.PAML',
    'Bio.UniGene',
    'Bio.UniProt',
    'Bio.Wise',
    'Bio._py3k',
    # Other top level packages,
    'BioSQL',
    ]

if os.name == 'jython':
    # Evil hack to work on Jython 2.7
    # This is to avoid java.lang.RuntimeException: Method code too large!
    # from Bio/Restriction/Restriction_Dictionary.py
    PACKAGES.remove('Bio.Restriction')


# packages that require Numeric Python
NUMPY_PACKAGES = [
    'Bio.Affy',
    'Bio.Cluster',
    'Bio.KDTree',
    'Bio.phenotype',
]

if os.name == 'java':
    # Jython doesn't support C extensions
    EXTENSIONS = []
elif is_ironpython():
    # Skip C extensions for now
    EXTENSIONS = []
elif is_pypy():
    # Two out of three ain't bad?
    EXTENSIONS = [
    Extension('Bio.cpairwise2',
              ['Bio/cpairwise2module.c'],
              ),
    # Bio.trie has a problem under PyPy2 v5.6 and 5.7
    Extension('Bio.Nexus.cnexus',
              ['Bio/Nexus/cnexus.c']
              ),
    ]
else:
    EXTENSIONS = [
    Extension('Bio.cpairwise2',
              ['Bio/cpairwise2module.c'],
              ),
    Extension('Bio.trie',
              ['Bio/triemodule.c',
               'Bio/trie.c'],
              include_dirs=["Bio"]
              ),
    Extension('Bio.Nexus.cnexus',
              ['Bio/Nexus/cnexus.c']
              ),
    ]

# Add extensions that requires NumPy to build
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
        Extension('Bio.motifs._pwm',
                  ["Bio/motifs/_pwm.c"],
                  include_dirs=[numpy_include_dir],
                  ))
    EXTENSIONS.append(
        Extension('Bio.PDB.QCPSuperimposer.qcprotmodule',
                  ["Bio/PDB/QCPSuperimposer/qcprotmodule.c"],
                  include_dirs=[numpy_include_dir],
                  ))


# We now define the Biopython version number in Bio/__init__.py
# Here we can't use "import Bio" then "Bio.__version__" as that would
# tell us the version of Biopython already installed (if any).
__version__ = "Undefined"
for line in open('Bio/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

# Simple trick to use the 2to3 converted source under Python 3,
# change the current directory before/after running setup.
# Note as a side effect there will be a build folder underneath
# the python3_source folder.
old_path = os.getcwd()
try:
    src_path = python3_source
except NameError:
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
os.chdir(src_path)
sys.path.insert(0, src_path)

setup_args = {
    "name": 'biopython',
    "version": __version__,
    "author": 'The Biopython Contributors',
    "author_email": 'biopython@biopython.org',
    "url": 'http://www.biopython.org/',
    "description": 'Freely available tools for computational molecular biology.',
    "download_url": 'http://biopython.org/DIST/',
    "cmdclass": {
        "install": install_biopython,
        "build_py": build_py_biopython,
        "build_ext": build_ext_biopython,
        "test": test_biopython,
        },
    "packages": PACKAGES,
    "ext_modules": EXTENSIONS,
    "package_data": {
        'Bio.Entrez': ['DTDs/*.dtd', 'DTDs/*.ent', 'DTDs/*.mod'],
        'Bio.PopGen': ['SimCoal/data/*.par'],
         },
   }

try:
    setup(**setup_args)
finally:
    del sys.path[0]
    os.chdir(old_path)
