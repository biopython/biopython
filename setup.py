"""setuptools based setup script for Biopython.

This uses setuptools which is now the standard python mechanism for
installing packages. If you have downloaded and uncompressed the
Biopython source code, or fetched it from git, for the simplest
installation just type the command::

    python setup.py install

However, you would normally install the latest Biopython release from
the PyPI archive with::

    pip install biopython

For more in-depth instructions, see the installation section of the
Biopython manual, linked to from:

http://biopython.org/wiki/Documentation

Or, if all else fails, feel free to write to the sign up to the Biopython
mailing list and ask for help.  See:

http://biopython.org/wiki/Mailing_lists
"""
from __future__ import print_function

import sys
import os

try:
    from setuptools import setup
    from setuptools import Command
    from setuptools.command.install import install
    from setuptools.command.build_py import build_py
    from setuptools.command.build_ext import build_ext
    from setuptools import Extension
except ImportError:
    sys.exit("We need the Python library setuptools to be installed. "
             "Try runnning: python -m ensurepip")

if "bdist_wheel" in sys.argv:
    try:
        import wheel  # noqa: F401
    except ImportError:
        sys.exit("We need both setuptools AND wheel packages installed "
                 "for bdist_wheel to work. Try running: pip install wheel")

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
    from distutils.ccompiler import new_compiler
    from distutils.sysconfig import customize_compiler
    # The compiler test should be made on the actual compiler that'll be used
    compiler = new_compiler()
    customize_compiler(compiler)
    cc = getoutput("{} -v".format(compiler.compiler[0]))
    if "gcc" in cc or "clang" not in cc:
        return
    for flag in ["CFLAGS", "CPPFLAGS"]:
        if flag not in os.environ:
            os.environ[flag] = "-Qunused-arguments"
        elif "-Qunused-arguments" not in os.environ[flag]:
            os.environ[flag] += " -Qunused-arguments"


osx_clang_fix()


def is_pypy():
    """Check if running under the PyPy implementation of Python."""
    import platform
    try:
        if platform.python_implementation() == 'PyPy':
            return True
    except AttributeError:
        # New in Python 2.6
        pass
    return False


def is_jython():
    """Check if running under the Jython implementation of Python."""
    import platform
    try:
        if platform.python_implementation() == 'Jython':
            return True
    except AttributeError:
        # This was missing prior to ~ Jython 2.7.0
        pass
    # Fall back which will work with older Jython:
    return os.name == "java"


def is_ironpython():
    """Check if running under the IronPython implementation of Python."""
    return sys.platform == "cli"
    # TODO - Use platform as in Pypy test?


# Make sure we have the right Python version.
if sys.version_info[:2] < (2, 7):
    sys.stderr.write("Biopython requires Python 2.7, or Python 3.4 or later. "
                     "Python %d.%d detected.\n" % sys.version_info[:2])
    sys.exit(1)
elif sys.version_info[0] == 3 and sys.version_info[:2] < (3, 4):
    sys.stderr.write("Biopython requires Python 3.4 or later (or Python 2.7). "
                     "Python %d.%d detected.\n" % sys.version_info[:2])
    sys.exit(1)

if is_jython():
    sys.stderr.write("WARNING: Biopython support for Jython is now deprecated.\n")


def check_dependencies_once():
    """Check dependencies, will cache and re-use the result."""
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

    # Currently there are no compile time dependencies
    return True


class install_biopython(install):
    """Override the standard install to check for dependencies.

    This will just run the normal install, and then print warning messages
    if packages are missing.
    """

    def run(self):
        """Run the installation."""
        if check_dependencies_once():
            # Run the normal install.
            install.run(self)


class build_py_biopython(build_py):
    """Biopython builder."""

    def run(self):
        """Run the build."""
        if not check_dependencies_once():
            return
        if is_jython() and "Bio.Restriction" in self.packages:
            # Evil hack to work on Jython 2.7 to avoid
            # java.lang.RuntimeException: Method code too large!
            # from Bio/Restriction/Restriction_Dictionary.py
            self.packages.remove("Bio.Restriction")
        # Add software that requires Numpy to be installed.
        if is_jython() or is_ironpython():
            pass
        else:
            self.packages.extend(NUMPY_PACKAGES)
        build_py.run(self)


class build_ext_biopython(build_ext):
    """Biopython extension builder."""

    def run(self):
        """Run the build."""
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
        """No-op, initialise options."""
        pass

    def finalize_options(self):
        """No-op, finalise options."""
        pass

    def run(self):
        """Run the tests."""
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("Tests")
        sys.path.insert(0, '')
        import run_tests
        run_tests.main([])

        # change back to the current directory
        os.chdir(this_dir)


def can_import(module_name):
    """Check we can import the requested module."""
    try:
        return __import__(module_name)
    except ImportError:
        return None


# Using requirements.txt is preferred for an application
# (and likely will pin specific version numbers), using
# setup.py's install_requires is preferred for a library
# (and should try not to be overly narrow with versions).
REQUIRES = [
    'numpy',
]

if is_jython() or is_ironpython():
    REQUIRES.remove("numpy")


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
    'Bio.GenBank',
    'Bio.Geo',
    'Bio.Graphics',
    'Bio.Graphics.GenomeDiagram',
    'Bio.HMM',
    'Bio.KEGG',
    'Bio.KEGG.Compound',
    'Bio.KEGG.Enzyme',
    'Bio.KEGG.Gene',
    'Bio.KEGG.Map',
    'Bio.PDB.mmtf',
    'Bio.KEGG.KGML',
    'Bio.Medline',
    'Bio.motifs',
    'Bio.motifs.applications',
    'Bio.motifs.jaspar',
    'Bio.Nexus',
    'Bio.NMR',
    'Bio.Pathway',
    'Bio.Pathway.Rep',
    'Bio.PDB',
    'Bio.PopGen',
    'Bio.PopGen.GenePop',
    'Bio.Restriction',
    'Bio.SCOP',
    'Bio.SearchIO',
    'Bio.SearchIO._legacy',
    'Bio.SearchIO._model',
    'Bio.SearchIO.BlastIO',
    'Bio.SearchIO.HmmerIO',
    'Bio.SearchIO.ExonerateIO',
    'Bio.SearchIO.InterproscanIO',
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

if is_jython():
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

if is_jython():
    # Jython doesn't support C extensions
    EXTENSIONS = []
elif is_ironpython():
    # Skip C extensions for now
    EXTENSIONS = []
else:
    EXTENSIONS = [
        Extension('Bio.Align._aligners',
                  ['Bio/Align/_aligners.c']),
        Extension('Bio.cpairwise2',
                  ['Bio/cpairwise2module.c']),
        Extension('Bio.Nexus.cnexus',
                  ['Bio/Nexus/cnexus.c']),
        Extension('Bio.PDB.QCPSuperimposer.qcprotmodule',
                  ["Bio/PDB/QCPSuperimposer/qcprotmodule.c"]),
        Extension('Bio.motifs._pwm',
                  ["Bio/motifs/_pwm.c"]),
        Extension('Bio.Cluster._cluster',
                  ['Bio/Cluster/cluster.c', 'Bio/Cluster/clustermodule.c']),
        Extension('Bio.PDB.kdtrees',
                  ["Bio/PDB/kdtrees.c"]),
        Extension('Bio.KDTree._CKDTree',
                  ["Bio/KDTree/KDTree.c", "Bio/KDTree/KDTreemodule.c"]),
        ]
    if not is_pypy():
        # Bio.trie has a problem under PyPy2 v5.6 and 5.7
        EXTENSIONS.extend([
                Extension('Bio.trie',
                          ['Bio/triemodule.c', 'Bio/trie.c'],
                          include_dirs=["Bio"]),
                ])


# We now define the Biopython version number in Bio/__init__.py
# Here we can't use "import Bio" then "Bio.__version__" as that would
# tell us the version of Biopython already installed (if any).
__version__ = "Undefined"
for line in open('Bio/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

# We now load in our reStructuredText README.rst file to pass
# explicitly in the metadata since at time of writing PyPI
# did not do this for us.
#
# Without declaring an encoding, if there was a problematic
# character in the file, it would work on Python 2 but might
# fail on Python 3 depending on the user's locale. By explicitly
# checking ASCII (could use latin1 or UTF8 if needed later),
# if any invalid character does appear in our README, this will
# fail and alert us immediately on either platform.
with open("README.rst", "rb") as handle:
    # Only Python 3's open has an encoding argument.
    # Opening in binary and doing decoding like this to work
    # on both Python 2 and 3.
    readme_rst = handle.read().decode("ascii")

setup(name='biopython',
      version=__version__,
      author='The Biopython Contributors',
      author_email='biopython@biopython.org',
      url='https://biopython.org/',
      description='Freely available tools for computational molecular biology.',
      long_description=readme_rst,
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: Freely Distributable',
          # Technically the "Biopython License Agreement" is not OSI approved,
          # but is almost https://opensource.org/licenses/HPND so might put:
          # 'License :: OSI Approved',
          # To resolve this we are moving to dual-licensing with 3-clause BSD:
          # 'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      cmdclass={
          "install": install_biopython,
          "build_py": build_py_biopython,
          "build_ext": build_ext_biopython,
          "test": test_biopython,
      },
      packages=PACKAGES,
      ext_modules=EXTENSIONS,
      include_package_data=True,  # done via MANIFEST.in under setuptools
      install_requires=REQUIRES,
      )
