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

import sys
import os

try:
    from setuptools import setup
    from setuptools import Command
    from setuptools import Extension
except ImportError:
    sys.exit(
        "We need the Python library setuptools to be installed. "
        "Try running: python -m ensurepip"
    )

if "bdist_wheel" in sys.argv:
    try:
        import wheel  # noqa: F401
    except ImportError:
        sys.exit(
            "We need both setuptools AND wheel packages installed "
            "for bdist_wheel to work. Try running: pip install wheel"
        )


# Make sure we have the right Python version.
MIN_PY_VER = (3, 7)
if sys.version_info[:2] < MIN_PY_VER:
    sys.stderr.write(
        ("ERROR: Biopython requires Python %i.%i or later. " % MIN_PY_VER)
        + ("Python %d.%d detected.\n" % sys.version_info[:2])
    )
    sys.exit(1)


class test_biopython(Command):
    """Run all of the tests for the package.

    This is a automatic test run class to make distutils kind of act like
    perl. With this you can do:

    python setup.py build
    python setup.py install
    python setup.py test

    """

    description = "Automatically run the test suite for Biopython."
    user_options = [("offline", None, "Don't run online tests")]

    def initialize_options(self):
        """No-op, initialise options."""
        self.offline = None

    def finalize_options(self):
        """No-op, finalise options."""
        pass

    def run(self):
        """Run the tests."""
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("Tests")
        sys.path.insert(0, "")
        import run_tests

        if self.offline:
            run_tests.main(["--offline"])
        else:
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
REQUIRES = ["numpy"]

# --- set up the packages we are going to install
# standard biopython packages
PACKAGES = [
    "Bio",
    "Bio.Affy",
    "Bio.Align",
    "Bio.Align.Applications",
    "Bio.Align.substitution_matrices",
    "Bio.AlignIO",
    "Bio.Alphabet",
    "Bio.Application",
    "Bio.Blast",
    "Bio.CAPS",
    "Bio.Cluster",
    "Bio.codonalign",
    "Bio.Compass",
    "Bio.Data",
    "Bio.Emboss",
    "Bio.Entrez",
    "Bio.ExPASy",
    "Bio.GenBank",
    "Bio.Geo",
    "Bio.Graphics",
    "Bio.Graphics.GenomeDiagram",
    "Bio.HMM",
    "Bio.KEGG",
    "Bio.KEGG.Compound",
    "Bio.KEGG.Enzyme",
    "Bio.KEGG.Gene",
    "Bio.KEGG.Map",
    "Bio.PDB.mmtf",
    "Bio.KEGG.KGML",
    "Bio.Medline",
    "Bio.motifs",
    "Bio.motifs.applications",
    "Bio.motifs.jaspar",
    "Bio.Nexus",
    "Bio.NMR",
    "Bio.Pathway",
    "Bio.Pathway.Rep",
    "Bio.PDB",
    "Bio.phenotype",
    "Bio.PopGen",
    "Bio.PopGen.GenePop",
    "Bio.Restriction",
    "Bio.SCOP",
    "Bio.SearchIO",
    "Bio.SearchIO._legacy",
    "Bio.SearchIO._model",
    "Bio.SearchIO.BlastIO",
    "Bio.SearchIO.HHsuiteIO",
    "Bio.SearchIO.HmmerIO",
    "Bio.SearchIO.ExonerateIO",
    "Bio.SearchIO.InterproscanIO",
    "Bio.SeqIO",
    "Bio.SeqUtils",
    "Bio.Sequencing",
    "Bio.Sequencing.Applications",
    "Bio.SVDSuperimposer",
    "Bio.PDB.QCPSuperimposer",
    "Bio.SwissProt",
    "Bio.TogoWS",
    "Bio.Phylo",
    "Bio.Phylo.Applications",
    "Bio.Phylo.PAML",
    "Bio.UniGene",
    "Bio.UniProt",
    "Bio.Wise",
    # Other top level packages,
    "BioSQL",
]

EXTENSIONS = [
    Extension("Bio.Align._aligners", ["Bio/Align/_aligners.c"]),
    Extension("Bio.cpairwise2", ["Bio/cpairwise2module.c"]),
    Extension("Bio.Nexus.cnexus", ["Bio/Nexus/cnexus.c"]),
    Extension(
        "Bio.PDB.QCPSuperimposer.qcprotmodule",
        ["Bio/PDB/QCPSuperimposer/qcprotmodule.c"],
    ),
    Extension("Bio.motifs._pwm", ["Bio/motifs/_pwm.c"]),
    Extension(
        "Bio.Cluster._cluster", ["Bio/Cluster/cluster.c", "Bio/Cluster/clustermodule.c"]
    ),
    Extension("Bio.PDB.kdtrees", ["Bio/PDB/kdtrees.c"]),
    Extension("Bio.PDB.ccealign", ["Bio/PDB/ccealignmodule.c"]),
    Extension("Bio.SeqIO._twoBitIO", ["Bio/SeqIO/_twoBitIO.c"]),
]

# We now define the Biopython version number in Bio/__init__.py
# Here we can't use "import Bio" then "Bio.__version__" as that would
# tell us the version of Biopython already installed (if any).
__version__ = "Undefined"
for line in open("Bio/__init__.py"):
    if line.startswith("__version__"):
        exec(line.strip())

# We now load in our reStructuredText README.rst file to pass explicitly in the
# metadata, since at time of writing PyPI did not do this for us.
#
# Must make encoding explicit to avoid any conflict with the local default.
# Currently keeping README as ASCII (might switch to UTF8 later if needed).
# If any invalid character does appear in README, this will fail and alert us.
with open("README.rst", encoding="ascii") as handle:
    readme_rst = handle.read()

setup(
    name="biopython",
    version=__version__,
    author="The Biopython Contributors",
    author_email="biopython@biopython.org",
    url="https://biopython.org/",
    description="Freely available tools for computational molecular biology.",
    long_description=readme_rst,
    project_urls={
        "Documentation": "https://biopython.org/wiki/Documentation",
        "Source": "https://github.com/biopython/biopython/",
        "Tracker": "https://github.com/biopython/biopython/issues",
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: Freely Distributable",
        # Technically the "Biopython License Agreement" is not OSI approved,
        # but is almost https://opensource.org/licenses/HPND so might put:
        # 'License :: OSI Approved',
        # To resolve this we are moving to dual-licensing with 3-clause BSD:
        # 'License :: OSI Approved :: BSD License',
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    cmdclass={"test": test_biopython},
    packages=PACKAGES,
    ext_modules=EXTENSIONS,
    include_package_data=True,  # done via MANIFEST.in under setuptools
    install_requires=REQUIRES,
    python_requires=">=%i.%i" % MIN_PY_VER,
)
