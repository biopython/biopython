"""Distutils based setup script for Martel.

This uses Distutils (http://python.org/sigs/distutils-sig/) the
standard python mechanism for installing packages. For the easiest
installation just type the command:

  python setup.py install

For more details about the options available from distutils, look at
the 'Installing Python Modules' distutils documentation, available
from:

  http://python.org/sigs/distutils-sig/doc/

Or, if all else fails, feel free to write to the biopython list
at biopython@biopython.org and ask for help.
"""

# This setup.py is a modified vesion of the standard Biopython
# setup.py, which explains why some of the code is overkill :)

import sys, os

try:
    from distutils.core import setup
    from distutils.command.install import install
    from distutils.core import Command
except ImportError:
    print "Martel installation requires distutils, available with python 2.0"
    print "or better, or from:"
    print "  http://python.org/sigs/distutils-sig/download.html"
    sys.exit(0)


# --- check for installed programs needed by Martel

def check_install(name, check_library, location, other_messages = None):
    """Check if a program is installed and print a warning message if not.

    This helps users at least know they are missing some installed stuff
    and where to get it when they install Martel.

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
        print "because otherwise Martel will not be useful"
        if other_messages:
            print other_messages

# -- functions to check for specific libraries and versions.

def check_mxTextTools():
    try:
        from mx import TextTools
        return 1
    except ImportError:
        pass
    try:
        import TextTools
        return 1
    except ImportError:
        pass
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
        check_install("mxTextTools", check_mxTextTools,
                      "http://www.lemburg.com/files/python/mxExtensions.html")

class run_local_tests(Command):
    """Run all of the tests for the package using uninstalled (local) files

    This is a automatic test run class to make distutils kind of act like
    perl. With this you can do:

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
        os.chdir("test")
        import run_tests
        run_tests.local_test_main([])

        # change back to the current directory
        os.chdir(this_dir)


class run_install_tests(run_local_tests):
    """Run all of the tests for the package using installed files

    This is a automatic test run class to make distutils kind of act like
    perl. With this you can do:

    python setup.py install
    python setup.py installtest
    
    """
    def run(self):
        this_dir = os.getcwd()

        # change to the test dir and run the tests
        os.chdir("test")
        import run_tests
        run_tests.install_test_main([])

        # change back to the current directory
        os.chdir(this_dir)

setup(name = "Martel",
      version = "1.43",
      description = "Parse flat-file formats as if they are in XML",
      author = "Dalke Scientific Software, LLC; " \
               "member of the The Biopython Consortium",
      author_email = "dalke@dalkescientific.com",
      url = "http://www.dalkescientific.com/Martel",
      license = "Biopython license",

      cmdclass = {"install" : my_install,
                  "test" : run_local_tests,
                  "installtest" : run_install_tests,
                  },
      package_dir = {"Martel": ""},
      packages = ["Martel"],
      )

