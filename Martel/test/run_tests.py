# test suite for the formats
import os, glob, imp, sys

def main(args):
    try:
        dirname = os.path.dirname(__file__)
    except NameError:
        dirname = "."
    files = glob.glob(os.path.join(dirname, "test_*.py"))
    files.sort()

    for file in files:
        name = os.path.basename(file)
        name = os.path.splitext(name)[0]
        print "###### running tests in", name, "###############"
        module = __import__(name)
        module.test()

# Run tests uses the local (uninstalled) files
def local_test_main(args):
    # Find the directory containing the Martel code.
    # In biopython CVS this is named 'Martel'.  In a standalone release
    # it could be named something like 'Martel-0.8'.  Discover which.
    sys.path.insert(0, "..")
    filename = os.path.abspath(__file__)
    dirname = os.path.dirname(os.path.dirname(filename))

    # Force module loading.
    imp.load_module("Martel", open("../__init__.py"), dirname,
                    (".py", "r", imp.PKG_DIRECTORY))

    # Do the normal tests
    return main(args)

# Run tests uses the installed files (not the ones in ..)
install_test_main = main

if __name__ == "__main__":
    install_test_main([])
