#!/usr/bin/env python

"""Generate documentation for the biopython libraries.

This uses HappyDoc (http://sourceforge.net/projects/happydoc) to do the
documentation.

Usage:
* Set your PYTHONPATH to include the HappyDoc directory.
* Run this file from the Doc directory as './make_api_docs.py'

Right now this is kludgy and rough and needs work."""

IMPORT_MESSAGE = \
"""***Error: make_docs.py requires HappyDoc to be on your PYTHONPATH.
You can get HappyDoc from http://sourceforge.net/projects/happydoc."""

import os
import sys

# change to the main biopython directory, it makes the docs look nicer if we
# run the program from there
if os.path.basename(os.getcwd()) == 'Doc':
    os.chdir(os.pardir)

# Kludge: manipulate sys.argv so it will look like happydoc is being called
# via the happydoc.py program.

try:
    # this is a HappyDoc file. We can't import happydoc_class or anything
    # that imports it until after we do this sys.argv hack.
    import StructuredText
except ImportError:
    print IMPORT_MESSAGE
    sys.exit(1)

# find where happydoc is via the imported StructuredText file
happydoc_location = os.path.dirname(StructuredText.__file__)
print 'happydoc_location:', happydoc_location

sys.argv = [os.path.join(happydoc_location, 'happydoc.py'), 'Bio']

try:
    from happydoc_class import HappyDoc
except ImportError:
    print IMPORT_MESSAGE
    sys.exit()


if __name__ == '__main__':
    
    try:
        doc_gen = HappyDoc()

        # set parameters for Biopython
        doc_gen.author_name = 'The Biopython folks (biopython@biopython.org)'
        doc_gen.output_directory = 'Doc/api'
        doc_gen.app_home = 'http://www.biopython.org'
        doc_gen.docset_title = 'Biopython API documentation'
        doc_gen.package_description_file = 'README'
        
        doc_gen.run()
    except HappyDoc.HelpRequested:
        pass
    
