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

# find out where happydoc is by checking through the imports for it
happydoc_location = None
for import_dir in sys.path:
    if import_dir != '':
        dir_contents = os.listdir(import_dir)
        if 'happydoc.py' in dir_contents:
            happydoc_location = import_dir

if happydoc_location is None:
    print IMPORT_MESSAGE
    sys.exit()

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
        package_description_file = 'README'
        
        doc_gen.run()
    except HappyDoc.HelpRequested:
        pass
    
