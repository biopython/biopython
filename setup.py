"""Distutils based setup script for Biopython.

This uses Distutils (http://python.org/sigs/distutils-sig/) the standard
python mechanism for installing packages. For the easiest installation
just type the command:

python setup.py install

For more in-depth instructions, see the installation section of the
biopython manual, linked to from:

http://biopython.org/wiki/html/BioPython/BiopythonCode.html

Or for more details about the options available from distutils, look at
the 'Installing Python Modules' distutils documentation, available from:

http://python.org/sigs/distutils-sig/doc/

Or, if all else, fails, feel free to write to the biopython list
(biopython@biopython.org) and ask for help.
"""

import sys
import os
try:
    from distutils.core import setup
except ImportError:
    print "Biopython installation requires distutils, avaiable with python 2.0"
    print "or from http://python.org/sigs/distutils-sig/download.html"
    sys.exit(0)

# check if the distutils has the new extension class stuff
# this is to support old distutils which do extensions differently
try:
    from distutils.extension import Extension
except ImportError:
    print "Your version of distutils is really old. You need to upgrade"
    print "to a newer version. The latest releases of distutils are available"
    print "from http://python.org/sigs/distutils-sig/download.html"
    sys.exit(0)

setup(name='biopython', 
      version='0.90d04',
      author='The Biopython Consortium',
      author_email='biopython@biopython.org',
      url='http://www.bipoython.org/',
      
      packages=['Bio',
                'Bio.Align',
                'Bio.Alphabet',
                'Bio.Blast',
                'Bio.Clustalw',
                'Bio.Data',
                'Bio.Encodings',
                'Bio.Entrez',
                'Bio.Enzyme',
                'Bio.Fasta',
                'Bio.GenBank',
                'Bio.Gobase',
                'Bio.Medline',
                'Bio.PDB',
                'Bio.Prosite',
                'Bio.Rebase',
                'Bio.SCOP',
                'Bio.SeqIO',
                'Bio.SubsMat',
                'Bio.SwissProt',
                'Bio.Tools',
                'Bio.Tools.Classification',
                'Bio.Tools.Parsers',
                'Bio.WWW'
                ],
      
      ext_modules = [Extension('Bio.Tools.Classification.cSVM',
                               ['Bio/Tools/Classification/cSVMmodule.c']
                               ),
                     Extension('Bio.Tools.clistfns',
                               ['Bio/Tools/clistfnsmodule.c']
                               ),
                     Extension('Bio.Tools.cmathfns',
                               ['Bio/Tools/cmathfnsmodule.c']
                               ),
                     Extension('Bio.Tools.cstringfns',
                               ['Bio/Tools/cstringfnsmodule.c']
                               )
                     ]
      )

