# XXX message about how distutils is required.
# how to run distutils to install the package

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
    new_extension = 1
except ImportError:
    print "Your version of distutils is really old. You should definately"
    print "upgrade to a newer version!"
    new_extension = 0

if new_extension:
    extensions = [Extension('Bio.Tools.Classification.cSVM',
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
else:
    extensions = [('Bio.Tools.Classification.cSVM',
                   { 'sources' : ['Bio/Tools/Classification/cSVMmodule.c'] }
                   ),
                  ('Bio.Tools.clistfns',
                   { 'sources' : ['Bio/Tools/clistfnsmodule.c'] }
                   ),
                  ('Bio.Tools.cmathfns',
                   { 'sources' : ['Bio/Tools/cmathfnsmodule.c'] }
                   ),
                  ('Bio.Tools.cstringfns',
                   { 'sources' : ['Bio/Tools/cstringfnsmodule.c'] }
                   )
                  ]


setup(name='biopython', 
      version='0.90d03',
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
                'Bio.Gobase',
                'Bio.Medline',
                'Bio.PDB',
                'Bio.Prosite',
                'Bio.Rebase',
                'Bio.SCOP',
                'Bio.SeqIO',
                'Bio.SubsMat',
                'Bio.SwissProt',
                'Bio.Tools', 'Bio.Tools.Classification',
                'Bio.WWW'
                ],
      
      ext_modules = extensions
      )

