# XXX message about how distutils is required.
# how to run distutils to install the package

import sys
try:
    from distutils.core import setup
except ImportError:
    # XXX improve this error message
    print "Requires distutils"
    sys.exit(0)

# check if the distutils has the new extension class stuff
# this is to support old distutils which do extensions differently
try:
    from distutils.extension import Extension
    new_extension = 1
except ImportError:
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
      version='0.90-d02',
      author='The Biopython Consortium',
      author_email='biopython@biopython.org',
      url='http://www.bipoython.org/',
      
      packages=['Bio',
                'Bio/SwissProt',
                'Bio/Blast',
                'Bio/Enzyme',
                'Bio/Fasta',
                'Bio/PDB',
                'Bio/Prosite',
                'Bio/Medline',
                'Bio/Data',
                'Bio/Entrez',
                'Bio/Tools', 'Bio/Tools/Classification',
                'Bio/WWW',
                'Bio/Alphabet',
                'Bio/Encodings',
                'Bio/SeqIO'
                ],

      ext_modules = extensions
      )

