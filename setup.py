# XXX message about how distutils is required.
# how to run distutils to install the package

import sys
try:
    from distutils.core import setup
except ImportError:
    # XXX improve this error message
    print "Requires distutils"
    sys.exit(0)

setup(name='biopython', 
      version='0.90-d01',
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
      ext_modules = [('Bio.Tools.Classification.cSVM',
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
      )

