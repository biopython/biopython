# Allow files in this directory to be imported

import warnings
from Bio import BiopythonDeprecationWarning
warnings.warn("Bio.NeuralNetwork has been deprecated, and we intend to remove"
              " it in a future release of Biopython. Please consider using"
              " scikit-learn or TensorFlow instead.  If you would like to"
              " continue using Bio.SomeModule, please contact the Biopython"
              " developers via the mailing list or GitHub.",
              BiopythonDeprecationWarning)
