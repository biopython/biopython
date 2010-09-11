"""Third party and other parsers useful internally to Biopython (DEPRECATED).
"""
import warnings
import Bio
warnings.warn("Bio.Parsers (including our copy of SPARK) and the only part "
              "of Biopython which used it, Bio.GenBank.LocationParser, are "
              "now deprecated and will be removed in a future release of "
              "Biopython. If you want to continue to use any of this code, "
              "please get in contact with the Biopython developers via "
              "the mailing lists to avoid its permanent removal from "
              "Biopython.", Bio.BiopythonDeprecationWarning)
