# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Deal with various biological databases and services on the web (DEPRECATED).

The functionality within Bio.WWW and its sub modules was deprecated in
Biopython 1.45, however most individual functions simply been moved:

Bio.WWW.ExPASy -> Bio.ExPASy
Bio.WWW.InterPro -> Bio.InterPro
Bio.WWW.NCBI -> Bio.Entrez
Bio.WWW.SCOP -> Bio.SCOP
"""
import time
import warnings
warnings.warn("Bio.WWW was deprecated.  Most of its functionality is now "\
              +"available from Bio.ExPASy, Bio.InterPro, Bio.Entrez and "\
              +"Bio.SCOP.")

class RequestLimiter:
    # This class implements a simple countdown timer for delaying WWW
    # requests.
    def __init__(self, delay):
        self.last_time = 0.0
        self.delay = delay
    def wait(self, delay=None):
        if delay is None:
            delay = self.delay
        how_long = self.last_time + delay - time.time()
        if how_long > 0:
            time.sleep(how_long)
        self.last_time = time.time()
