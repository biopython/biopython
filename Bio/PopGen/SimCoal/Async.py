# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.

"""
This module allows to cache Simcoal2 results, and return on the fly
in case the calculation was done. Async version

This version will run Sincoal2 (if necessary) Asynchrously.

"""

from logging import debug
from sys import exit
import os
import tarfile
import tempfile
from PopGen.SimCoal.Controller import SimCoalController
from PopGen.SimCoal import Cache
from PopGen import Config
from PopGen import Async

class SimCoalCache(Cache.SimCoalCache):
    def __init__(self, simcoalDir = None):
        Cache.SimCoalCache.__init__(self, simcoalDir)

    def runJob(self, parameters, inputFiles):
        parFile = parameters['parFile']
        numSims = parameters['numSims']
        ploydi = parameters.get('ploydi', '1')
        f = inputFiles[parFile]
        text = f.read()
        f.close()
        w = open (os.sep.join([Config.dataDir, 'SimCoal', 'runs', parFile]), 'w')
        w.write(text)
        w.close()
        self.run_simcoal(parFile, numSims, ploydi)
        return 0, None
