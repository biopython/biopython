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

from Controller import SimCoalController
import Cache

class SimCoalCache(Cache.SimCoalCache):
    def __init__(self, data_dir, simcoal_dir):
        self.data_dir = data_dir
        Cache.SimCoalCache.__init__(self, data_dir, simcoal_dir)

    def runJob(self, parameters, inputFiles):
        parFile = parameters['parFile']
        numSims = parameters['numSims']
        ploydi = parameters.get('ploydi', '1')
        f = inputFiles[parFile]
        text = f.read()
        f.close()
        w = open (os.sep.join([self.data_dir, 'SimCoal', 'runs', parFile]), 'w')
        w.write(text)
        w.close()
        self.run_simcoal(parFile, numSims, ploydi)
        return 0, None
