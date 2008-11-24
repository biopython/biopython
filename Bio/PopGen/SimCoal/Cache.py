# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.

"""
This module allows to cache Simcoal2 results, and return on the fly
in case the calculation was done.

"""

from logging import debug
from sys import exit
import os
import tarfile
import tempfile
import SimCoalController
import Config

class SimCoalCache:
    def __init__(self, simcoalDir = None):
        """Initializes the cache.
        
        simcoalDir - where the binaries are

        IMPORTANT: The cache only makes sense if the file nameunivically
        identifies the model.
          For now use use the model name as key,
          and it will probably stay like that.
        """
        self.cacheDir = os.sep.join([Config.dataDir, 'SimCoal', 'cache'])
        if simcoalDir == None:
            self.simcoalDir = Config.simcoalDir
        else:
            self.simcoalDir = simcoalDir


    def run_simcoal(self, par_file, num_sims, ploydi = '1', parDir = None):
        if parDir == None:
            parDir = os.sep.join([Config.dataDir, 'SimCoal', 'runs'])
        par_file_root = par_file[:-4]
        tar_name = os.sep.join([self.cacheDir, ploydi, par_file_root + 
            '.tar.bz2'])
        if os.access(tar_name, os.R_OK):
            tf = tarfile.open(tar_name)
            tar_num_sims = len(tf.getmembers()) - 3
        else:
            tar_num_sims = 0
        if tar_num_sims >= num_sims:
            tf.extractall(parDir)
            tf.close()
            return
        else:
            try:
                tf.close()
            except NameError:
                pass #not opened in the first place, OK.
        scc = SimCoalController(self.simcoalDir)
        scc.run_simcoal(par_file, num_sims, ploydi, parDir)
        tf = tarfile.open(tar_name, 'w:bz2')
        tf.add(os.sep.join([parDir, par_file_root]), par_file_root)
        tf.close()

    def listSimulations(self, ploidy = '1'):
        '''
           Lists available simulations.
        '''
        files = os.listdir(self.cacheDir + os.sep + ploidy)
        sims = []
        for file in files:
            if file.endswith('.tar.bz2'):
                sims.append(file[:-8])
        return sims

    def getSimulation(self, sim_name, ploidy = '1', parDir = None):
        '''
           Makes available a cached simulation.

           @param sim_name simulation name.

           This mainly means untaring a file.
        '''
        if parDir == None:
            parDir = os.sep.join([Config.dataDir, 'SimCoal', 'runs'])
        tar_name = os.sep.join([self.cacheDir, ploidy, sim_name + 
            '.tar.bz2'])
        tf = tarfile.open(tar_name)
        tf.extractall(parDir)
        tf.close()


if __name__ == '__main__':
  cache = Cache('/home/work/werk/consolidator/sc_cache',
      '/home/work/software/simcoal')
  cache.run_simcoal('.', 'island_snp-50_0.0025_10_0.083_100_60.par', 102)
