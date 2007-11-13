# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.

"""
This module allows to control Simcoal2.

"""

import os
import tempfile
from shutil import copyfile
from logging import debug

from PopGen import Config

class SimCoalController:
    def __init__(self, simcoalDir = None):
        """Initializes the controller.
        
        simcoalDir is the directory where simcoal is.
        
        The initializer checks for existance and executability of binaries.
        """
        if simcoalDir == None:
            self.simcoalDir = Config.simcoalDir
        else:
            self.simcoalDir = simcoalDir
        self.os_name = os.name
        if self.os_name=='nt':
            self.bin_name = 'simcoal.exe'
            #this is wrong (the exe name), most probably
        else:
            self.bin_name = 'simcoal2_1_2'
            #This name is too specific
        dir_contents = os.listdir(self.simcoalDir)
        if self.bin_name in dir_contents:
            if not os.access(self.simcoalDir + os.sep +
                self.bin_name, os.X_OK):
                raise IOError, "SimCoal not executable"
        else:
            raise IOError, "SimCoal not available"

    def run_simcoal(self, par_file, num_sims, ploydi = '1', parDir = None):
        """Executes SimCoal.
        """
        if parDir == None:
            parDir = os.sep.join([Config.dataDir, 'SimCoal', 'runs'])
        curr_dir = os.getcwd()
        os.chdir(parDir)
        os.system(self.simcoalDir + os.sep + self.bin_name + ' ' +
          par_file + ' ' + str(num_sims) + ' ' + ploydi + ' >/dev/null 2>&1')
        os.chdir(curr_dir)
    

