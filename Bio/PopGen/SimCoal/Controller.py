# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.

"""
This module allows to control Simcoal2.

"""

import os
import sys
import tempfile
from shutil import copyfile
from logging import debug

class SimCoalController:
    def __init__(self, simcoal_dir):
        """Initializes the controller.
        
        simcoal_dir is the directory where simcoal is.
        
        The initializer checks for existance and executability of binaries.
        """
        self.simcoal_dir = simcoal_dir
        self.os_name = os.name
        if self.os_name=='nt' or sys.platform=='cygwin':
            self.bin_name = 'simcoal2.exe'
            #this is wrong (the exe name), most probably
        else:
            self.bin_name = 'simcoal2'
            #This name is too specific
        dir_contents = os.listdir(self.simcoal_dir)
        if self.bin_name in dir_contents:
            if not os.access(self.simcoal_dir + os.sep +
                self.bin_name, os.X_OK):
                raise IOError, "SimCoal not executable"
        else:
            raise IOError, "SimCoal not available"

    def run_simcoal(self, par_file, num_sims, ploydi = '1', par_dir = '.'):
        """Executes SimCoal.
        """
        if par_dir == None:
            par_dir = os.sep.join([Config.dataDir, 'SimCoal', 'runs'])
        curr_dir = os.getcwd()
        os.chdir(par_dir)
        os.system(self.simcoal_dir + os.sep + self.bin_name + ' ' +
          par_file + ' ' + str(num_sims) + ' ' + ploydi + ' >/dev/null 2>&1')
        os.chdir(curr_dir)
    

