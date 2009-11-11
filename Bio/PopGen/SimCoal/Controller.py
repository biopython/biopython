# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
This module allows you to control Simcoal2.

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
        dir_contents = os.listdir(self.simcoal_dir)
        #We expect the tool to be installed as simcoal2(.exe)
        #without any trailing version number.
        if self.os_name=='nt' or sys.platform=='cygwin':
            self.bin_name = 'simcoal2.exe'
            #Windows is case insenstive
            dir_contents = [x.lower() for x in dir_contents]
        else:
            self.bin_name = 'simcoal2'
        if self.bin_name in dir_contents:
            if not os.access(self.simcoal_dir + os.sep +
                self.bin_name, os.X_OK):
                raise IOError("SimCoal not executable")
        else:
            raise IOError("SimCoal not available")

    def run_simcoal(self, par_file, num_sims, ploydi = '1', par_dir = '.'):
        """Executes SimCoal.
        """
        if par_dir == None:
            par_dir = os.sep.join([Config.dataDir, 'SimCoal', 'runs'])
        curr_dir = os.getcwd()
        #TODO - Make sure we change drive on Windows as well?
        os.chdir(par_dir)
        cmd = self.simcoal_dir + os.sep + self.bin_name + ' ' + \
              par_file + ' ' + str(num_sims) + ' ' + ploydi
        if sys.platform=="win32":
            #There is no /dev/nul on Windows
            cmd += ' > nul 2>nul'
        else:
            cmd += ' >/dev/null 2>&1'
        os.system(cmd)
        os.chdir(curr_dir)
    

