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

class SimCoalController(object):
    def __init__(self, simcoal_dir):
        """Initializes the controller.
        
        simcoal_dir is the directory where simcoal is.
        
        The initializer checks for existance and executability of binaries.
        """
        self.simcoal_dir = simcoal_dir
        self.os_name = os.name #remove this?
        dir_contents = os.listdir(self.simcoal_dir)
        #We expect the tool to be installed as simcoal2(.exe)
        #without any trailing version number.
        self.bin_name = "simcoal2"
        if self.bin_name not in dir_contents:
            #Try case insensitive,
            dir_contents = [x.lower() for x in dir_contents]
        if self.bin_name not in dir_contents:
            #Try with .exe
            self.bin_name += '.exe'
        if self.bin_name not in dir_contents:
            raise IOError("SimCoal not available")
        if not os.access(os.path.join(self.simcoal_dir, self.bin_name),
                         os.X_OK):
            raise IOError("SimCoal not executable")

    def run_simcoal(self, par_file, num_sims, ploydi = '1', par_dir = '.'):
        """Executes SimCoal.
        """
        if par_dir == None:
            par_dir = os.sep.join([".", 'SimCoal', 'runs'])
        curr_dir = os.getcwd()
        #TODO - Make sure we change drive on Windows as well?
        os.chdir(par_dir)
        exe = os.path.join(self.simcoal_dir, self.bin_name)
        if " " in exe:
            exe = '"' + exe + '"'
        cmd = exe + ' ' + par_file + ' ' + str(num_sims) + ' ' + ploydi
        #TODO - Better way to spot if on Jython on Windows?
        if sys.platform=="win32" or self.bin_name.endswith(".exe"):
            #There is no /dev/nul on Windows
            cmd += ' > nul 2>nul'
        else:
            cmd += ' >/dev/null 2>&1'
        os.system(cmd)
        os.chdir(curr_dir)
    

