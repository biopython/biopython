# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import os
import sys
from Bio import MissingExternalDependencyError

def is_exe(filepath):
        return os.path.exists(filepath) and os.access(filepath, os.X_OK)

def which(program):
    filepath, filename = os.path.split(program)
    os_path = os.environ["PATH"].split(os.pathsep)
    if sys.platform == "win32":
        try:
            #This can vary depending on the Windows language.
            prog_files = os.environ["PROGRAMFILES"]
        except KeyError:
            prog_files = r"C:\Program Files"
        #For Windows, the user is instructed to move the programs to a folder
        #and then to add the folder to the system path. Just in case they didn't
        #do that, we can check for it in Program Files.
        likely_dirs = ["", #Current dir
                       prog_files,
                       os.path.join(prog_files,"paml41"),
                       os.path.join(prog_files,"paml43"),
                       os.path.join(prog_files,"paml44")] + sys.path
        os_path.extend(likely_dirs)
    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None

if sys.platform == "win32":
    binaries = ["codeml.exe", "baseml.exe", "yn00.exe"]
else:
    binaries = ["codeml", "baseml", "yn00"]
for binary in binaries:    
    if which(binary) is None:
        raise MissingExternalDependencyError(\
            "Install PAML if you want to use the Bio.Phylo.PAML wrapper.")

