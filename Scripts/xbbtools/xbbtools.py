#!/usr/bin/env python
# Created: Thu Jul 13 11:09:00 2000
# Last changed: Time-stamp: <00/08/08 23:32:18 thomas>
# Thomas.Sicheritz@molbio.uu.se, http://evolution.bmc.uu.se/~thomas
# File: xbbtools.py

import string, re, regsub
import posixpath, posix
import os, sys  # os.system, sys.argv
sys.path.insert(0, '.')
from Tkinter import *

from xbb_widget import xbb_widget
win = Tk()

xbbtools = xbb_widget()
xbbtools.main_frame.option_add('*frame.background', 'dimgrey')

xbbtools.open(sys.argv[1])

win.mainloop()

