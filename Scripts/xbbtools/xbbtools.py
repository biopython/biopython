#!/usr/bin/env python
# Created: Thu Jul 13 11:09:00 2000
# Last changed: Time-stamp: <00/12/02 15:57:45 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xbbtools.py

# Copyright 2000 by Thomas Sicheritz-Ponten.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import sys
sys.path.insert(0, '.')
from Tkinter import *

from xbb_widget import xbb_widget
win = Tk()

xbbtools = xbb_widget()
xbbtools.main_frame.option_add('*frame.background', 'dimgrey')

try:
    xbbtools.open(sys.argv[1])
except:
    pass

win.mainloop()

