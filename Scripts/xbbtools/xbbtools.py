#!/usr/bin/env python
# Copyright 2000 by Thomas Sicheritz-Ponten.
# Copyrigth 2016 by Markus Piotrowski.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Created: Thu Jul 13 11:09:00 2000
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas

"""Main code entry point for graphical Xbbtools tool."""

import sys

try:
    import Tkinter as tk  # Python 2
except ImportError:
    import tkinter as tk  # Python 3

from xbb_widget import xbb_widget


win = tk.Tk()
win.title('xbb tools')
xbbtools = xbb_widget()

try:
    xbbtools.open(sys.argv[1])
except IndexError:  # Started script without specifying a filename
    pass

win.mainloop()
