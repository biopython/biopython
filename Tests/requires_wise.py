#!/usr/bin/env python

# Copyright 2004 by Michael Hoffman.  All rights reserved.  This code is
# part of the Biopython distribution and governed by its license.
# Please see the LICENSE file that should have been included as part
# of this package.

import commands
import re

if commands.getoutput("dnal").find("command not found") != -1:
    raise ImportError
