#!/usr/bin/env python

# Copyright 2004 by Michael Hoffman.  All rights reserved.  This code is
# part of the Biopython distribution and governed by its license.
# Please see the LICENSE file that should have been included as part
# of this package.

import commands
import re

not_found_types = ["command not found", "dnal: not found"]
dnal_output = commands.getoutput("dnal")

for not_found in not_found_types:
    if dnal_output.find(not_found) != -1:
        raise ImportError(dnal_output)
