# Copyright 2002 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import Martel as _M   # make sure not imported into client's namespaces
import string as _string

blank_expr = _M.AssertNot(_M.Re("."))
html_expr = _M.Rep(_M.Any(_string.whitespace)) + _M.NoCase(_M.Str("<html>"))

def has_expr(expr):
    return _M.Rep(_M.AssertNot(expr) + _M.Alt(_M.Re("."), _M.AnyEol())) + expr

def has_str(str):
    return has_expr(_M.Str(str))
