#!/usr/bin/env python

import glob, sys, os

try:
    import EUtils
except ImportError:
    raise SystemExit("""\
The EUtils package is not on your PYTHONPATH.

Here are a couple ways to test the EUtils distribution
before doing an install.

 1. make a symbolic link to the distribution directory
 
      ln -s ../EUtils .

    This is what I do for my testing.

 2. Use env to set the PYTHONPATH environment variable for
    this test

     env PYTHONPATH=.. python runtests.py

     or if you already have a PYTHONPATH

     env PYTHONPATH=..:$PYTHONPATH python runtests.py

If you are on a non-Unix-like system, consult the documentation for
setting the PYTHONPATH.

You should not need either of these options to test an installed
version of EUtils.
""")


filenames = glob.glob("test_*.py")
error_count = 0
run_count = 0
for filename in filenames:
    print "Running tests in:", filename
    retval = os.system("%s %s" % (sys.executable, filename))
    error_count = error_count + (retval != 0)
    run_count = run_count + 1
    
if error_count:
    sys.exit("Test failures: %d of %d" % (error_count, run_count))

print "All %d regression test files ran successfully." % (run_count,)
