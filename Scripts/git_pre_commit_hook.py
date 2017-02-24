#!/usr/bin/env python
"""Biopython's pre-commit git hook.

This is intended to be installed as part of setting up a development
copy of Biopyton, e.g.::

    $ git clone git@github.com:biopython/biopython.git
    $ cd biopython
    $ ln -s Scripts/git_pre_commit_hook.py .git/hooks/pre-commit

This script will by default (such as when run by git when you run the
``git commit`` command) determine all new of changed files according
to git, and apply our style checks to them.

The script can alternatively be run with filenames or folders to check,
in which case the git status of the files is ignored. This can be used
for example when working on this script itself (e.g. as part of work to
tighten the style checks), or for continous integration testing via
TravisCI and/or Tox.
"""
from __future__ import with_statement, print_function

import os
import re
import shutil
import subprocess
import sys
import tempfile

try:
    from subprocess import getoutput
except ImportError:
    # Python 2
    from commands import getoutput

# These assume we can append one or more filenames/directories to the end...
python_check_templates = (
    ("pep8", "--ignore", "E121,E122,E123,E124,E125,E126,E127,E128,E129,E131,E501"),
    ("flake8", "--ignore", "E121,E122,E123,E124,E125,E126,E127,E128,E129,E131,E501"),
    ("pydocstyle", "--ignore", "D100,D101,D102,D103,D104,D105,D200,D203,D204,D205,D207,D208,D209,D210,D211,D212,D213,D301,D302,D400,D401,D402,D403,D404"),
)

# TODO - Implement a proper API for this?
if len(sys.argv) > 1:
    # assume called with list of local directories and/or files to lint
    files = tuple(sys.argv[1:])
    tempdir = None
    cwd = "."
else:
    # Assume called via 'git commit' and ask git for changed/new files
    modified = re.compile("^[AM]+\s+(?P<name>.*\.py)", re.MULTILINE)
    files = getoutput("git status --porcelain")
    files = tuple(modified.findall(files))
    if not files:
        sys.stderr.write("Seems nothing has changed which needs linting...\n")
        sys.exit(0)
    elif len(files) <= 5:
        sys.stderr.write("Linting files: %s\n" % " ".join(files))
    else:
        sys.stderr.write("Linting %i files...\n" % len(files))

    # Prepare temp directory of the files as staged in git
    # Can't just lint in place as would see changes pending git add
    tempdir = tempfile.mkdtemp()
    for name in files:
        filename = os.path.join(tempdir, name)
        filepath = os.path.dirname(filename)

        if not os.path.exists(filepath):
            os.makedirs(filepath)
        with open(filename, "w") as f:
            subprocess.check_call(("git", "show", ":" + name), stdout=f)
    cwd = tempdir

# We run all the checks up front so user doesn't get a false sense of
# what is wrong if we just showed problems from the first too.
failed = False
for cmd_template in python_check_templates:
    child = subprocess.Popen(cmd_template + files, cwd=cwd)
    child.communicate()
    return_code = child.returncode
    if return_code:
        failed = True

if tempdir:
    shutil.rmtree(tempdir)
if failed:
    sys.stderr.write("Style validation failed. Please fix, or force with 'git commit --no-verify'\n")
    sys.exit(1)
