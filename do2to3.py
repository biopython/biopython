"""Helper script for building and installing Biopython on Python 3.

Note that we can't just use distutils.command.build_py function build_py_2to3
in setup.py since (as far as I can see) that does not allow us to alter the
2to3 options. In particular, we need to turn off the long fixer for some of
our files.

This code is intended to be called from setup.py automatically under Python 3,
and is not intended for end users. The basic idea follows the approach taken
by NumPy with their setup.py file calling tools/py3tool.py to do the 2to3
conversion automatically.

This calls the lib2to3 library functions to convert the Biopython source code
from Python 2 to Python 3, tracking changes to files so that unchanged files
need not be reconverted making development much easier (i.e. if you edit one
source file, doing 'python setup.py install' will only reconvert the one file).
This is done by the last modified date stamps (which will be updated by git if
you switch branches).

NOTE - This is intended to be run under Python 3 (not under Python 2), but
care has been taken to make it run under Python 2 enough to give a clear error
message. In particular, this meant avoiding with statements etc.
"""
import sys
if sys.version_info[0] < 3:
    sys.stderr.write("Please run this under Python 3\n")
    sys.exit(1)

import shutil
import os
import lib2to3.main
from io import StringIO


def run2to3(filenames):
    stderr = sys.stderr
    handle = StringIO()
    try:
        #Want to capture stderr (otherwise too noisy)
        sys.stderr = handle
        while filenames:
            filename = filenames.pop(0)
            print("Converting %s" % filename)
            #TODO - Configurable options per file?
            args = ["--nofix=long", "--no-diffs", "-n", "-w"]
            e = lib2to3.main.main("lib2to3.fixes", args + [filename])
            if e != 0:
                sys.stderr = stderr
                sys.stderr.write(handle.getvalue())
                os.remove(filename) #Don't want a half edited file!
                raise RuntimeError("Error %i from 2to3 on %s" \
                                   % (e, filename))
            #And again for any doctests,
            e = lib2to3.main.main("lib2to3.fixes", args + ["-d", filename])
            if e != 0:
                sys.stderr = stderr
                sys.stderr.write(handle.getvalue())
                os.remove(filename) #Don't want a half edited file!
                raise RuntimeError("Error %i from 2to3 (doctests) on %s" \
                                   % (e, filename))
    except KeyboardInterrupt:
        sys.stderr = stderr
        sys.stderr.write("Interrupted during %s\n" % filename)
        os.remove(filename) #Don't want a half edited file!
        for filename in filenames:
            if os.path.isfile(filename):
                #Don't want uncoverted files left behind:
                os.remove(filename)
        sys.exit(1)
    finally:
        #Restore stderr
        sys.stderr = stderr


def do_update(py2folder, py3folder, verbose=False):
    if not os.path.isdir(py2folder):
        raise ValueError("Python 2 folder %r does not exist" % py2folder)
    if not os.path.isdir(py3folder):
        os.mkdir(py3folder)
    #First remove any files from the 3to2 conversion which no
    #longer existing the Python 2 origin (only expected to happen
    #on a development machine).
    for dirpath, dirnames, filenames in os.walk(py3folder):
        relpath = os.path.relpath(dirpath, py3folder)
        for d in dirnames:
            new = os.path.join(py3folder, relpath, d)
            old = os.path.join(py2folder, relpath, d)
            if not os.path.isdir(old):
                print("Removing %s" % new)
                shutil.rmtree(new)
        for f in filenames:
            new = os.path.join(py3folder, relpath, f)
            old = os.path.join(py2folder, relpath, f)
            if not os.path.isfile(old):
                print("Removing %s" % new)
                os.remove(new)
    #Check all the Python 2 original files have been copied/converted
    #Note we need to do all the conversions *after* copying the files
    #so that 2to3 can detect local imports successfully.
    to_convert = []
    for dirpath, dirnames, filenames in os.walk(py2folder):
        if verbose: print("Processing %s" % dirpath)
        relpath = os.path.relpath(dirpath, py2folder)
        #This is just to give cleaner filenames
        if relpath[:2] == "/.":
            relpath = relpath[2:]
        elif relpath == ".":
            relpath = ""
        for d in dirnames:
            new = os.path.join(py3folder, relpath, d)
            if not os.path.isdir(new):
                os.mkdir(new)
        for f in filenames:
            if f.startswith("."):
                #Ignore hidden files
                continue
            elif f.endswith("~") or f.endswith(".bak") \
            or f.endswith(".swp"):
                #Ignore backup files
                continue
            elif f.endswith(".pyc") or f.endswith("$py.class"):
                #Ignore compiled python
                continue
            old = os.path.join(py2folder, relpath, f)
            new = os.path.join(py3folder, relpath, f)
            #The filesystem can (in Linux) record nanoseconds, but
            #when copying only microsecond accuracy is used.
            #See http://bugs.python.org/issue10148
            #Compare modified times down to milliseconds only. In theory
            #might able to use times down to microseconds (10^-6), but
            #that doesn't work on this Windows machine I'm testing on.
            if os.path.isfile(new) \
            and round(os.stat(new).st_mtime*1000) >= \
                round(os.stat(old).st_mtime*1000):
                if verbose: print("Current: %s" % new)
                continue
            #Python, C code, data files, etc - copy with date stamp etc
            shutil.copy2(old, new)
            assert abs(os.stat(old).st_mtime-os.stat(new).st_mtime)<0.0001, \
                   "Modified time not copied! %0.8f vs %0.8f, diff %f" \
                   % (os.stat(old).st_mtime, os.stat(new).st_mtime,
                      abs(os.stat(old).st_mtime-os.stat(new).st_mtime))
            if f.endswith(".py"):
                #Also run 2to3 on it
                to_convert.append(new)
                if verbose: print("Will convert %s" % new)
            else:
                if verbose: print("Updated %s" % new)
    if to_convert:
        print("Have %i python files to convert" % len(to_convert))
        run2to3(to_convert)

            
def main(python2_source, python3_source,
         children=["Bio", "BioSQL", "Tests", "Scripts", "Doc"]):
    #Note want to use different folders for Python 3.1, 3.2, etc
    #since the 2to3 libraries have changed so the conversion
    #may differ slightly.
    print("The 2to3 library will be called automatically now,")
    print("and the converted files cached under %s" % python3_source)
    if not os.path.isdir(python3_source):
        os.mkdir(python3_source)
    for child in children:
        print("Processing %s" % child)
        do_update(os.path.join(python2_source, child),
                  os.path.join(python3_source, child))
    print("Python 2to3 processing done.")
              
if __name__ == "__main__":
    python2_source = "."
    python3_source = "build/py%i.%i" % sys.version_info[:2]
    main(python2_source, python3_source)
