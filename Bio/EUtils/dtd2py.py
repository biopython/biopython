#!/usr/bin/env python
"""
Internal class to convert DTDs to python form.

dtd2py [-m mixinmodule] [-u] <DTDfile>...

Converts the given DTD to a python form, and stores it in the dtds package
directory. That directory should be writable by the user running this. 
The new dtds module will have the same base name as the DTD, but with dashes and dots converted to underscores.

Options:
	-m <modulename> -- specify a mixin module to load. If a class in this
	module matches the DTD element name plus "Mixin" appended, then that class
	will be inherited from in the generated element class.
	-u  -- convert all element tag names to uppercase. Otherwise, case is
	preserved. Some XML applications specify using uppercase, others case
	preserving. Some XML element names may match python keyword or builtin
	names, and those DTDs should also be converted to uppercase to avoid syntax
	errors in the generated file.

"""

import sys, os
import getopt

from Bio.EUtils import POM

def compile_file(dtdfile, mixin, toupper):
	outfilename = POM.get_mod_file(dtdfile)
##	try:
	if 1:
		outfile = open(outfilename, "w")
		try:
			dtdp = POM.get_dtd_compiler(outfile, mixin, toupper)
			dtdp.parse_resource(dtdfile)
		finally:
			outfile.close()
##	except:
##		os.unlink(outfilename)
##		print "COMPILE ERROR: while parsing:", dtdfile
##		import traceback, pdb
##		t, v, tb = sys.exc_info()
##		traceback.print_exception(t, v, tb)
#		pdb.pm()

def do_it(filelist, mixin, toupper):
	if len(filelist) == 1:
		compile_file(filelist[0], mixin, toupper)
	else:
		for dtdfile in filelist:
			print "Compiling:", dtdfile
			compile_file(dtdfile, mixin, toupper)


def main(argv):
	if len(argv) > 1:
		mixin = None
		toupper = 0
		args, files = getopt.getopt(argv[1:], "m:u")
		for opt, value in args:
			if opt == "-m":
				try:
					mixin = __import__(value)
				except:
					print "Warning: could not import mixin module"
					mixin = None
			if opt == "-u":
				toupper = 1
		do_it(files, mixin, toupper)
	else:
		print __doc__


if __name__ == "__main__":
	main(sys.argv)

