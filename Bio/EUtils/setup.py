#!/usr/bin/env python

import sys
from distutils.core import setup
try:
    import EUtils
except ImportError:
    import __init__ as EUtils

def _dict(**kwargs):
    return kwargs

d = _dict(
    name = "EUtils",
    version = EUtils.__version__,
    description = "Client interface to NCBI's EUtils/Entrez server",
    author = "Andrew Dalke",
    author_email = "dalke@dalkescientific.com",
    maintainer = "Dalke Scientific Software, LLC",
    maintainer_email = "dalke@dalkescientific.com",
      
    url = "http://www.dalkescientific.com/EUtils/",

    long_description = """\
EUtils is a client library for the Entrez databases at NCBI.

NCBI provides the EUtils web service so that software can query Entrez
directly, rather than going through the web interface and dealing with
the hassles of web scraping.  For more information see

  http://www.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html

This package provides two levels of interface.  The lowest one makes a
programmatic interface to construct the query URL and make the
request.  The higher level ones support history tracking and parsing
of query results.  These greatly simplify working with the EUtils
server.
""",

    package_dir = {"": ".."},     
    packages = ["EUtils", "EUtils.DTDs"],
      
    classifiers = [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: Freely Distributable",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Bio-Informatics", # a '-'? !
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Internet",
        ],
    )
if sys.version_info < (2,2,4):
    del d["classifiers"]


if __name__ == "__main__":
    setup(**d)
