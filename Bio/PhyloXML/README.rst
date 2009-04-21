+++++++++++++++++++++++++++++
PhyloXML module for Biopython
+++++++++++++++++++++++++++++

PhyloXML is an XML format for phylogenetic trees, designed to allow storing
information about the trees themselves (such as branch lengths and multiple
support values) along with data such as taxonomic and genomic annotations.
Connecting these pieces of evolutionary information in a standard format is key
for comparative genomics.

A Bioperl driver for phyloXML was created during the 2008 Summer of Code; this
project aims to build a similar module for the popular Biopython package.


Core Elements
-------------

[todo]


See Also
--------

PhyloXML homepage: http://www.phyloxml.org/

Biopython wiki: http://biopython.org/wiki/Active_projects

This GitHub branch: http://github.com/etal/biopython/tree/phyloxml

GSoC project page: http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798969

BioPerl implementation (2008): http://www.bioperl.org/wiki/Phyloxml_Project_Demo

BioRuby implementation (current): http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022800294



Timeline
--------

:4/20:
    Getting started:

    - Add sample phyloXML files from phyloxml.org to this repository
    - Add a list of "core" elements to this file
    - Start a PhyloXML page on the Biopython wiki

        - http://biopython.org/wiki/Active_projects

:4/27:
    Community bonding:

    - Review existing modules for ideas and conventions
    - Discuss some topics on biopython-dev:

        - repository layout
        - xml.sax vs. xml.etree and Python 2.4 compatibility timeline
        - warning/logging/strictness system for parsing odd or malformed files

:5/4:
    - Hash out these topics on biopython-dev
    - Document decisions here and on the wiki page
    - Start writing unit tests for basic XML parsing

:5/11:
    Coding should be underway now.

    - Write unit tests that fail:

        - no-op parsing of phyloXML files with xml.sax or xml.etree
        - basic loading of the root node into a Python object
        - loading core elements

:6/1:
    Parsing from file:

    - Write a simple base class for all phyloXML elements ::

        class Tree(object): ...

    - Write an easy wrapper function for loading local files by name ::

        def load(path): ...

    - Create a specific exception to raise for invalid phyloXML files ::

        class ParseException(Exception): pass

:6/8:
    Map nodes to classes:

    - Write classes corresponding to "core" XML elements

    - Instantiate from XML nodes (sax events or etree elements)

:6/15:
    Verification and documentation:

    - Finish unittests for parsing and instantiating core elements
    - Test and check parser performance versus Bioperl and Archaeopterix loading
      time
    - Document results of parser testing and performance (on wiki or here)

:6/22:
    Serialization back to file:

    - Write unit tests for serialization
    - Write serialization methods for each class
    - Write a top-level function for triggering serialization of the whole
      hierarchy

:6/29:
    Verification and documentation:

    - Test and benchmark serialization of object hierarchy
    - Document results of serialization testing and benchmarking

:7/6:
    Make it pretty:

    - Write unit tests for Pythonic syntax sugar (e.g.  __getitem__)
    - Add the corresponding magic methods to the base class

:7/13:
    Extend the core to the rest of the spec:

    - Begin adding unit tests and classes to support additional (non-core)
      phyloXML elements

:7/20:
    - Continue adding support for the rest of the phyloXML spec and testing
    - Hopefully, the entire phyloXML spec is covered by the end of this week

:7/27:
    Document all completed functionality.

:8/3:
    - Run tests and benchmarks on alternate platforms and document results
    - Discuss merging back upstream

