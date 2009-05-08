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

Tier 1:
    branch_length
    clade
    code
    confidence
    name
    phylogeny
    phyloxml
    taxonomy

Tier 2:
    domain
    domain_architecture
    duplications
    events
    scientific_name
    sequence
    speciations


See Also
--------

Project page on NESCent:
https://www.nescent.org/wg_phyloinformatics/PhyloSoC:Biopython_support_for_parsing_and_writing_phyloXML

Biopython wiki doc:
http://biopython.org/wiki/PhyloXML

PhyloXML homepage:
http://www.phyloxml.org/

GSoC project proposal:
http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798969

This GitHub branch:
http://github.com/etal/biopython/tree/phyloxml


Timeline
--------

:Backlog:
    Documentation on Biopython wiki:

    - explain use cases

:Done:
    - Add sample phyloXML files from phyloxml.org to this repository
    - Add a list of "core" elements to this file
    - Review existing modules for ideas and conventions [*ongoing*]
    - Choose an XML parser [*xml.etree, or lxml/elementtree for Py2.4*]
    - Choose a warning system [*will match Bio.PDB, using warnings module*]
    - Decide repository layout [*code in Bio/PhyloXML/, tests in Tests/*]

    Unit tests:

    - no-op parsing of small example phyloXML files with xml.etree
    - no-op parsing of a large zipped phyloXML file
    - basic loading of the root node into a Python object
    - loading core elements

    Documentation (Biopython wiki):

    - Explain xml.etree, list 3rd-party equivalents for Py2.4


:6/1:
    Parsing from file:

    - Write a simple base class for all phyloXML elements ::

        class PhyloElement(object): ...

    - Write an easy wrapper function for loading local files by name ::

        def read(path): ...

    - Create a specific exception to raise for invalid phyloXML files ::

        class PhyloXMLError(Exception): pass

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

