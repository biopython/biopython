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


Links
=====

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
========

:Backlog:
    Documentation on Biopython wiki:

    - explain use cases

    Integration:

    - Play nicely with Nexus, Newick
    - Identify more Biopython objects to reuse or export to

:Done:
    - Add sample phyloXML files from phyloxml.org to this repository
    - Add a list of "core" elements to this file
    - Review existing modules for ideas and conventions [*ongoing*]
    - Choose an XML parser [*xml.etree, or lxml/elementtree for Py2.4*]
    - Choose a warning system [*will match Bio.PDB, using warnings module*]
    - Decide repository layout [*code in Bio/PhyloXML/, tests in Tests/*]

    Code:

    - Wrapper function for loading local files by name [*__init__.read()*]
    - Specific exception and warning to raise for invalid phyloXML files
      [*Exceptions.py*]
    - Simple base class for all phyloXML elements [*Parser.PhyloElement*]
    - Instantiate Tier 0, 1 and 2 elements from an XML stream ("end" events)

    Unit tests:

    - no-op parsing of small example phyloXML files with xml.etree
    - no-op parsing of a large zipped phyloXML file
    - basic loading of the root node into a Python object
    - creating a list of phylogeny objects under the root node
    - creating trees of at least 3 clades deep

    Documentation (Biopython wiki):

    - Explain xml.etree, list 3rd-party equivalents for Py2.4

:6/15:
    Verification and documentation:

    - Finish unittests for parsing and instantiating core elements
    - Test and check parser performance versus Bioperl and Archaeopterix loading
      time
    - Document results of parser testing and performance (on wiki or here)
    - Document basic usage of the parser on the Biopython wiki

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

    - Write unit tests for Pythonic syntax sugar (e.g. __getattr__, __getitem__,
      __contains__)

        - PhyloXML.__getitem__(): get the phylogeny with matching name or id
        - \*.__str__(): pretty representation for printing nodes

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


Notes
=====

Core Elements
-------------

See:
    * http://www.phyloxml.org/documentation/version_100/phyloxml.xsd.html
    * http://www.phyloxml.org/examples/phyloxml_examples.xml

Tier 0 (essential tree structure):

    - done: phyloxml, phylogeny, clade

Tier 1 (used in all example files):

    - done: branch_length, confidence, name, taxonomy, code

Tier 2 (used in at least one example file, but not all):

    - done: accession, alt, annotation, clade_relation, common_name, date, desc,
        description, distribution, domain, domain_architecture, duplications,
        events, id, lat, long, mol_seq, point, property, rank, scientific_name,
        sequence, sequence_relation, speciations, symbol, uri value,

Tier 3 (not found in example files):

    - to do:
        absent, 
        bc, 
        binary_characters,
        gained,
        lost,
        polygon,
        present,
        reference,
        width

    - done:
        blue,
        color,
        green,
        location,
        losses,
        node_id,
        red,
        type,

Namespaces:

    :phy:   http://www.phyloxml.org
    :xml:   http://www.w3.org/XML/1998/namespace
    :xs:    http://www.w3.org/2001/XMLSchema


Diagram
-------

::

    phyloxml
        { xsi:schemaLocation="..."}
        phylogeny * (none)
            { rooted=bool
              rerootable=bool
              branch_length_unit=token
              type=token
            }
            name ? (token)
            id ? (token)
                { type=token }
            description ? (token)
            date ? (token or number)
                { unit=token
                  range=
                }
                desc ? (token)
                value ? (token?)
            confidence * (double)
                { type=token }
            clade ? (none)
                { branch_length=number
                  id_source=identifier
                }
                name ^
                branch_length ?     # same as using the attribute
                confidence ^
                width ?
                color ?
                    red (byte)
                    green (byte)
                    blue (byte)
                node_id ?           # see id
                taxonomy *
                    { type=
                      id_source=
                    }
                    id ^
                    code ? ( [a-zA-Z0-9_]{2,10} )   # see TaxonomyCode
                    scientific_name ? (token)
                    common_name * (token)
                    rank ? (one of:
                        ['domain', 'kingdom', 'subkingdom', 'branch',
                        'infrakingdom', 'superphylum', 'phylum', 'subphylum',
                        'infraphylum', 'microphylum', 'superdivision',
                        'division', 'subdivision', 'infradivision',
                        'superclass', 'class', 'subclass', 'infraclass',
                        'superlegion', 'legion', 'sublegion', 'infralegion',
                        'supercohort', 'cohort', 'subcohort', 'infracohort',
                        'superorder', 'order', 'suborder', 'superfamily',
                        'family', 'subfamily', 'supertribe', 'tribe',
                        'subtribe', 'infratribe', 'genus', 'subgenus',
                        'superspecies', 'species', 'subspecies', 'variety',
                        'subvariety', 'form', 'subform', 'cultivar', 'unknown',
                        'other'] )
                    uri ? (token, generally URL)
                        { desc=token
                          type=token
                        }
                    OTHER *
                sequence *
                    { type=token
                      id_source=token
                      id_ref=identifier
                    }
                    symbol ? ( \S{1,10} )
                    accession ? (token)
                        { source=token }
                    name ^
                    location ?
                    mol_seq ? ( [a-zA-Z\.\-\?\*_]+ )
                    uri ^
                    annotation +
                        { ref=[a-zA-Z0-9_]+:[a-zA-Z0-9_\.\-\s]+
                          source=token
                          evidence=
                          type=
                        }
                        desc ^
                        confidence ^
                        property * (none)
                            { ref=^
                              unit=a-zA-Z0-9_]+:[a-zA-Z0-9_\.\-\s]+
                              datatype=
                                ['xsd:string', 'xsd:boolean', 'xsd:decimal',
                                'xsd:float', 'xsd:double', 'xsd:duration',
                                'xsd:dateTime', 'xsd:time', 'xsd:date',
                                'xsd:gYearMonth', 'xsd:gYear', 'xsd:gMonthDay',
                                'xsd:gDay', 'xsd:gMonth', 'xsd:hexBinary',
                                'xsd:base64Binary', 'xsd:anyURI',
                                'xsd:normalizedString', 'xsd:token',
                                'xsd:integer', 'xsd:nonPositiveInteger',
                                'xsd:negativeInteger', 'xsd:long', 'xsd:int',
                                'xsd:short', 'xsd:byte',
                                'xsd:nonNegativeInteger', 'xsd:unsignedLong',
                                'xsd:unsignedInt', 'xsd:unsignedShort',
                                'xsd:unsignedByte', 'xsd:positiveInteger']
                              applies_to=
                                ['phylogeny', 'clade', 'node', 'annotation',
                                'parent_branch', 'other']
                              id_ref=identifier
                            }
                        uri ^
                    domain_architecture ?
                        { length=int }
                        domain + (token)
                            { from=int >0
                              to=int >0
                              confidence=double
                              id=token
                            }
                    OTHER *
                events ?
                    type ? (one of:
                        ['transfer', 'fusion', 'speciation_or_duplication',
                        'other', 'mixed', 'unassigned'] )   # see EventType
                    duplications ?
                    speciations ?
                    losses ?
                    confidence ^
                binary_characters ? (none)
                    { type=
                      gained_count=
                      lost_count=
                      present_count=
                      absent_count=
                    }
                    gained ?            # see BinaryCharacterList for these
                        bc + (token)
                    lost ?
                    present ?
                    absent ?
                distribution * (none)
                    desc ^
                    point * (none)
                        { geodetic_datum="WGS84" }
                        lat (double)
                        long (double)
                        alt ? (int?)
                    polygon * (none)    # list of at least 3 points
                        point {3,} ^
                date ^
                reference * (none)
                    { doi=[a-zA-Z0-9_\.]+/[a-zA-Z0-9_\.]+ }
                    desc ^
                property ^
                clade ^
            clade_relation * (none)
                { id_ref_0=identifier
                  id_ref_1=identifier
                  distance=
                  type=token
                }
                confidence ^
            sequence_relation *
                { id_ref_0=identifier
                  id_ref_1=identifier
                  distance=
                  type=SequenceRelationType
                }
                confidence ^
            property ^
            OTHER *
        NOT *                       # arbitrary elements from other namespaces

