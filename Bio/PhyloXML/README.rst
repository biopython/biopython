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

PhyloXML homepage:
    http://www.phyloxml.org/

Biopython wiki documentation:
    http://biopython.org/wiki/PhyloXML

Project page on NESCent:
    https://www.nescent.org/wg_phyloinformatics/PhyloSoC:Biopython_support_for_parsing_and_writing_phyloXML

GSoC project proposal:
    http://socghop.appspot.com/student_project/show/google/gsoc2009/nescent/t124022798969

This GitHub branch:
    http://github.com/etal/biopython/tree/phyloxml


Timeline
========

:Done:
    - Add sample phyloXML files from phyloxml.org to this repository
    - Add a list of "core" elements to this file
    - Review existing modules for ideas and conventions [*ongoing*]
    - Choose an XML parser [*xml.etree, or lxml/elementtree for Py2.4*]
    - Choose a warning system [*warnings module, module-specific warning types*]
    - Decide repository layout [*code in Bio/Tree/ and Bio/TreeIO/, tests in
      Tests/*]
    - Test XSD validation for phyloXML v1.10 [*it passes*]

    Code:

    - Simple base class for all phyloXML elements [*Tree.PhyloElement*]
    - Class definitions for all phyloXML types [*PhyloXML.py*]
    - Specific exception and warning to raise for noncompliant phyloXML files
    - Parse individual phylogenies or read a complete PhyloXML object from an
      XML stream [*PhyloXMLIO.py*]
    - Write a complete PhyloXML object to a file or stream [*PhyloXMLIO.py*]
    - XML utility: dump_tags [*Parser.py*]
    - Tree utilies: pretty_print, to_networkx, draw_graphviz [*Utils.py*]
    - Methods to convert PhyloXML objects to core Biopython types:
        - PhyloXML.Sequence to/from SeqRecord (*in progress*)
        - PhyloXML.BranchColor to HTML/CSS-friendly hex string [*to_rgb*]
    - Methods for navigating and using PhyloXML objects:
        - Extract sub-trees as new trees [*Phylogeny.to_phyloxml,
          Clade.to_phylogeny*]
        - Search all sub-nodes for matching attributes
          [*Phylogeny/Clade.find()*]

    - Sugar:
        - \*.__str__(): pretty representation for displaying nodes
            - Customized for Taxonomy and Date
        - \*.__repr__(): simplified object instantiation code
        - Phyloxml.__getitem__(): access trees by index or name
        - Clade.__getitem__(): extended index access to sub-clades
        - Tree.Events acts like a mapping type (dictionary)
        - For plural attributes that are usually single (taxonomies,
          confidences), a singular property retrieves a single item

    Unit tests:

    - no-op parsing of example phyloXML files with xml.etree [*test_dump_tags*]
    - write a complete, readable representation of a phyloXML object
      [*test_pretty_print*]
    - build a complete phyloXML object from each example file
    - parse individual phylogenies as needed
    - create trees of at least 3 clades deep
    - instantiate all phyloXML element types
    - round-trip parsing and serialization of three example files

    Documentation (Biopython wiki):

    - PhyloXML:
        - Explain xml.etree, list 3rd-party equivalents for Py2.4
        - Usage: parsing, writing, object navigation, Bio integration, utilities
        - Performance: read, parse, write

    - TreeIO:
        - stub

    - Tree:
        - stub

:8/3:
    Enhancements (time permitting):

    - Parser class ENH:
        - parse sequence and taxonomy incrementally, too? (handle 'other')
        - drop classmethod decorators
        - drop "to" and "parse" from method names

    - Port common methods to Bio.Tree.BaseTree -- see Bio.Nexus.Tree, Bioperl
      node objects, PyCogent, p4-phylogenetics, lagrange, newick

        - Tree method: update_nested_set_index
            - calculate left_idx, right_idx for nested-set representation
            - see http://www.oreillynet.com/pub/a/network/2002/11/27/bioconf.html

        - 'external' boolean kwarg to Tree.find(): None=all nodes (default),
          True=external nodes only, False=internal nodes only
        - Also: get_terminals, is_identical, distance

    Automated testing:

    - Re-run performance benchmarks
    - Run tests and benchmarks on alternate platforms
    - Check epydoc's generated API documentation

    Update wiki documentation with new features:

    - Tree: find(), base classes
    - TreeIO: 'phyloxml' and 'nexus' wrappers; PhyloXMLIO extras; warn that
      Nexus wrappers don't return Bio.Tree objects yet
    - PhyloXML: singular properties, improved str()

    Discuss merging back upstream.

:8/10:
    Soft "pencils down":

    - Scrub wiki documentation -- PhyloXML, Tree, TreeIO
    - Check unit tests for complete coverage
    - NB: Deadline is Aug. 17


Notes
=====

Core Elements
-------------

See:
    * http://www.phyloxml.org/documentation/version_100/phyloxml.xsd.html
    * http://www.phyloxml.org/examples/phyloxml_examples.xml

Tier 0 (essential tree structure):

    phyloxml, phylogeny, clade

Tier 1 (used in all example files):

    branch_length, confidence, name, taxonomy, code

Tier 2 (used in at least one example file, but not all):

    accession, alt, annotation, bc, binary_characters, clade_relation,
    common_name, date, desc, description, distribution, domain,
    domain_architecture, duplications, events, gained, id, lat, long, lost,
    mol_seq, point, present, property, rank, reference, scientific_name,
    sequence, sequence_relation, speciations, symbol, uri, value

Tier 3 (not found in example files):

    absent, color, red, blue, green, location, losses, polygon, node_id, width

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

