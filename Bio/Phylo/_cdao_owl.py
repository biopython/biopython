# Copyright (C) 2013 by Ben Morris (ben@bendmorris.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

import xml.etree.ElementTree as ET

cdao_namespaces = {
    'cdao': 'http://purl.obolibrary.org/obo/cdao.owl#',
    'obo': 'http://purl.obolibrary.org/obo/',
}


def resolve_uri(s, namespaces=cdao_namespaces, cdao_to_obo=True, xml_style=False):
    """Converts prefixed URIs to full URIs.

    Optionally, converts CDAO named identifiers to OBO numeric identifiers.
    """

    if cdao_to_obo and s.startswith('cdao:'):
        return resolve_uri('obo:%s' % cdao_elements[s[5:]], namespaces, cdao_to_obo)

    for prefix in namespaces:
        if xml_style:
            s = s.replace(prefix + ':', '{%s}' % namespaces[prefix])
        else:
            s = s.replace(prefix + ':', namespaces[prefix])

    return s


cdao_owl = '''<?xml version="1.0"?>


<!DOCTYPE rdf:RDF [
    <!ENTITY owl "http://www.w3.org/2002/07/owl#" >
    <!ENTITY obo "http://purl.obolibrary.org/obo/" >
    <!ENTITY dc "http://purl.org/dc/elements/1.1/" >
    <!ENTITY xsd "http://www.w3.org/2001/XMLSchema#" >
    <!ENTITY rdfs "http://www.w3.org/2000/01/rdf-schema#" >
    <!ENTITY rdf "http://www.w3.org/1999/02/22-rdf-syntax-ns#" >
]>


<rdf:RDF xmlns="http://www.evolutionaryontology.org/cdao/1.0/cdao.owl#"
     xml:base="http://www.evolutionaryontology.org/cdao/1.0/cdao.owl"
     xmlns:dc="http://purl.org/dc/elements/1.1/"
     xmlns:obo="http://purl.obolibrary.org/obo/"
     xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
     xmlns:owl="http://www.w3.org/2002/07/owl#"
     xmlns:xsd="http://www.w3.org/2001/XMLSchema#"
     xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">
    <owl:Ontology rdf:about="&obo;cdao.owl">
        <dc:coverage rdf:datatype="&xsd;string">Comparison of two or more biological entities of the same class when the similarities and differences of the entities are treated explicitly as the product of an evolutionary process of descent with modification.</dc:coverage>
        <dc:description rdf:datatype="&xsd;string">The Comparative Data Analysis Ontology (CDAO) provides a framework for understanding data in the context of evolutionary-comparative analysis.  This comparative approach is used commonly in bioinformatics and other areas of biology to draw inferences from a comparison of differently evolved versions of something, such as differently evolved versions of a protein.  In this kind of analysis, the things-to-be-compared typically are classes called &#39;OTUs&#39; (Operational Taxonomic Units).  The OTUs can represent biological species, but also may be drawn from higher or lower in a biological hierarchy, anywhere from molecules to communities.  The features to be compared among OTUs are rendered in an entity-attribute-value model sometimes referred to as the &#39;character-state data model&#39;.  For a given character, such as &#39;beak length&#39;, each OTU has a state, such as &#39;short&#39; or &#39;long&#39;.  The differences between states are understood to emerge by a historical process of evolutionary transitions in state, represented by a model (or rules) of transitions along with a phylogenetic tree.  CDAO provides the framework for representing OTUs, trees, transformations, and characters.  The representation of characters and transformations may depend on imported ontologies for a specific type of character.</dc:description>
        <dc:creator xml:lang="en">CDAO Team</dc:creator>
        <dc:title xml:lang="en">Comparative Data Analysis Ontology</dc:title>
        <dc:subject xml:lang="en">comparative analysis; comparative data analysis; evolutionary comparative analysis; evolution;  phylogeny; phylogenetics</dc:subject>
        <dc:rights rdf:resource="http://creativecommons.org/publicdomain/zero/1.0/"/>
        <owl:versionIRI rdf:resource="&obo;cdao/2012-06-06/cdao.owl"/>
        <owl:imports rdf:resource="&obo;iao/ontology-metadata.owl"/>
    </owl:Ontology>
    


    <!--
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // Annotation properties
    //
    ///////////////////////////////////////////////////////////////////////////////////////
     -->

    <owl:AnnotationProperty rdf:about="&dc;creator"/>
    <owl:AnnotationProperty rdf:about="&dc;subject"/>
    <owl:AnnotationProperty rdf:about="&dc;description"/>
    <owl:AnnotationProperty rdf:about="&dc;coverage"/>
    <owl:AnnotationProperty rdf:about="&dc;language"/>
    <owl:AnnotationProperty rdf:about="&dc;identifier"/>
    <owl:AnnotationProperty rdf:about="&dc;date"/>
    <owl:AnnotationProperty rdf:about="&dc;source"/>
    <owl:AnnotationProperty rdf:about="&dc;title"/>
    <owl:AnnotationProperty rdf:about="&dc;rights"/>
    


    


    <!--
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // Object Properties
    //
    ///////////////////////////////////////////////////////////////////////////////////////
     -->

    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000142 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000142">
        <rdfs:label rdf:datatype="&xsd;string">has_Character</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property associates a character data matrix with a character (a column) represented in the matrix.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000056"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000071"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000143 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000143">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Edge_as_Child</rdfs:label>
        <dc:description>The property links a Node to the Edge it belongs to in the child position.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000139"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000146"/>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000209"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000144 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000144">
        <rdf:type rdf:resource="&owl;TransitiveProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">has_Ancestor</rdfs:label>
        <dc:description>The property links a node to any of the other nodes that are its ancestors in a rooted tree.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000174"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000145 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000145">
        <rdfs:label rdf:datatype="&xsd;string">has_Nucleotide_State</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property associates a nucleotide character-state instance with a state value from the domain of nucleotide states.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000002"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000184"/>
        <rdfs:range>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000015"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000133"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:range>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000146 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000146">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Edge</rdfs:label>
        <dc:description>The property links a Node to one of the edges that are incident on such node.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000099"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000190"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000147 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000147">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Character_State_Data_Matrix</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000056"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000190"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000071"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000148 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000148">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">has_Root</rdfs:label>
        <dc:description>The property links a rooted tree to the specific node that represents the unique root of the tree.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000012"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000149 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000149">
        <rdfs:label rdf:datatype="&xsd;string">has_Child</rdfs:label>
        <dc:description>The property links a node to a node that is an immediate descendant in the tree.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000174"/>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000179"/>
        <owl:propertyChainAxiom rdf:parseType="Collection">
            <rdf:Description rdf:about="&obo;CDAO_0000177"/>
            <rdf:Description rdf:about="&obo;CDAO_0000209"/>
        </owl:propertyChainAxiom>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000150 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000150">
        <rdfs:label rdf:datatype="&xsd;string">has_First_Coordinate_Item</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">The property that relates a coordinate list to the first item in the list.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000092"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
        <rdfs:range>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000003"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000095"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:range>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000151 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000151">
        <rdfs:label rdf:datatype="&xsd;string">has_Coordinate</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000022"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000098"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000152 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000152">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Continuous_Character</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000019"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000068"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000205"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000153 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000153">
        <rdfs:label rdf:datatype="&xsd;string">has_Datum</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a character to a state datum for the character.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000098"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000071"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000154 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000154">
        <rdfs:label rdf:datatype="&xsd;string">has_Standard_Datum</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000008"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000183"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000075"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000155 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000155">
        <rdfs:label rdf:datatype="&xsd;string">subtree_of</rdfs:label>
        <dc:description>This property links two networks where the latter is a substructure of the former</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000006"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000006"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000156 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000156">
        <rdfs:label rdf:datatype="&xsd;string">has_Amino_Acid_State</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property associates a amino acid character-state instance with a state value from the domain of amino acid states.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000112"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000184"/>
        <rdfs:range>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000015"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000076"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:range>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000157 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000157">
        <rdfs:label rdf:datatype="&xsd;string">is_annotation_of</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000040"/>
        <rdfs:range rdf:resource="&owl;Thing"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000158 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000158">
        <rdfs:label rdf:datatype="&xsd;string">has_RNA_Datum</rdfs:label>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000206"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000159 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000159">
        <rdfs:label rdf:datatype="&xsd;string">has_Left_State</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a transformation to a &#39;left&#39; state (the state associated with the &#39;left&#39; node).</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000091"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000097"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000182"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000160 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000160">
        <rdfs:label rdf:datatype="&xsd;string">precedes</rdfs:label>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000161 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000161">
        <rdfs:label rdf:datatype="&xsd;string">exclude</rdfs:label>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000162 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000162">
        <rdfs:label rdf:datatype="&xsd;string">has_Node</rdfs:label>
        <dc:description>Property that associates to each Edge the Nodes it connects.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000099"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000163 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000163">
        <rdfs:label rdf:datatype="&xsd;string">nca_node_of</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000059"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000164 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000164">
        <rdfs:label rdf:datatype="&xsd;string">has_External_Reference</rdfs:label>
        <rdfs:comment rdf:datatype="&rdfs;Literal">Associates a TU to some external taxonomy reference.</rdfs:comment>
        <rdfs:domain rdf:resource="&obo;CDAO_0000138"/>
        <rdfs:range rdf:resource="&owl;Thing"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000165 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000165">
        <rdfs:label rdf:datatype="&xsd;string">has_Coordinate_System</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property links a coordinate to the coordinate system it references.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000022"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000104"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000166 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000166">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Nucleotide_Character</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000002"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000094"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000205"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000167 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000167">
        <rdf:type rdf:resource="&owl;SymmetricProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">connects_to</rdfs:label>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000167"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000168 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000168">
        <rdfs:label rdf:datatype="&xsd;string">has_Amino_Acid_Datum</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates an amino acid character (a column in a protein sequence alignment) to a state datum for the character (an individual cell in the alignment column).</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000112"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000206"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000131"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000169 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000169">
        <rdfs:label rdf:datatype="&xsd;string">hereditary_change_of</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a type of evolutionary change (an Edge_Transformation) to the character that undergoes the change.  The change is a transformation_of the affected character.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000071"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000170 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000170">
        <rdfs:label rdf:datatype="&xsd;string">has_Compound_Datum</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a compound character (a character with some states that are subdividable) to a state datum for the character.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000136"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000183"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000078"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000171 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000171">
        <rdfs:label rdf:datatype="&xsd;string">has_Descendants</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000059"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000080"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000172 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000172">
        <rdfs:label rdf:datatype="&xsd;string">reconciliation_of</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000030"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000110"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000173 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000173">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Amino_Acid_Character</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000112"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000131"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000205"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000174 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000174">
        <rdf:type rdf:resource="&owl;TransitiveProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">has_Descendant</rdfs:label>
        <dc:description>A property that links a node to any of its descendants in a rooted tree.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000175 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000175">
        <rdfs:label rdf:datatype="&xsd;string">has_Continuous_State</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property associates a character-state instance with a state value on a continuous numeric scale.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000019"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000031"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000184"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000176 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000176">
        <rdfs:label rdf:datatype="&xsd;string">has_Type</rdfs:label>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000177 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000177">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Edge_as_Parent</rdfs:label>
        <dc:description>The property links a Node to one of the Edges where the node appears in the parent position (i.e., closer to the root).</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000139"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000146"/>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000201"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000178 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000178">
        <rdfs:label rdf:datatype="&xsd;string">has</rdfs:label>
        <dc:description>Generic &#39;has&#39; property.</dc:description>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000179 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000179">
        <rdfs:label rdf:datatype="&xsd;string">has_Parent</rdfs:label>
        <dc:description>The property that links a node to its unique parent in a rooted tree.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000144"/>
        <owl:propertyChainAxiom rdf:parseType="Collection">
            <rdf:Description rdf:about="&obo;CDAO_0000143"/>
            <rdf:Description rdf:about="&obo;CDAO_0000201"/>
        </owl:propertyChainAxiom>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000180 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000180">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Compound_Character</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000078"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000136"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000205"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000181 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000181">
        <rdfs:label rdf:datatype="&xsd;string">homologous_to</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This propery relates different instances of the same character, including the case when the states of the character differ (e.g., large_beak of beak_size_character of TU A is homologous_to small_beak of beak_size_character of TU B).</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000098"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000098"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000182 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000182">
        <rdfs:label rdf:datatype="&xsd;string">has_Change_Component</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a transformation to the components that compose it.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000097"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000183 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000183">
        <rdfs:label rdf:datatype="&xsd;string">has_Categorical_Datum</rdfs:label>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000153"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000184 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000184">
        <rdfs:label rdf:datatype="&xsd;string">has_State</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property associates a character-state instance with its state value, e.g., a state value expressed in terms of an imported domain ontology.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000091"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000098"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000185 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000185">
        <rdfs:label rdf:datatype="&xsd;string">has_Left_Node</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a transformation to a &#39;left&#39; node (the node that has the &#39;left&#39; state).</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000097"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000182"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000186 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000186">
        <rdfs:label rdf:datatype="&xsd;string">has_Right_State</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a transformation to a &#39;right&#39; state (the state associated with the &#39;right&#39; node).</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000091"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000097"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000182"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000187 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000187">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">represents_TU</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a TU or taxonomic unit (typically associated with character data) to a phylogenetic history (Tree).</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000138"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000188 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000188">
        <rdfs:label rdf:datatype="&xsd;string">exclude_Node</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000161"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000189 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000189">
        <rdfs:label rdf:datatype="&xsd;string">has_Compound_State</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property associates a compound character-state instance with its compound state value.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000136"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000184"/>
        <rdfs:range>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000015"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000055"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:range>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000190 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000190">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to</rdfs:label>
        <dc:description>Generic property that links a concept to another concept it is a constituent of. The property is a synonym of part_of.</dc:description>
        <owl:equivalentProperty rdf:resource="&obo;CDAO_0000194"/>
        <rdfs:range rdf:resource="&owl;Thing"/>
        <rdfs:domain rdf:resource="&owl;Thing"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000191 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000191">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_TU</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a character-state datum to its TU.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000098"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000138"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000190"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000192 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000192">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Network</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000006"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000190"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000099"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000140"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000193 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000193">
        <rdfs:label rdf:datatype="&xsd;string">has_Annotation</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000040"/>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000157"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
        <rdfs:domain rdf:resource="&owl;Thing"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000194 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000194">
        <rdfs:label rdf:datatype="&xsd;string">part_of</rdfs:label>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000195 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000195">
        <rdfs:label rdf:datatype="&xsd;string">has_Nucleotide_Datum</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a nucleotide character (a column in a nucleotide alignment) to a state datum for the character (an individual cell in the alignment column).</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000002"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000206"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000094"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000196 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000196">
        <rdfs:label rdf:datatype="&xsd;string">represented_by_Node</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a TU to a node that represents it in a network.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000138"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000187"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000197 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000197">
        <rdfs:label rdf:datatype="&xsd;string">has_Remaining_Coordinate_List</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">The property that relates a coordinate list to the item in the list beyond the first item.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000092"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000092"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000198 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000198">
        <rdfs:label rdf:datatype="&xsd;string">has_Element</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000118"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000199 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000199">
        <rdfs:label rdf:datatype="&xsd;string">exclude_Subtree</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000070"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000161"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000200 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000200">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Tree</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000110"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000190"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000099"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000140"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000201 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000201">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">has_Parent_Node</rdfs:label>
        <dc:description>Associates to a Directed Edge the Node that is in the parent position in the edge (i.e., the node touched by the edge and closer to the root of the tree)</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000139"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000162"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000202 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000202">
        <rdfs:label rdf:datatype="&xsd;string">has_Lineage_node</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000004"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000203 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000203">
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Tree_as_Root</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000110"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000140"/>
        <owl:inverseOf rdf:resource="&obo;CDAO_0000148"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000190"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000204 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000204">
        <rdfs:label rdf:datatype="&xsd;string">has_Hereditary_Change</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000097"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000099"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000205 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000205">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">belongs_to_Character</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000071"/>
        <rdfs:domain rdf:resource="&obo;CDAO_0000098"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000190"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000206 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000206">
        <rdfs:label rdf:datatype="&xsd;string">has_Molecular_Datum</rdfs:label>
        <rdfs:range rdf:resource="&obo;CDAO_0000050"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000183"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000115"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000207 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000207">
        <rdfs:label rdf:datatype="&xsd;string">has_Continuous_Datum</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a continuous character to a state datum for the character.</dc:description>
        <rdfs:range rdf:resource="&obo;CDAO_0000019"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000153"/>
        <rdfs:domain>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000068"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000138"/>
                </owl:unionOf>
            </owl:Class>
        </rdfs:domain>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000208 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000208">
        <rdfs:label rdf:datatype="&xsd;string">has_TU</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property associates a character data matrix with a TU (a row) represented in the matrix.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000056"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000138"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000178"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000209 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000209">
        <rdf:type rdf:resource="&owl;FunctionalProperty"/>
        <rdfs:label rdf:datatype="&xsd;string">has_Child_Node</rdfs:label>
        <dc:description>The property associates to a Directed Edge the Node that is in the child position in the edge, i.e., the node touched by the edge and closer to the leaves of the tree.</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000139"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000162"/>
    </owl:ObjectProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000210 -->

    <owl:ObjectProperty rdf:about="&obo;CDAO_0000210">
        <rdfs:label rdf:datatype="&xsd;string">has_Right_Node</rdfs:label>
        <dc:description rdf:datatype="&xsd;string">This property relates a transformation to a &#39;right&#39; node (the node that has the &#39;right&#39; state).</dc:description>
        <rdfs:domain rdf:resource="&obo;CDAO_0000097"/>
        <rdfs:range rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000182"/>
    </owl:ObjectProperty>
    


    <!--
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // Data properties
    //
    ///////////////////////////////////////////////////////////////////////////////////////
     -->

    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000211 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000211">
        <rdfs:label rdf:datatype="&xsd;string">has_Precision</rdfs:label>
        <rdfs:range rdf:resource="&xsd;float"/>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000212 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000212">
        <rdfs:label rdf:datatype="&xsd;string">has_Point_Coordinate_Value</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000003"/>
        <rdfs:range rdf:resource="&xsd;integer"/>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000213 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000213">
        <rdfs:label rdf:datatype="&xsd;string">has_Int_Value</rdfs:label>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000215"/>
        <rdfs:range rdf:resource="&xsd;int"/>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000214 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000214">
        <rdfs:label rdf:datatype="&xsd;string">has_Support_Value</rdfs:label>
        <rdfs:range rdf:resource="&xsd;float"/>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000215 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000215">
        <rdfs:label rdf:datatype="&xsd;string">has_Value</rdfs:label>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000216 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000216">
        <rdfs:label rdf:datatype="&xsd;string">has_Uncertainty_Factor</rdfs:label>
        <rdfs:range rdf:resource="&xsd;float"/>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000217 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000217">
        <rdfs:label rdf:datatype="&xsd;string">has_Range_End_Value</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000095"/>
        <rdfs:range rdf:resource="&xsd;integer"/>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000218 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000218">
        <rdfs:label rdf:datatype="&xsd;string">has_Float_Value</rdfs:label>
        <rdfs:subPropertyOf rdf:resource="&obo;CDAO_0000215"/>
        <rdfs:range rdf:resource="&xsd;float"/>
    </owl:DatatypeProperty>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000219 -->

    <owl:DatatypeProperty rdf:about="&obo;CDAO_0000219">
        <rdfs:label rdf:datatype="&xsd;string">has_Range_Start_Value</rdfs:label>
        <rdfs:domain rdf:resource="&obo;CDAO_0000095"/>
        <rdfs:range rdf:resource="&xsd;integer"/>
    </owl:DatatypeProperty>
    


    <!--
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // Classes
    //
    ///////////////////////////////////////////////////////////////////////////////////////
     -->

    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000002 -->

    <owl:Class rdf:about="&obo;CDAO_0000002">
        <rdfs:label rdf:datatype="&xsd;string">DesoxiRibonucleotideResidueStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000050"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000003 -->

    <owl:Class rdf:about="&obo;CDAO_0000003">
        <rdfs:label rdf:datatype="&xsd;string">CoordinatePoint</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000022"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000004 -->

    <owl:Class rdf:about="&obo;CDAO_0000004">
        <rdfs:label rdf:datatype="&xsd;string">Lineage</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000012"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000202"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000140"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000005 -->

    <owl:Class rdf:about="&obo;CDAO_0000005">
        <rdfs:label rdf:datatype="&xsd;string">Phylo4Tree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000074"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000006 -->

    <owl:Class rdf:about="&obo;CDAO_0000006">
        <rdfs:label rdf:datatype="&xsd;string">Network</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000178"/>
                <owl:allValuesFrom>
                    <owl:Class>
                        <owl:unionOf rdf:parseType="Collection">
                            <rdf:Description rdf:about="&obo;CDAO_0000099"/>
                            <rdf:Description rdf:about="&obo;CDAO_0000140"/>
                        </owl:unionOf>
                    </owl:Class>
                </owl:allValuesFrom>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000007 -->

    <owl:Class rdf:about="&obo;CDAO_0000007">
        <rdfs:label rdf:datatype="&xsd;string">ModelDescription</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000040"/>
        <dc:description>Description of a model of transformations.</dc:description>
        <rdfs:comment>This is a non-computible description of a model, not the fully specified mathematical model, which typically relates the probability of a transformation to various parameters.</rdfs:comment>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000008 -->

    <owl:Class rdf:about="&obo;CDAO_0000008">
        <rdfs:label rdf:datatype="&xsd;string">StandardStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000089"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000009 -->

    <owl:Class rdf:about="&obo;CDAO_0000009">
        <rdfs:label rdf:datatype="&xsd;string">ContinuousCharacterLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000063"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000010 -->

    <owl:Class rdf:about="&obo;CDAO_0000010">
        <rdfs:label rdf:datatype="&xsd;string">ContinuousCharBayesianLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000009"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000011 -->

    <owl:Class rdf:about="&obo;CDAO_0000011">
        <rdfs:label rdf:datatype="&xsd;string">NEXUSTreeBlock</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000074"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000012 -->

    <owl:Class rdf:about="&obo;CDAO_0000012">
        <rdfs:label rdf:datatype="&xsd;string">RootedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000155"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000012"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000148"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000140"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000178"/>
                <owl:allValuesFrom>
                    <owl:Class>
                        <owl:unionOf rdf:parseType="Collection">
                            <rdf:Description rdf:about="&obo;CDAO_0000139"/>
                            <rdf:Description rdf:about="&obo;CDAO_0000140"/>
                        </owl:unionOf>
                    </owl:Class>
                </owl:allValuesFrom>
            </owl:Restriction>
        </rdfs:subClassOf>
        <owl:disjointWith rdf:resource="&obo;CDAO_0000088"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000013 -->

    <owl:Class rdf:about="&obo;CDAO_0000013">
        <rdfs:label rdf:datatype="&xsd;string">Kimura2Parameters</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000020"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000014 -->

    <owl:Class rdf:about="&obo;CDAO_0000014">
        <rdfs:label rdf:datatype="&xsd;string">TreeProcedure</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000044"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000015 -->

    <owl:Class rdf:about="&obo;CDAO_0000015">
        <rdfs:label rdf:datatype="&xsd;string">Generic_State</rdfs:label>
        <owl:equivalentClass>
            <owl:Class>
                <owl:oneOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000222"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000221"/>
                    <rdf:Description rdf:about="&obo;CDAO_0000223"/>
                </owl:oneOf>
            </owl:Class>
        </owl:equivalentClass>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000039"/>
        <rdfs:comment>This class should be renamed.  These are not generic states but non-concrete states including gap, unknown and missing.</rdfs:comment>
        <dc:description>This concept is tied to the verbally ambiguous &#39;gap&#39; concept and to the use of a gap character (often the en dash &#39;-&#39;) in text representations of sequence alignments. In general, this represents the absence of any positively diagnosed Character-State. As such, the gap may be interpreted as an additional Character-State, as the absence of the Character, or as an unknown value.  In some cases it is helpful to separate these.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000016 -->

    <owl:Class rdf:about="&obo;CDAO_0000016">
        <rdfs:label rdf:datatype="&xsd;string">UnrootedSubtree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000070"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000017 -->

    <owl:Class rdf:about="&obo;CDAO_0000017">
        <rdfs:label rdf:datatype="&xsd;string">UnresolvedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000018 -->

    <owl:Class rdf:about="&obo;CDAO_0000018">
        <rdfs:label rdf:datatype="&xsd;string">BifurcatingTree</rdfs:label>
        <owl:equivalentClass rdf:resource="&obo;CDAO_0000130"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000019 -->

    <owl:Class rdf:about="&obo;CDAO_0000019">
        <rdfs:label rdf:datatype="&xsd;string">ContinuousStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000098"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000020 -->

    <owl:Class rdf:about="&obo;CDAO_0000020">
        <rdfs:label rdf:datatype="&xsd;string">SubstitutionModel</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000007"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000021 -->

    <owl:Class rdf:about="&obo;CDAO_0000021">
        <rdfs:label rdf:datatype="&xsd;string">JukesKantor</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000020"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000022 -->

    <owl:Class rdf:about="&obo;CDAO_0000022">
        <rdfs:label rdf:datatype="&xsd;string">DatumCoordinate</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000165"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000104"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>A positional coordinate giving the source of a character state, used for molecular sequences.</dc:description>
        <rdfs:comment>drawing from seqloc categories from NCBI at http://www.ncbi.nlm.nih.gov/IEB/ToolBox/SDKDOCS/SEQLOC.HTML#_Seq-loc:_Locations_on</rdfs:comment>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000023 -->

    <owl:Class rdf:about="&obo;CDAO_0000023">
        <rdfs:label rdf:datatype="&xsd;string">UnresolvedRootedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000012"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000017"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000024 -->

    <owl:Class rdf:about="&obo;CDAO_0000024">
        <rdfs:label rdf:datatype="&xsd;string">Branch</rdfs:label>
        <owl:equivalentClass rdf:resource="&obo;CDAO_0000099"/>
        <dc:description>&#39;Branch&#39; is the domain-specific synonym for an edge of a (Phylogenetic) Tree or Network.  Branches may have properties such as length and degree of support.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000025 -->

    <owl:Class rdf:about="&obo;CDAO_0000025">
        <rdfs:label rdf:datatype="&xsd;string">CharacterStateDataMatrixAnnotation</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000040"/>
        <dc:description>Meta-information associated with a character matrix, such as, for the case of a sequence alignment, the method of alignment.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000026 -->

    <owl:Class rdf:about="&obo;CDAO_0000026">
        <rdfs:label rdf:datatype="&xsd;string">AncestralNode</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000140"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000174"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000140"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000194"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000012"/>
                <owl:minQualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:minQualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000027 -->

    <owl:Class rdf:about="&obo;CDAO_0000027">
        <rdfs:label rdf:datatype="&xsd;string">UnresolvedUnrootedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000017"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000088"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000029 -->

    <owl:Class rdf:about="&obo;CDAO_0000029">
        <rdfs:label rdf:datatype="&xsd;string">UncertainStateDomain</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000091"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000216"/>
                <owl:someValuesFrom rdf:resource="&xsd;float"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000030 -->

    <owl:Class rdf:about="&obo;CDAO_0000030">
        <rdfs:label rdf:datatype="&xsd;string">ReconcileTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000172"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000110"/>
                <owl:minQualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">2</owl:minQualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000031 -->

    <owl:Class rdf:about="&obo;CDAO_0000031">
        <rdfs:label rdf:datatype="&xsd;string">Continuous</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000091"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000215"/>
                <owl:cardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:cardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>This class describes a continuous value. The link to the actual float value is through the property has_Value. It could have also other properties attached (e.g., has_Precision).</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000032 -->

    <owl:Class rdf:about="&obo;CDAO_0000032">
        <rdfs:label rdf:datatype="&xsd;string">AlignmentProcedure</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000025"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000033 -->

    <owl:Class rdf:about="&obo;CDAO_0000033">
        <rdfs:label rdf:datatype="&xsd;string">Dichotomy</rdfs:label>
        <owl:equivalentClass rdf:resource="&obo;CDAO_0000124"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000026"/>
        <owl:disjointWith rdf:resource="&obo;CDAO_0000042"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000034 -->

    <owl:Class rdf:about="&obo;CDAO_0000034">
        <rdfs:label rdf:datatype="&xsd;string">Molecular</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000039"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000035 -->

    <owl:Class rdf:about="&obo;CDAO_0000035">
        <rdfs:label rdf:datatype="&xsd;string">ContinuousCharParsimonyLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000009"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000039 -->

    <owl:Class rdf:about="&obo;CDAO_0000039">
        <rdfs:label rdf:datatype="&xsd;string">Categorical</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000091"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000040 -->

    <owl:Class rdf:about="&obo;CDAO_0000040">
        <rdfs:label rdf:datatype="&xsd;string">CDAOAnnotation</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:comment>Its possible that this base class should be discarded and that annotations should inherit from an imported base class if one exists.</rdfs:comment>
        <dc:description>The base class of annotations in CDAO.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000041 -->

    <owl:Class rdf:about="&obo;CDAO_0000041">
        <rdfs:label rdf:datatype="&xsd;string">originationEvent</rdfs:label>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000190"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000097"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000042 -->

    <owl:Class rdf:about="&obo;CDAO_0000042">
        <rdfs:label rdf:datatype="&xsd;string">Polytomy</rdfs:label>
        <owl:equivalentClass>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000177"/>
                <owl:minCardinality rdf:datatype="&xsd;nonNegativeInteger">3</owl:minCardinality>
            </owl:Restriction>
        </owl:equivalentClass>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000026"/>
        <owl:disjointWith rdf:resource="&obo;CDAO_0000124"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000043 -->

    <owl:Class rdf:about="&obo;CDAO_0000043">
        <rdfs:label rdf:datatype="&xsd;string">PolymorphicStateDomain</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000091"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000216"/>
                <owl:hasValue rdf:datatype="&xsd;float">1.0</owl:hasValue>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000044 -->

    <owl:Class rdf:about="&obo;CDAO_0000044">
        <rdfs:label rdf:datatype="&xsd;string">TreeAnnotation</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000040"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000157"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000110"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000045 -->

    <owl:Class rdf:about="&obo;CDAO_0000045">
        <rdfs:label rdf:datatype="&xsd;string">Standard</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000039"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000046 -->

    <owl:Class rdf:about="&obo;CDAO_0000046">
        <rdfs:label rdf:datatype="&xsd;string">EdgeLength</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000101"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000176"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000063"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000215"/>
                <owl:someValuesFrom rdf:resource="&rdfs;Literal"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:comment>Its possible that this should not be classed as an &#39;annotation&#39; since it contains data rather than meta-data.</rdfs:comment>
        <dc:description>The length of an edge (branch) of a Tree or Network, typically in units of evolutionary changes in character-state per character.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000047 -->

    <owl:Class rdf:about="&obo;CDAO_0000047">
        <rdfs:label rdf:datatype="&xsd;string">RibonucleotideResidue</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000034"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000048 -->

    <owl:Class rdf:about="&obo;CDAO_0000048">
        <rdfs:label rdf:datatype="&xsd;string">Clade</rdfs:label>
        <owl:equivalentClass rdf:resource="&obo;CDAO_0000129"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000049 -->

    <owl:Class rdf:about="&obo;CDAO_0000049">
        <rdfs:label rdf:datatype="&xsd;string">DiscreteCharParsimonyLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000100"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000050 -->

    <owl:Class rdf:about="&obo;CDAO_0000050">
        <rdfs:label rdf:datatype="&xsd;string">MolecularStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000089"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000051 -->

    <owl:Class rdf:about="&obo;CDAO_0000051">
        <rdfs:label rdf:datatype="&xsd;string">PolyphyleticGroup</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000006"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000188"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000140"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <owl:disjointWith rdf:resource="&obo;CDAO_0000127"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000052 -->

    <owl:Class rdf:about="&obo;CDAO_0000052">
        <rdfs:label rdf:datatype="&xsd;string">NexusDataBlock</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000107"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000053 -->

    <owl:Class rdf:about="&obo;CDAO_0000053">
        <rdfs:label rdf:datatype="&xsd;string">BranchingNode</rdfs:label>
        <owl:equivalentClass>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000177"/>
                <owl:minCardinality rdf:datatype="&xsd;nonNegativeInteger">2</owl:minCardinality>
            </owl:Restriction>
        </owl:equivalentClass>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000026"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000055 -->

    <owl:Class rdf:about="&obo;CDAO_0000055">
        <rdfs:label rdf:datatype="&xsd;string">Compound</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000039"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000056 -->

    <owl:Class rdf:about="&obo;CDAO_0000056">
        <rdfs:label rdf:datatype="&xsd;string">CharacterStateDataMatrix</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000178"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000025"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000208"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000138"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000142"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000071"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>A matrix of character-state data, typically containing observed data, though in some cases the states in the matrix might be simulated or hypothetical. Synonyms: character Data matrix, character-state matrix</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000057 -->

    <owl:Class rdf:about="&obo;CDAO_0000057">
        <rdfs:label rdf:datatype="&xsd;string">RibonucleotideResidueStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000050"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000058 -->

    <owl:Class rdf:about="&obo;CDAO_0000058">
        <rdfs:label rdf:datatype="&xsd;string">TimeCalibratedLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000063"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000059 -->

    <owl:Class rdf:about="&obo;CDAO_0000059">
        <rdfs:label rdf:datatype="&xsd;string">SetOfNodes</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000118"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000198"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000140"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000060 -->

    <owl:Class rdf:about="&obo;CDAO_0000060">
        <rdfs:label rdf:datatype="&xsd;string">MRCANode</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000080"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000163"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000118"/>
                <owl:minQualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:minQualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000061 -->

    <owl:Class rdf:about="&obo;CDAO_0000061">
        <rdfs:label rdf:datatype="&xsd;string">FASTADataMatrix</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000107"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000062 -->

    <owl:Class rdf:about="&obo;CDAO_0000062">
        <rdfs:label rdf:datatype="&xsd;string">evolutionaryTransition</rdfs:label>
        <owl:equivalentClass rdf:resource="&obo;CDAO_0000065"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000097"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000159"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000091"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000186"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000091"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000169"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000071"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000063 -->

    <owl:Class rdf:about="&obo;CDAO_0000063">
        <rdfs:label rdf:datatype="&xsd;string">EdgeLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000064 -->

    <owl:Class rdf:about="&obo;CDAO_0000064">
        <rdfs:label rdf:datatype="&xsd;string">cladogeneticChange</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000097"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000065 -->

    <owl:Class rdf:about="&obo;CDAO_0000065">
        <rdfs:label rdf:datatype="&xsd;string">anageneticChange</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000097"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000066 -->

    <owl:Class rdf:about="&obo;CDAO_0000066">
        <rdfs:label rdf:datatype="&xsd;string">TUAnnotation</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000040"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000067 -->

    <owl:Class rdf:about="&obo;CDAO_0000067">
        <rdfs:label rdf:datatype="&xsd;string">PhyloTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000074"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000068 -->

    <owl:Class rdf:about="&obo;CDAO_0000068">
        <rdfs:label rdf:datatype="&xsd;string">ContinuousCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000071"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000207"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000019"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000153"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000019"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000069 -->

    <owl:Class rdf:about="&obo;CDAO_0000069">
        <rdfs:label rdf:datatype="&xsd;string">PHYLIPTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000074"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000070 -->

    <owl:Class rdf:about="&obo;CDAO_0000070">
        <rdfs:label rdf:datatype="&xsd;string">Subtree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000155"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000110"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000071 -->

    <owl:Class rdf:about="&obo;CDAO_0000071">
        <rdfs:label rdf:datatype="&xsd;string">Character</rdfs:label>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000153"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000098"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:comment rdf:datatype="&xsd;string">Traits shown to be relevant for phylogenetic classification</rdfs:comment>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000072 -->

    <owl:Class rdf:about="&obo;CDAO_0000072">
        <rdfs:label rdf:datatype="&xsd;string">GalledTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000006"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000073 -->

    <owl:Class rdf:about="&obo;CDAO_0000073">
        <rdfs:label rdf:datatype="&xsd;string">SpeciesTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000074 -->

    <owl:Class rdf:about="&obo;CDAO_0000074">
        <rdfs:label rdf:datatype="&xsd;string">TreeFormat</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000044"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000075 -->

    <owl:Class rdf:about="&obo;CDAO_0000075">
        <rdfs:label rdf:datatype="&xsd;string">StandardCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000111"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000076 -->

    <owl:Class rdf:about="&obo;CDAO_0000076">
        <rdfs:label rdf:datatype="&xsd;string">AminoAcidResidue</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000034"/>
        <dc:description>This class will be declared equivalent ot the amino acid class description imported</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000077 -->

    <owl:Class rdf:about="&obo;CDAO_0000077">
        <rdfs:label rdf:datatype="&xsd;string">geneDuplication</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000064"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000078 -->

    <owl:Class rdf:about="&obo;CDAO_0000078">
        <rdfs:label rdf:datatype="&xsd;string">CompoundCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000111"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000170"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000136"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000142"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000071"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000153"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000136"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>A character that could be divided into separate characters but is not due to the non-independence of changes that would result, e.g., as in the case of a subsequence that is either present or absent as a block.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000079 -->

    <owl:Class rdf:about="&obo;CDAO_0000079">
        <rdfs:label rdf:datatype="&xsd;string">SIMMAPTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000074"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000080 -->

    <owl:Class rdf:about="&obo;CDAO_0000080">
        <rdfs:label rdf:datatype="&xsd;string">CommonAncestralNode</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000026"/>
        <rdfs:subClassOf>
            <owl:Class>
                <owl:unionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000053"/>
                    <owl:Restriction>
                        <owl:onProperty rdf:resource="&obo;CDAO_0000174"/>
                        <owl:someValuesFrom rdf:resource="&obo;CDAO_0000053"/>
                    </owl:Restriction>
                </owl:unionOf>
            </owl:Class>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000081 -->

    <owl:Class rdf:about="&obo;CDAO_0000081">
        <rdfs:label rdf:datatype="&xsd;string">NewickTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000074"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000082 -->

    <owl:Class rdf:about="&obo;CDAO_0000082">
        <rdfs:label rdf:datatype="&xsd;string">TimeProportionalLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000063"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000083 -->

    <owl:Class rdf:about="&obo;CDAO_0000083">
        <rdfs:label rdf:datatype="&xsd;string">DiscreteCharDistanceLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000100"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000084 -->

    <owl:Class rdf:about="&obo;CDAO_0000084">
        <rdfs:label rdf:datatype="&xsd;string">StarTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000012"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000149"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000108"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000085 -->

    <owl:Class rdf:about="&obo;CDAO_0000085">
        <rdfs:label rdf:datatype="&xsd;string">FullyResolvedUnrootedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000018"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000088"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000086 -->

    <owl:Class rdf:about="&obo;CDAO_0000086">
        <rdfs:label rdf:datatype="&xsd;string">ParaphyleticGroup</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000127"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000199"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000070"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <owl:disjointWith rdf:resource="&obo;CDAO_0000129"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000087 -->

    <owl:Class rdf:about="&obo;CDAO_0000087">
        <rdfs:label rdf:datatype="&xsd;string">geneticEvent</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000041"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000088 -->

    <owl:Class rdf:about="&obo;CDAO_0000088">
        <rdfs:label rdf:datatype="&xsd;string">UnrootedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000089 -->

    <owl:Class rdf:about="&obo;CDAO_0000089">
        <rdfs:label rdf:datatype="&xsd;string">CategoricalStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000098"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000090 -->

    <owl:Class rdf:about="&obo;CDAO_0000090">
        <rdfs:label rdf:datatype="&xsd;string">DiscreteCharLikelihoodLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000100"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000091 -->

    <owl:Class rdf:about="&obo;CDAO_0000091">
        <rdfs:label rdf:datatype="&xsd;string">CharacterStateDomain</rdfs:label>
        <dc:description>The universe of possible states for a particular type of character, e.g., the states of an Amino_Acid character come from the Amino_Acid domain.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000092 -->

    <owl:Class rdf:about="&obo;CDAO_0000092">
        <rdfs:label rdf:datatype="&xsd;string">CoordinateList</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000022"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000093 -->

    <owl:Class rdf:about="&obo;CDAO_0000093">
        <rdfs:label rdf:datatype="&xsd;string">GammaDistribution</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000020"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000094 -->

    <owl:Class rdf:about="&obo;CDAO_0000094">
        <rdfs:label rdf:datatype="&xsd;string">DesoxiRibonucleotideResidueCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000115"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000195"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000002"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000153"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000002"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000095 -->

    <owl:Class rdf:about="&obo;CDAO_0000095">
        <rdfs:label rdf:datatype="&xsd;string">CoordinateRange</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000022"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000096 -->

    <owl:Class rdf:about="&obo;CDAO_0000096">
        <rdfs:label rdf:datatype="&xsd;string">ReticulateEvolution</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000006"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000097 -->

    <owl:Class rdf:about="&obo;CDAO_0000097">
        <rdfs:label rdf:datatype="&xsd;string">hereditaryChange</rdfs:label>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000186"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000091"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000159"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000091"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000169"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000071"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000098 -->

    <owl:Class rdf:about="&obo;CDAO_0000098">
        <rdfs:label rdf:datatype="&xsd;string">CharacterStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000205"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000071"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000191"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000138"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>The instance of a given character for a given TU.  Its state is an object property drawn from a particular character state domain, e.g., the state of an Amino_Acid_State_Datum is an object property drawn from the domain Amino_Acid.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000099 -->

    <owl:Class rdf:about="&obo;CDAO_0000099">
        <rdfs:label rdf:datatype="&xsd;string">Edge</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000193"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000101"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000162"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000140"/>
                <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">2</owl:qualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>An edge connecting two nodes in a (Phylogenetic) Tree or Network, also known as a &#39;branch&#39;.  Edges may have attributes such as length, degree of support, and direction.  An edge can be a surrogate for a &#39;split&#39; or bipartition, since each edge in a tree divides the terminal nodes into two sets.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000100 -->

    <owl:Class rdf:about="&obo;CDAO_0000100">
        <rdfs:label rdf:datatype="&xsd;string">DiscreteCharacterLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000063"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000101 -->

    <owl:Class rdf:about="&obo;CDAO_0000101">
        <rdfs:label rdf:datatype="&xsd;string">EdgeAnnotation</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000040"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000102 -->

    <owl:Class rdf:about="&obo;CDAO_0000102">
        <rdfs:label rdf:datatype="&xsd;string">FullyResolvedRootedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000012"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000018"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000178"/>
                <owl:allValuesFrom>
                    <owl:Class>
                        <owl:unionOf rdf:parseType="Collection">
                            <rdf:Description rdf:about="&obo;CDAO_0000033"/>
                            <rdf:Description rdf:about="&obo;CDAO_0000099"/>
                            <rdf:Description rdf:about="&obo;CDAO_0000108"/>
                        </owl:unionOf>
                    </owl:Class>
                </owl:allValuesFrom>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000103 -->

    <owl:Class rdf:about="&obo;CDAO_0000103">
        <rdfs:label rdf:datatype="&xsd;string">GrafenLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000063"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000104 -->

    <owl:Class rdf:about="&obo;CDAO_0000104">
        <rdfs:label rdf:datatype="&xsd;string">CoordinateSystem</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <dc:description>A reference to an external coordinate system.  Coordinates for data must refer to some such external coordinate system.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000105 -->

    <owl:Class rdf:about="&obo;CDAO_0000105">
        <rdfs:label rdf:datatype="&xsd;string">GenBankDataMatrix</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000107"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000107 -->

    <owl:Class rdf:about="&obo;CDAO_0000107">
        <rdfs:label rdf:datatype="&xsd;string">DataMatrixFormat</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000025"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000108 -->

    <owl:Class rdf:about="&obo;CDAO_0000108">
        <rdfs:label rdf:datatype="&xsd;string">TerminalNode</rdfs:label>
        <owl:equivalentClass>
            <owl:Class>
                <owl:intersectionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000140"/>
                    <owl:Restriction>
                        <owl:onProperty rdf:resource="&obo;CDAO_0000149"/>
                        <owl:allValuesFrom>
                            <owl:Class>
                                <owl:complementOf rdf:resource="&obo;CDAO_0000140"/>
                            </owl:Class>
                        </owl:allValuesFrom>
                    </owl:Restriction>
                </owl:intersectionOf>
            </owl:Class>
        </owl:equivalentClass>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000140"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000109 -->

    <owl:Class rdf:about="&obo;CDAO_0000109">
        <rdfs:label rdf:datatype="&xsd;string">RibonucleotideResidueCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000115"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000153"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000057"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000110 -->

    <owl:Class rdf:about="&obo;CDAO_0000110">
        <rdfs:label rdf:datatype="&xsd;string">Tree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000006"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000111 -->

    <owl:Class rdf:about="&obo;CDAO_0000111">
        <rdfs:label rdf:datatype="&xsd;string">CategoricalCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000071"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000112 -->

    <owl:Class rdf:about="&obo;CDAO_0000112">
        <rdfs:label rdf:datatype="&xsd;string">AminoAcidResidueStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000050"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000113 -->

    <owl:Class rdf:about="&obo;CDAO_0000113">
        <rdfs:label rdf:datatype="&xsd;string">PHYLIPDataMatrix</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000107"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000114 -->

    <owl:Class rdf:about="&obo;CDAO_0000114">
        <rdfs:label rdf:datatype="&xsd;string">ContinuousCharLikelihoodLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000009"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000115 -->

    <owl:Class rdf:about="&obo;CDAO_0000115">
        <rdfs:label rdf:datatype="&xsd;string">MolecularCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000111"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000116 -->

    <owl:Class rdf:about="&obo;CDAO_0000116">
        <rdfs:label rdf:datatype="&xsd;string">hereditaryPersistance</rdfs:label>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000190"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000097"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000117 -->

    <owl:Class rdf:about="&obo;CDAO_0000117">
        <rdfs:label rdf:datatype="&xsd;string">SetOfCharacters</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000118"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000118 -->

    <owl:Class rdf:about="&obo;CDAO_0000118">
        <rdfs:label rdf:datatype="&xsd;string">SetOfThings</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000198"/>
                <owl:allValuesFrom>
                    <owl:Class>
                        <owl:unionOf rdf:parseType="Collection">
                            <rdf:Description rdf:about="&obo;CDAO_0000071"/>
                            <rdf:Description rdf:about="&obo;CDAO_0000117"/>
                        </owl:unionOf>
                    </owl:Class>
                </owl:allValuesFrom>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>The class is used to describe either colletions of characters or higher order grouping (e.g., groups of groups of characters). This extends the CharSet block of NEXUS.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000120 -->

    <owl:Class rdf:about="&obo;CDAO_0000120">
        <rdfs:label rdf:datatype="&xsd;string">Sequence</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000178"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000098"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000178"/>
                <owl:allValuesFrom>
                    <owl:Restriction>
                        <owl:onProperty rdf:resource="&obo;CDAO_0000151"/>
                        <owl:onClass rdf:resource="&obo;CDAO_0000022"/>
                        <owl:minQualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:minQualifiedCardinality>
                    </owl:Restriction>
                </owl:allValuesFrom>
            </owl:Restriction>
        </rdfs:subClassOf>
        <dc:description>A set of ordered states, typically the residues in a macromolecular sequence.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000121 -->

    <owl:Class rdf:about="&obo;CDAO_0000121">
        <rdfs:label rdf:datatype="&xsd;string">speciation</rdfs:label>
        <owl:equivalentClass rdf:resource="&obo;CDAO_0000122"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000064"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000122 -->

    <owl:Class rdf:about="&obo;CDAO_0000122">
        <rdfs:label rdf:datatype="&xsd;string">cladogenesis</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000064"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000124 -->

    <owl:Class rdf:about="&obo;CDAO_0000124">
        <rdfs:label rdf:datatype="&xsd;string">Bifurcation</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000026"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000177"/>
                <owl:cardinality rdf:datatype="&xsd;nonNegativeInteger">2</owl:cardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000125 -->

    <owl:Class rdf:about="&obo;CDAO_0000125">
        <rdfs:label rdf:datatype="&xsd;string">DiscreteCharBayesianLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000100"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000126 -->

    <owl:Class rdf:about="&obo;CDAO_0000126">
        <rdfs:label rdf:datatype="&xsd;string">TaxonomicLink</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000066"/>
        <dc:description>Link to an externally defined taxonomic hierarchy.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000127 -->

    <owl:Class rdf:about="&obo;CDAO_0000127">
        <rdfs:label rdf:datatype="&xsd;string">MonophyleticGroup</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000006"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000128 -->

    <owl:Class rdf:about="&obo;CDAO_0000128">
        <rdfs:label rdf:datatype="&xsd;string">molecularRecombination</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000132"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000129 -->

    <owl:Class rdf:about="&obo;CDAO_0000129">
        <rdfs:label rdf:datatype="&xsd;string">HolophyleticGroup</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000127"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000130 -->

    <owl:Class rdf:about="&obo;CDAO_0000130">
        <rdfs:label rdf:datatype="&xsd;string">FullyResolvedTree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000110"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000131 -->

    <owl:Class rdf:about="&obo;CDAO_0000131">
        <rdfs:label rdf:datatype="&xsd;string">AminoAcidResidueCharacter</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000115"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000168"/>
                <owl:someValuesFrom rdf:resource="&obo;CDAO_0000112"/>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000153"/>
                <owl:allValuesFrom rdf:resource="&obo;CDAO_0000112"/>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000132 -->

    <owl:Class rdf:about="&obo;CDAO_0000132">
        <rdfs:label rdf:datatype="&xsd;string">recombination</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000087"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000133 -->

    <owl:Class rdf:about="&obo;CDAO_0000133">
        <rdfs:label rdf:datatype="&xsd;string">DesoxiRibonucleotideResidue</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000034"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000134 -->

    <owl:Class rdf:about="&obo;CDAO_0000134">
        <rdfs:label rdf:datatype="&xsd;string">RootedSubtree</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000012"/>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000070"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000136 -->

    <owl:Class rdf:about="&obo;CDAO_0000136">
        <rdfs:label rdf:datatype="&xsd;string">CompoundStateDatum</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000089"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000137 -->

    <owl:Class rdf:about="&obo;CDAO_0000137">
        <rdfs:label rdf:datatype="&xsd;string">GapCost</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000007"/>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000138 -->

    <owl:Class rdf:about="&obo;CDAO_0000138">
        <rdfs:label rdf:datatype="&xsd;string">TU</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <dc:description>A unit of analysis that may be tied to a node in a tree and to a row in a character matrix.  It subsumes the traditional concepts of &#39;OTU&#39; and &#39;HTU&#39;.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000139 -->

    <owl:Class rdf:about="&obo;CDAO_0000139">
        <rdfs:label rdf:datatype="&xsd;string">DirectedEdge</rdfs:label>
        <owl:equivalentClass>
            <owl:Class>
                <owl:intersectionOf rdf:parseType="Collection">
                    <rdf:Description rdf:about="&obo;CDAO_0000099"/>
                    <owl:Restriction>
                        <owl:onProperty rdf:resource="&obo;CDAO_0000201"/>
                        <owl:onClass rdf:resource="&obo;CDAO_0000140"/>
                        <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
                    </owl:Restriction>
                    <owl:Restriction>
                        <owl:onProperty rdf:resource="&obo;CDAO_0000209"/>
                        <owl:onClass rdf:resource="&obo;CDAO_0000140"/>
                        <owl:qualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:qualifiedCardinality>
                    </owl:Restriction>
                </owl:intersectionOf>
            </owl:Class>
        </owl:equivalentClass>
        <dc:description>A directed edge. Rooted trees have directed edges. The direction is specified by way of the parent and child relationships of nodes that the edge connects.</dc:description>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000140 -->

    <owl:Class rdf:about="&obo;CDAO_0000140">
        <rdfs:label rdf:datatype="&xsd;string">Node</rdfs:label>
        <rdfs:subClassOf rdf:resource="&owl;Thing"/>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000143"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000139"/>
                <owl:maxQualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:maxQualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
        <rdfs:subClassOf>
            <owl:Restriction>
                <owl:onProperty rdf:resource="&obo;CDAO_0000194"/>
                <owl:onClass rdf:resource="&obo;CDAO_0000006"/>
                <owl:minQualifiedCardinality rdf:datatype="&xsd;nonNegativeInteger">1</owl:minQualifiedCardinality>
            </owl:Restriction>
        </rdfs:subClassOf>
    </owl:Class>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000141 -->

    <owl:Class rdf:about="&obo;CDAO_0000141">
        <rdfs:label rdf:datatype="&xsd;string">ContinuousCharDistanceLengthType</rdfs:label>
        <rdfs:subClassOf rdf:resource="&obo;CDAO_0000009"/>
    </owl:Class>
    


    <!--
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // Individuals
    //
    ///////////////////////////////////////////////////////////////////////////////////////
     -->

    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000220 -->

    <owl:Thing rdf:about="&obo;CDAO_0000220">
        <rdf:type rdf:resource="&obo;CDAO_0000133"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">dA</rdfs:label>
    </owl:Thing>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000221 -->

    <owl:Thing rdf:about="&obo;CDAO_0000221">
        <rdf:type rdf:resource="&obo;CDAO_0000015"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">absent</rdfs:label>
    </owl:Thing>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000222 -->

    <owl:Thing rdf:about="&obo;CDAO_0000222">
        <rdf:type rdf:resource="&obo;CDAO_0000015"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">unknown</rdfs:label>
    </owl:Thing>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000223 -->

    <owl:Thing rdf:about="&obo;CDAO_0000223">
        <rdf:type rdf:resource="&obo;CDAO_0000015"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">gap</rdfs:label>
    </owl:Thing>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000224 -->

    <owl:Thing rdf:about="&obo;CDAO_0000224">
        <rdf:type rdf:resource="&obo;CDAO_0000133"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">dG</rdfs:label>
    </owl:Thing>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000225 -->

    <owl:Thing rdf:about="&obo;CDAO_0000225">
        <rdf:type rdf:resource="&obo;CDAO_0000057"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">rU</rdfs:label>
    </owl:Thing>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000226 -->

    <owl:Thing rdf:about="&obo;CDAO_0000226">
        <rdf:type rdf:resource="&obo;CDAO_0000133"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">dC</rdfs:label>
    </owl:Thing>
    


    <!-- http://purl.obolibrary.org/obo/CDAO_0000227 -->

    <owl:Thing rdf:about="&obo;CDAO_0000227">
        <rdf:type rdf:resource="&obo;CDAO_0000133"/>
        <rdf:type rdf:resource="&owl;NamedIndividual"/>
        <rdfs:label rdf:datatype="&xsd;string">dT</rdfs:label>
    </owl:Thing>
    


    <!--
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // General axioms
    //
    ///////////////////////////////////////////////////////////////////////////////////////
     -->

    <rdf:Description>
        <rdf:type rdf:resource="&owl;AllDisjointClasses"/>
        <owl:members rdf:parseType="Collection">
            <rdf:Description rdf:about="&obo;CDAO_0000068"/>
            <rdf:Description rdf:about="&obo;CDAO_0000094"/>
            <rdf:Description rdf:about="&obo;CDAO_0000131"/>
        </owl:members>
    </rdf:Description>
    <rdf:Description>
        <rdf:type rdf:resource="&owl;AllDisjointClasses"/>
        <owl:members rdf:parseType="Collection">
            <rdf:Description rdf:about="&obo;CDAO_0000002"/>
            <rdf:Description rdf:about="&obo;CDAO_0000019"/>
            <rdf:Description rdf:about="&obo;CDAO_0000112"/>
        </owl:members>
    </rdf:Description>
</rdf:RDF>



<!-- Generated by the OWL API (version 3.2.3.1824) http://owlapi.sourceforge.net -->

'''

cdao_elements = {}

root = ET.fromstring(cdao_owl)
for node_type in 'ObjectProperty', 'Class', 'DatatypeProperty':
    for element in root.findall('{http://www.w3.org/2002/07/owl#}%s' % node_type):
        obo = element.attrib[
            '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about'].split('/')[-1]
        cdao = element.find(
            '{http://www.w3.org/2000/01/rdf-schema#}label').text
        cdao_elements[cdao] = obo
