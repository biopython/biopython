# Copyright (C) 2009 by Eric Talevich (eric.talevich@gmail.com)
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""PhyloXML reader/parser, writer, and associated functions.

Instantiates tree elements from a parsed PhyloXML file, and constructs an XML
file from a ``Bio.Phylo.PhyloXML`` object.

About capitalization:
 - phyloXML means the file format specification
 - PhyloXML means the Biopython module ``Bio.Phylo.PhyloXML`` and its classes
 - Phyloxml means the top-level class used by ``PhyloXMLIO.read`` (but not
   ``Bio.Phylo.read``!), containing a list of Phylogenies (objects derived from
   ``BaseTree.Tree``)

"""

from xml.etree import ElementTree

from Bio.Phylo import PhyloXML as PX


# Recognize the phyloXML namespace when parsing
# See http://effbot.org/zone/element-namespaces.htm
NAMESPACES = {"phy": "http://www.phyloxml.org"}

try:
    register_namespace = ElementTree.register_namespace
except AttributeError:
    if not hasattr(ElementTree, "_namespace_map"):
        # cElementTree needs the pure-Python xml.etree.ElementTree
        from xml.etree import ElementTree as ET_py

        ElementTree._namespace_map = ET_py._namespace_map

    def register_namespace(prefix, uri):
        """Set the namespace for ElementTree."""
        ElementTree._namespace_map[uri] = prefix


for prefix, uri in NAMESPACES.items():
    register_namespace(prefix, uri)

# Tell ElementTree how to write to text handles
DEFAULT_ENCODING = "unicode"


class PhyloXMLError(Exception):
    """Exception raised when PhyloXML object construction cannot continue.

    XML syntax errors will be found and raised by the underlying ElementTree
    module; this exception is for valid XML that breaks the phyloXML
    specification.
    """

    pass


# ---------------------------------------------------------
# Public API


def read(file):
    """Parse a phyloXML file or stream and build a tree of Biopython objects.

    The children of the root node are phylogenies and possibly other arbitrary
    (non-phyloXML) objects.

    :returns: a single ``Bio.Phylo.PhyloXML.Phyloxml`` object.

    """
    return Parser(file).read()


def parse(file):
    """Iterate over the phylogenetic trees in a phyloXML file.

    This ignores any additional data stored at the top level, but may be more
    memory-efficient than the ``read`` function.

    :returns: a generator of ``Bio.Phylo.PhyloXML.Phylogeny`` objects.

    """
    return Parser(file).parse()


def write(obj, file, encoding=DEFAULT_ENCODING, indent=True):
    """Write a phyloXML file.

    :Parameters:
        obj
            an instance of ``Phyloxml``, ``Phylogeny`` or ``BaseTree.Tree``,
            or an iterable of either of the latter two. The object will be
            converted to a Phyloxml object before serialization.
        file
            either an open handle or a file name.

    """

    def fix_single(tree):
        if isinstance(tree, PX.Phylogeny):
            return tree
        if isinstance(tree, PX.Clade):
            return tree.to_phylogeny()
        if isinstance(tree, PX.BaseTree.Tree):
            return PX.Phylogeny.from_tree(tree)
        if isinstance(tree, PX.BaseTree.Clade):
            return PX.Phylogeny.from_tree(PX.BaseTree.Tree(root=tree))
        else:
            raise ValueError("iterable must contain Tree or Clade types")

    if isinstance(obj, PX.Phyloxml):
        pass
    elif isinstance(obj, (PX.BaseTree.Tree, PX.BaseTree.Clade)):
        obj = fix_single(obj).to_phyloxml()
    elif hasattr(obj, "__iter__"):
        obj = PX.Phyloxml({}, phylogenies=(fix_single(t) for t in obj))
    else:
        raise ValueError(
            "First argument must be a Phyloxml, Phylogeny, "
            "Tree, or iterable of Trees or Phylogenies."
        )
    return Writer(obj).write(file, encoding=encoding, indent=indent)


# ---------------------------------------------------------
# Functions I wish ElementTree had


def _local(tag):
    """Extract the local tag from a namespaced tag name (PRIVATE)."""
    if tag[0] == "{":
        return tag[tag.index("}") + 1 :]
    return tag


def _split_namespace(tag):
    """Split a tag into namespace and local tag strings (PRIVATE)."""
    try:
        return tag[1:].split("}", 1)
    except ValueError:
        return ("", tag)


def _ns(tag, namespace=NAMESPACES["phy"]):
    """Format an XML tag with the given namespace (PRIVATE)."""
    return f"{{{namespace}}}{tag}"


def _get_child_as(parent, tag, construct):
    """Find a child node by tag, and pass it through a constructor (PRIVATE).

    Returns None if no matching child is found.
    """
    child = parent.find(_ns(tag))
    if child is not None:
        return construct(child)


def _get_child_text(parent, tag, construct=str):
    """Find a child node by tag; pass its text through a constructor (PRIVATE).

    Returns None if no matching child is found.
    """
    child = parent.find(_ns(tag))
    if child is not None and child.text:
        return construct(child.text)


def _get_children_as(parent, tag, construct):
    """Find child nodes by tag; pass each through a constructor (PRIVATE).

    Returns an empty list if no matching child is found.
    """
    return [construct(child) for child in parent.findall(_ns(tag))]


def _get_children_text(parent, tag, construct=str):
    """Find child nodes by tag; pass each node's text through a constructor (PRIVATE).

    Returns an empty list if no matching child is found.
    """
    return [construct(child.text) for child in parent.findall(_ns(tag)) if child.text]


def _indent(elem, level=0):
    """Add line breaks and indentation to ElementTree in-place (PRIVATE).

    Sources:
     - http://effbot.org/zone/element-lib.htm#prettyprint
     - http://infix.se/2007/02/06/gentlemen-indent-your-xml

    """
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        for e in elem:
            _indent(e, level + 1)
            if not e.tail or not e.tail.strip():
                e.tail = i + "  "
        if not e.tail or not e.tail.strip():
            e.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


# ---------------------------------------------------------
# INPUT
# ---------------------------------------------------------


def _str2bool(text):
    """Convert string to boolean (PRIVATE)."""
    if text == "true" or text == "1":
        return True
    if text == "false" or text == "0":
        return False
    raise ValueError("String could not be converted to boolean: " + text)


def _dict_str2bool(dct, keys):
    """Return a new dictionary where string values are replaced with booleans (PRIVATE)."""
    out = dct.copy()
    for key in keys:
        if key in out:
            out[key] = _str2bool(out[key])
    return out


def _int(text):
    """Return text as an integer (PRIVATE)."""
    if text is not None:
        try:
            return int(text)
        except Exception:
            return None


def _float(text):
    """Return text as a float (PRIVATE)."""
    if text is not None:
        try:
            return float(text)
        except Exception:
            return None


def _collapse_wspace(text):
    """Replace all spans of whitespace with a single space character (PRIVATE).

    Also remove leading and trailing whitespace. See "Collapse Whitespace
    Policy" in the phyloXML spec glossary:
    http://phyloxml.org/documentation/version_100/phyloxml.xsd.html#Glossary
    """
    if text is not None:
        return " ".join(text.split())


# NB: Not currently used
def _replace_wspace(text):
    """Replace tab, LF and CR characters with spaces, but don't collapse (PRIVATE).

    See "Replace Whitespace Policy" in the phyloXML spec glossary:
    http://phyloxml.org/documentation/version_100/phyloxml.xsd.html#Glossary
    """
    for char in ("\t", "\n", "\r"):
        if char in text:
            text = text.replace(char, " ")
    return text


class Parser:
    """Methods for parsing all phyloXML nodes from an XML stream.

    To minimize memory use, the tree of ElementTree parsing events is cleared
    after completing each phylogeny, clade, and top-level 'other' element.
    Elements below the clade level are kept in memory until parsing of the
    current clade is finished -- this shouldn't be a problem because clade is
    the only recursive element, and non-clade nodes below this level are of
    bounded size.
    """

    def __init__(self, file):
        """Initialize the class."""
        # Get an iterable context for XML parsing events
        context = iter(ElementTree.iterparse(file, events=("start", "end")))
        event, root = next(context)
        self.root = root
        self.context = context

    def read(self):
        """Parse the phyloXML file and create a single Phyloxml object."""
        phyloxml = PX.Phyloxml({_local(key): val for key, val in self.root.items()})
        other_depth = 0
        for event, elem in self.context:
            namespace, localtag = _split_namespace(elem.tag)
            if event == "start":
                if namespace != NAMESPACES["phy"]:
                    other_depth += 1
                    continue
                if localtag == "phylogeny":
                    phylogeny = self._parse_phylogeny(elem)
                    phyloxml.phylogenies.append(phylogeny)
            if event == "end" and namespace != NAMESPACES["phy"]:
                # Deal with items not specified by phyloXML
                other_depth -= 1
                if other_depth == 0:
                    # We're directly under the root node -- evaluate
                    otr = self.other(elem, namespace, localtag)
                    phyloxml.other.append(otr)
                    self.root.clear()
        return phyloxml

    def parse(self):
        """Parse the phyloXML file incrementally and return each phylogeny."""
        phytag = _ns("phylogeny")
        for event, elem in self.context:
            if event == "start" and elem.tag == phytag:
                yield self._parse_phylogeny(elem)

    # Special parsing cases -- incremental, using self.context

    def _parse_phylogeny(self, parent):
        """Parse a single phylogeny within the phyloXML tree (PRIVATE).

        Recursively builds a phylogenetic tree with help from parse_clade, then
        clears the XML event history for the phylogeny element and returns
        control to the top-level parsing function.
        """
        phylogeny = PX.Phylogeny(
            **_dict_str2bool(parent.attrib, ["rooted", "rerootable"])
        )
        list_types = {
            # XML tag, plural attribute
            "confidence": "confidences",
            "property": "properties",
            "clade_relation": "clade_relations",
            "sequence_relation": "sequence_relations",
        }
        for event, elem in self.context:
            namespace, tag = _split_namespace(elem.tag)
            if event == "start" and tag == "clade":
                if phylogeny.root is not None:
                    raise ValueError("Phylogeny object should only have 1 clade")
                phylogeny.root = self._parse_clade(elem)
                continue
            if event == "end":
                if tag == "phylogeny":
                    parent.clear()
                    break
                # Handle the other non-recursive children
                if tag in list_types:
                    getattr(phylogeny, list_types[tag]).append(getattr(self, tag)(elem))
                # Complex types
                elif tag in ("date", "id"):
                    setattr(phylogeny, tag, getattr(self, tag)(elem))
                # Simple types
                elif tag in ("name", "description"):
                    setattr(phylogeny, tag, _collapse_wspace(elem.text))
                # Unknown tags
                elif namespace != NAMESPACES["phy"]:
                    phylogeny.other.append(self.other(elem, namespace, tag))
                    parent.clear()
                else:
                    # NB: This shouldn't happen in valid files
                    raise PhyloXMLError("Misidentified tag: " + tag)
        return phylogeny

    _clade_complex_types = ["color", "events", "binary_characters", "date"]
    _clade_list_types = {
        "confidence": "confidences",
        "distribution": "distributions",
        "reference": "references",
        "property": "properties",
    }
    _clade_tracked_tags = (
        set(_clade_complex_types)
        .union(_clade_list_types.keys())
        .union(["branch_length", "name", "node_id", "width"])
    )

    def _parse_clade(self, parent):
        """Parse a Clade node and its children, recursively (PRIVATE)."""
        clade = PX.Clade(**parent.attrib)
        if clade.branch_length is not None:
            clade.branch_length = float(clade.branch_length)
        # NB: Only evaluate nodes at the current level
        tag_stack = []
        for event, elem in self.context:
            namespace, tag = _split_namespace(elem.tag)
            if event == "start":
                if tag == "clade":
                    clade.clades.append(self._parse_clade(elem))
                    continue
                if tag == "taxonomy":
                    clade.taxonomies.append(self._parse_taxonomy(elem))
                    continue
                if tag == "sequence":
                    clade.sequences.append(self._parse_sequence(elem))
                    continue
                if tag in self._clade_tracked_tags:
                    tag_stack.append(tag)
            if event == "end":
                if tag == "clade":
                    elem.clear()
                    break
                if tag != tag_stack[-1]:
                    continue
                tag_stack.pop()
                # Handle the other non-recursive children
                if tag in self._clade_list_types:
                    getattr(clade, self._clade_list_types[tag]).append(
                        getattr(self, tag)(elem)
                    )
                elif tag in self._clade_complex_types:
                    setattr(clade, tag, getattr(self, tag)(elem))
                elif tag == "branch_length":
                    # NB: possible collision with the attribute
                    if clade.branch_length is not None:
                        raise PhyloXMLError(
                            "Attribute branch_length was already set for this Clade."
                        )
                    clade.branch_length = _float(elem.text)
                elif tag == "width":
                    clade.width = _float(elem.text)
                elif tag == "name":
                    clade.name = _collapse_wspace(elem.text)
                elif tag == "node_id":
                    clade.node_id = PX.Id(
                        elem.text.strip(), elem.attrib.get("provider")
                    )
                elif namespace != NAMESPACES["phy"]:
                    clade.other.append(self.other(elem, namespace, tag))
                    elem.clear()
                else:
                    raise PhyloXMLError("Misidentified tag: " + tag)
        return clade

    def _parse_sequence(self, parent):
        """Parse a molecular sequence (PRIVATE)."""
        sequence = PX.Sequence(**parent.attrib)
        for event, elem in self.context:
            namespace, tag = _split_namespace(elem.tag)
            if event == "end":
                if tag == "sequence":
                    parent.clear()
                    break
                if tag in ("accession", "mol_seq", "uri", "domain_architecture"):
                    setattr(sequence, tag, getattr(self, tag)(elem))
                elif tag == "annotation":
                    sequence.annotations.append(self.annotation(elem))
                elif tag == "name":
                    sequence.name = _collapse_wspace(elem.text)
                elif tag in ("symbol", "location"):
                    setattr(sequence, tag, elem.text)
                elif namespace != NAMESPACES["phy"]:
                    sequence.other.append(self.other(elem, namespace, tag))
                    parent.clear()
        return sequence

    def _parse_taxonomy(self, parent):
        """Parse taxonomic information for a clade (PRIVATE)."""
        taxonomy = PX.Taxonomy(**parent.attrib)
        for event, elem in self.context:
            namespace, tag = _split_namespace(elem.tag)
            if event == "end":
                if tag == "taxonomy":
                    parent.clear()
                    break
                if tag in ("id", "uri"):
                    setattr(taxonomy, tag, getattr(self, tag)(elem))
                elif tag == "common_name":
                    taxonomy.common_names.append(_collapse_wspace(elem.text))
                elif tag == "synonym":
                    taxonomy.synonyms.append(elem.text)
                elif tag in ("code", "scientific_name", "authority", "rank"):
                    # ENH: check_str on rank
                    setattr(taxonomy, tag, elem.text)
                elif namespace != NAMESPACES["phy"]:
                    taxonomy.other.append(self.other(elem, namespace, tag))
                    parent.clear()
        return taxonomy

    def other(self, elem, namespace, localtag):
        """Create an Other object, a non-phyloXML element."""
        return PX.Other(
            localtag,
            namespace,
            elem.attrib,
            value=elem.text and elem.text.strip() or None,
            children=[
                self.other(child, *_split_namespace(child.tag)) for child in elem
            ],
        )

    # Complex types

    def accession(self, elem):
        """Create accession object."""
        return PX.Accession(elem.text.strip(), elem.get("source"))

    def annotation(self, elem):
        """Create annotation object."""
        return PX.Annotation(
            desc=_collapse_wspace(_get_child_text(elem, "desc")),
            confidence=_get_child_as(elem, "confidence", self.confidence),
            properties=_get_children_as(elem, "property", self.property),
            uri=_get_child_as(elem, "uri", self.uri),
            **elem.attrib,
        )

    def binary_characters(self, elem):
        """Create binary characters object."""

        def bc_getter(elem):
            """Get binary characters from subnodes."""
            return _get_children_text(elem, "bc")

        return PX.BinaryCharacters(
            type=elem.get("type"),
            gained_count=_int(elem.get("gained_count")),
            lost_count=_int(elem.get("lost_count")),
            present_count=_int(elem.get("present_count")),
            absent_count=_int(elem.get("absent_count")),
            # Flatten BinaryCharacterList sub-nodes into lists of strings
            gained=_get_child_as(elem, "gained", bc_getter),
            lost=_get_child_as(elem, "lost", bc_getter),
            present=_get_child_as(elem, "present", bc_getter),
            absent=_get_child_as(elem, "absent", bc_getter),
        )

    def clade_relation(self, elem):
        """Create clade relationship object."""
        return PX.CladeRelation(
            elem.get("type"),
            elem.get("id_ref_0"),
            elem.get("id_ref_1"),
            distance=elem.get("distance"),
            confidence=_get_child_as(elem, "confidence", self.confidence),
        )

    def color(self, elem):
        """Create branch color object."""
        red, green, blue = (
            _get_child_text(elem, color, int) for color in ("red", "green", "blue")
        )
        return PX.BranchColor(red, green, blue)

    def confidence(self, elem):
        """Create confidence object."""
        return PX.Confidence(_float(elem.text), elem.get("type"))

    def date(self, elem):
        """Create date object."""
        return PX.Date(
            unit=elem.get("unit"),
            desc=_collapse_wspace(_get_child_text(elem, "desc")),
            value=_get_child_text(elem, "value", float),
            minimum=_get_child_text(elem, "minimum", float),
            maximum=_get_child_text(elem, "maximum", float),
        )

    def distribution(self, elem):
        """Create geographic distribution object."""
        return PX.Distribution(
            desc=_collapse_wspace(_get_child_text(elem, "desc")),
            points=_get_children_as(elem, "point", self.point),
            polygons=_get_children_as(elem, "polygon", self.polygon),
        )

    def domain(self, elem):
        """Create protein domain object."""
        return PX.ProteinDomain(
            elem.text.strip(),
            int(elem.get("from")) - 1,
            int(elem.get("to")),
            confidence=_float(elem.get("confidence")),
            id=elem.get("id"),
        )

    def domain_architecture(self, elem):
        """Create domain architecture object."""
        return PX.DomainArchitecture(
            length=int(elem.get("length")),
            domains=_get_children_as(elem, "domain", self.domain),
        )

    def events(self, elem):
        """Create events object."""
        return PX.Events(
            type=_get_child_text(elem, "type"),
            duplications=_get_child_text(elem, "duplications", int),
            speciations=_get_child_text(elem, "speciations", int),
            losses=_get_child_text(elem, "losses", int),
            confidence=_get_child_as(elem, "confidence", self.confidence),
        )

    def id(self, elem):
        """Create identifier object."""
        provider = elem.get("provider") or elem.get("type")
        return PX.Id(elem.text.strip(), provider)

    def mol_seq(self, elem):
        """Create molecular sequence object."""
        is_aligned = elem.get("is_aligned")
        if is_aligned is not None:
            is_aligned = _str2bool(is_aligned)
        return PX.MolSeq(elem.text.strip(), is_aligned=is_aligned)

    def point(self, elem):
        """Create point object, coordinates of a point."""
        return PX.Point(
            elem.get("geodetic_datum"),
            _get_child_text(elem, "lat", float),
            _get_child_text(elem, "long", float),
            alt=_get_child_text(elem, "alt", float),
            alt_unit=elem.get("alt_unit"),
        )

    def polygon(self, elem):
        """Create polygon object, list of points."""
        return PX.Polygon(points=_get_children_as(elem, "point", self.point))

    def property(self, elem):
        """Create properties from external resources."""
        return PX.Property(
            elem.text.strip(),
            elem.get("ref"),
            elem.get("applies_to"),
            elem.get("datatype"),
            unit=elem.get("unit"),
            id_ref=elem.get("id_ref"),
        )

    def reference(self, elem):
        """Create literature reference object."""
        return PX.Reference(doi=elem.get("doi"), desc=_get_child_text(elem, "desc"))

    def sequence_relation(self, elem):
        """Create sequence relationship object, relationship between two sequences."""
        return PX.SequenceRelation(
            elem.get("type"),
            elem.get("id_ref_0"),
            elem.get("id_ref_1"),
            distance=_float(elem.get("distance")),
            confidence=_get_child_as(elem, "confidence", self.confidence),
        )

    def uri(self, elem):
        """Create uri object, expected to be a url."""
        return PX.Uri(
            elem.text.strip(),
            desc=_collapse_wspace(elem.get("desc")),
            type=elem.get("type"),
        )


# ---------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------


def _serialize(value):
    """Convert a Python primitive to a phyloXML-compatible string (PRIVATE)."""
    if isinstance(value, float):
        return str(value).upper()
    elif isinstance(value, bool):
        return str(value).lower()
    return str(value)


def _clean_attrib(obj, attrs):
    """Create a dictionary from an object's specified, non-None attributes (PRIVATE)."""
    out = {}
    for key in attrs:
        val = getattr(obj, key)
        if val is not None:
            out[key] = _serialize(val)
    return out


def _handle_complex(tag, attribs, subnodes, has_text=False):
    """Handle to serialize nodes with subnodes (PRIVATE)."""

    def wrapped(self, obj):
        """Wrap nodes and subnodes as elements."""
        elem = ElementTree.Element(tag, _clean_attrib(obj, attribs))
        for subn in subnodes:
            if isinstance(subn, str):
                # singular object: method and attribute names are the same
                if getattr(obj, subn) is not None:
                    elem.append(getattr(self, subn)(getattr(obj, subn)))
            else:
                # list: singular method, pluralized attribute name
                method, plural = subn
                for item in getattr(obj, plural):
                    elem.append(getattr(self, method)(item))
        if has_text:
            elem.text = _serialize(obj.value)
        return elem

    wrapped.__doc__ = f"Serialize a {tag} and its subnodes, in order."
    return wrapped


def _handle_simple(tag):
    """Handle to serialize simple nodes (PRIVATE)."""

    def wrapped(self, obj):
        """Wrap node as element."""
        elem = ElementTree.Element(tag)
        elem.text = _serialize(obj)
        return elem

    wrapped.__doc__ = f"Serialize a simple {tag} node."
    return wrapped


class Writer:
    """Methods for serializing a PhyloXML object to XML."""

    def __init__(self, phyloxml):
        """Build an ElementTree from a PhyloXML object."""
        assert isinstance(phyloxml, PX.Phyloxml), "Not a Phyloxml object"
        self._tree = ElementTree.ElementTree(self.phyloxml(phyloxml))

    def write(self, file, encoding=DEFAULT_ENCODING, indent=True):
        """Write PhyloXML to a file."""
        if indent:
            _indent(self._tree.getroot())
        self._tree.write(file, encoding)
        return len(self._tree.getroot())

    # Convert classes to ETree elements

    def phyloxml(self, obj):
        """Convert phyloxml to Etree element."""
        elem = ElementTree.Element("phyloxml", obj.attributes)  # Namespaces
        for tree in obj.phylogenies:
            elem.append(self.phylogeny(tree))
        for otr in obj.other:
            elem.append(self.other(otr))
        return elem

    def other(self, obj):
        """Convert other to Etree element."""
        elem = ElementTree.Element(_ns(obj.tag, obj.namespace), obj.attributes)
        elem.text = obj.value
        for child in obj.children:
            elem.append(self.other(child))
        return elem

    phylogeny = _handle_complex(
        "phylogeny",
        ("rooted", "rerootable", "branch_length_unit", "type"),
        (
            "name",
            "id",
            "description",
            "date",
            ("confidence", "confidences"),
            "clade",
            ("clade_relation", "clade_relations"),
            ("sequence_relation", "sequence_relations"),
            ("property", "properties"),
            ("other", "other"),
        ),
    )

    clade = _handle_complex(
        "clade",
        ("id_source",),
        (
            "name",
            "branch_length",
            ("confidence", "confidences"),
            "width",
            "color",
            "node_id",
            ("taxonomy", "taxonomies"),
            ("sequence", "sequences"),
            "events",
            "binary_characters",
            ("distribution", "distributions"),
            "date",
            ("reference", "references"),
            ("property", "properties"),
            ("clade", "clades"),
            ("other", "other"),
        ),
    )

    accession = _handle_complex("accession", ("source",), (), has_text=True)

    annotation = _handle_complex(
        "annotation",
        ("ref", "source", "evidence", "type"),
        ("desc", "confidence", ("property", "properties"), "uri"),
    )

    def binary_characters(self, obj):
        """Serialize a binary_characters node and its subnodes."""
        elem = ElementTree.Element(
            "binary_characters",
            _clean_attrib(
                obj,
                ("type", "gained_count", "lost_count", "present_count", "absent_count"),
            ),
        )
        for subn in ("gained", "lost", "present", "absent"):
            subelem = ElementTree.Element(subn)
            for token in getattr(obj, subn):
                subelem.append(self.bc(token))
            elem.append(subelem)
        return elem

    clade_relation = _handle_complex(
        "clade_relation", ("id_ref_0", "id_ref_1", "distance", "type"), ("confidence",)
    )

    color = _handle_complex("color", (), ("red", "green", "blue"))

    confidence = _handle_complex("confidence", ("type",), (), has_text=True)

    date = _handle_complex("date", ("unit",), ("desc", "value", "minimum", "maximum"))

    distribution = _handle_complex(
        "distribution", (), ("desc", ("point", "points"), ("polygon", "polygons"))
    )

    def domain(self, obj):
        """Serialize a domain node."""
        elem = ElementTree.Element(
            "domain", {"from": str(obj.start + 1), "to": str(obj.end)}
        )
        if obj.confidence is not None:
            elem.set("confidence", _serialize(obj.confidence))
        if obj.id is not None:
            elem.set("id", obj.id)
        elem.text = _serialize(obj.value)
        return elem

    domain_architecture = _handle_complex(
        "domain_architecture", ("length",), (("domain", "domains"),)
    )

    events = _handle_complex(
        "events", (), ("type", "duplications", "speciations", "losses", "confidence")
    )

    id = _handle_complex("id", ("provider",), (), has_text=True)

    mol_seq = _handle_complex("mol_seq", ("is_aligned",), (), has_text=True)

    node_id = _handle_complex("node_id", ("provider",), (), has_text=True)

    point = _handle_complex(
        "point", ("geodetic_datum", "alt_unit"), ("lat", "long", "alt")
    )

    polygon = _handle_complex("polygon", (), (("point", "points"),))

    property = _handle_complex(
        "property",
        ("ref", "unit", "datatype", "applies_to", "id_ref"),
        (),
        has_text=True,
    )

    reference = _handle_complex("reference", ("doi",), ("desc",))

    sequence = _handle_complex(
        "sequence",
        ("type", "id_ref", "id_source"),
        (
            "symbol",
            "accession",
            "name",
            "location",
            "mol_seq",
            "uri",
            ("annotation", "annotations"),
            "domain_architecture",
            ("other", "other"),
        ),
    )

    sequence_relation = _handle_complex(
        "sequence_relation",
        ("id_ref_0", "id_ref_1", "distance", "type"),
        ("confidence",),
    )

    taxonomy = _handle_complex(
        "taxonomy",
        ("id_source",),
        (
            "id",
            "code",
            "scientific_name",
            "authority",
            ("common_name", "common_names"),
            ("synonym", "synonyms"),
            "rank",
            "uri",
            ("other", "other"),
        ),
    )

    uri = _handle_complex("uri", ("desc", "type"), (), has_text=True)

    # Primitive types

    # Floating point
    alt = _handle_simple("alt")
    branch_length = _handle_simple("branch_length")
    lat = _handle_simple("lat")
    long = _handle_simple("long")
    maximum = _handle_simple("maximum")
    minimum = _handle_simple("minimum")
    value = _handle_simple("value")
    width = _handle_simple("width")

    # Integers
    blue = _handle_simple("blue")
    duplications = _handle_simple("duplications")
    green = _handle_simple("green")
    losses = _handle_simple("losses")
    red = _handle_simple("red")
    speciations = _handle_simple("speciations")

    # Strings
    bc = _handle_simple("bc")
    code = _handle_simple("code")
    common_name = _handle_simple("common_name")
    desc = _handle_simple("desc")
    description = _handle_simple("description")
    location = _handle_simple("location")
    name = _handle_simple("name")
    rank = _handle_simple("rank")
    scientific_name = _handle_simple("scientific_name")
    symbol = _handle_simple("symbol")
    synonym = _handle_simple("synonym")
    type = _handle_simple("type")
