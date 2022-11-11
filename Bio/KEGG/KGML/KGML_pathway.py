# Copyright 2013 by Leighton Pritchard.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Classes to represent a KGML Pathway Map.

The KGML definition is as of release KGML v0.7.2
(http://www.kegg.jp/kegg/xml/docs/)

Classes:
 - Pathway - Specifies graph information for the pathway map
 - Relation - Specifies a relationship between two proteins or KOs,
   or protein and compound. There is an implied direction to the
   relationship in some cases.
 - Reaction - A specific chemical reaction between a substrate and
   a product.
 - Entry - A node in the pathway graph
 - Graphics - Entry subelement describing its visual representation

"""

import time
from itertools import chain
from xml.dom import minidom
import xml.etree.ElementTree as ET


# Pathway
class Pathway:
    """Represents a KGML pathway from KEGG.

    Specifies graph information for the pathway map, as described in
    release KGML v0.7.2 (http://www.kegg.jp/kegg/xml/docs/)

    Attributes:
     - name - KEGGID of the pathway map
     - org - ko/ec/[org prefix]
     - number - map number (integer)
     - title - the map title
     - image - URL of the image map for the pathway
     - link - URL of information about the pathway
     - entries - Dictionary of entries in the pathway, keyed by node ID
     - reactions - Set of reactions in the pathway

    The name attribute has a restricted format, so we make it a property and
    enforce the formatting.

    The Pathway object is the only allowed route for adding/removing
    Entry, Reaction, or Relation elements.

    Entries are held in a dictionary and keyed by the node ID for the
    pathway graph - this allows for ready access via the Reaction/Relation
    etc. elements.  Entries must be added before reference by any other
    element.

    Reactions are held in a dictionary, keyed by node ID for the path.
    The elements referred to in the reaction must be added before the
    reaction itself.

    """

    def __init__(self):
        """Initialize the class."""
        self._name = ""
        self.org = ""
        self._number = None
        self.title = ""
        self.image = ""
        self.link = ""
        self.entries = {}
        self._reactions = {}
        self._relations = set()

    def get_KGML(self):
        """Return the pathway as a string in prettified KGML format."""
        header = "\n".join(
            [
                '<?xml version="1.0"?>',
                "<!DOCTYPE pathway SYSTEM "
                '"http://www.genome.jp/kegg/xml/'
                'KGML_v0.7.2_.dtd">',
                f"<!-- Created by KGML_Pathway.py {time.asctime()} -->",
            ]
        )
        rough_xml = header + ET.tostring(self.element, "utf-8").decode()
        reparsed = minidom.parseString(rough_xml)
        return reparsed.toprettyxml(indent="  ")

    def add_entry(self, entry):
        """Add an Entry element to the pathway."""
        # We insist that the node ID is an integer
        if not isinstance(entry.id, int):
            raise TypeError(
                f"Node ID must be an integer, got {type(entry.id)} ({entry.id})"
            )
        entry._pathway = self  # Let the entry know about the pathway
        self.entries[entry.id] = entry

    def remove_entry(self, entry):
        """Remove an Entry element from the pathway."""
        if not isinstance(entry.id, int):
            raise TypeError(
                f"Node ID must be an integer, got {type(entry.id)} ({entry.id})"
            )
        # We need to remove the entry from any other elements that may
        # contain it, which means removing those elements
        # TODO
        del self.entries[entry.id]

    def add_reaction(self, reaction):
        """Add a Reaction element to the pathway."""
        # We insist that the node ID is an integer and corresponds to an entry
        if not isinstance(reaction.id, int):
            raise ValueError(
                f"Node ID must be an integer, got {type(reaction.id)} ({reaction.id})"
            )
        if reaction.id not in self.entries:
            raise ValueError("Reaction ID %d has no corresponding entry" % reaction.id)
        reaction._pathway = self  # Let the reaction know about the pathway
        self._reactions[reaction.id] = reaction

    def remove_reaction(self, reaction):
        """Remove a Reaction element from the pathway."""
        if not isinstance(reaction.id, int):
            raise TypeError(
                f"Node ID must be an integer, got {type(reaction.id)} ({reaction.id})"
            )
        # We need to remove the reaction from any other elements that may
        # contain it, which means removing those elements
        # TODO
        del self._reactions[reaction.id]

    def add_relation(self, relation):
        """Add a Relation element to the pathway."""
        relation._pathway = self  # Let the reaction know about the pathway
        self._relations.add(relation)

    def remove_relation(self, relation):
        """Remove a Relation element from the pathway."""
        self._relations.remove(relation)

    def __str__(self):
        """Return a readable summary description string."""
        outstr = [
            f"Pathway: {self.title}",
            f"KEGG ID: {self.name}",
            f"Image file: {self.image}",
            f"Organism: {self.org}",
            "Entries: %d" % len(self.entries),
            "Entry types:",
        ]
        for t in ["ortholog", "enzyme", "reaction", "gene", "group", "compound", "map"]:
            etype = [e for e in self.entries.values() if e.type == t]
            if len(etype):
                outstr.append("\t%s: %d" % (t, len(etype)))
        return "\n".join(outstr) + "\n"

    # Assert correct formatting of the pathway name, and other attributes
    def _getname(self):
        return self._name

    def _setname(self, value):
        if not value.startswith("path:"):
            raise ValueError(f"Pathway name should begin with 'path:', got {value}")
        self._name = value

    def _delname(self):
        del self._name

    name = property(_getname, _setname, _delname, "The KEGGID for the pathway map.")

    def _getnumber(self):
        return self._number

    def _setnumber(self, value):
        self._number = int(value)

    def _delnumber(self):
        del self._number

    number = property(_getnumber, _setnumber, _delnumber, "The KEGG map number.")

    @property
    def compounds(self):
        """Get a list of entries of type compound."""
        return [e for e in self.entries.values() if e.type == "compound"]

    @property
    def maps(self):
        """Get a list of entries of type map."""
        return [e for e in self.entries.values() if e.type == "map"]

    @property
    def orthologs(self):
        """Get a list of entries of type ortholog."""
        return [e for e in self.entries.values() if e.type == "ortholog"]

    @property
    def genes(self):
        """Get a list of entries of type gene."""
        return [e for e in self.entries.values() if e.type == "gene"]

    @property
    def reactions(self):
        """Get a list of reactions in the pathway."""
        return self._reactions.values()

    @property
    def reaction_entries(self):
        """List of entries corresponding to each reaction in the pathway."""
        return [self.entries[i] for i in self._reactions]

    @property
    def relations(self):
        """Get a list of relations in the pathway."""
        return list(self._relations)

    @property
    def element(self):
        """Return the Pathway as a valid KGML element."""
        # The root is this Pathway element
        pathway = ET.Element("pathway")
        pathway.attrib = {
            "name": self._name,
            "org": self.org,
            "number": str(self._number),
            "title": self.title,
            "image": self.image,
            "link": self.link,
        }
        # We add the Entries in node ID order
        for eid, entry in sorted(self.entries.items()):
            pathway.append(entry.element)
        # Next we add Relations
        for relation in self._relations:
            pathway.append(relation.element)
        for eid, reaction in sorted(self._reactions.items()):
            pathway.append(reaction.element)
        return pathway

    @property
    def bounds(self):
        """Coordinate bounds for all Graphics elements in the Pathway.

        Returns the [(xmin, ymin), (xmax, ymax)] coordinates for all
        Graphics elements in the Pathway
        """
        xlist, ylist = [], []
        for b in [g.bounds for g in self.entries.values()]:
            xlist.extend([b[0][0], b[1][0]])
            ylist.extend([b[0][1], b[1][1]])
        return [(min(xlist), min(ylist)), (max(xlist), max(ylist))]


# Entry
class Entry:
    """Represent an Entry from KGML.

    Each Entry element is a node in the pathway graph, as described in
    release KGML v0.7.2 (http://www.kegg.jp/kegg/xml/docs/)

    Attributes:
     - id - The ID of the entry in the pathway map (integer)
     - names - List of KEGG IDs for the entry
     - type - The type of the entry
     - link - URL of information about the entry
     - reaction - List of KEGG IDs of the corresponding reactions
       (integer)
     - graphics -    List of Graphics objects describing the Entry's visual
       representation
     - components - List of component node ID for this Entry ('group')
     - alt - List of alternate names for the Entry

    NOTE: The alt attribute represents a subelement of the substrate and
    product elements in the KGML file

    """

    def __init__(self):
        """Initialize the class."""
        self._id = None
        self._names = []
        self.type = ""
        self.image = ""
        self.link = ""
        self.graphics = []
        self.components = set()
        self.alt = []
        self._pathway = None
        self._reactions = []

    def __str__(self):
        """Return readable descriptive string."""
        outstr = [
            "Entry node ID: %d" % self.id,
            f"Names: {self.name}",
            f"Type: {self.type}",
            f"Components: {self.components}",
            f"Reactions: {self.reaction}",
            "Graphics elements: %d %s" % (len(self.graphics), self.graphics),
        ]
        return "\n".join(outstr) + "\n"

    def add_component(self, element):
        """Add an element to the entry.

        If the Entry is already part of a pathway, make sure
        the component already exists.
        """
        if self._pathway is not None:
            if element.id not in self._pathway.entries:
                raise ValueError(
                    f"Component {element.id} is not an entry in the pathway"
                )
        self.components.add(element)

    def remove_component(self, value):
        """Remove the entry with the passed ID from the group."""
        self.components.remove(value)

    def add_graphics(self, entry):
        """Add the Graphics entry."""
        self.graphics.append(entry)

    def remove_graphics(self, entry):
        """Remove the Graphics entry with the passed ID from the group."""
        self.graphics.remove(entry)

    # Names may be given as a space-separated list of KEGG identifiers
    def _getname(self):
        return " ".join(self._names)

    def _setname(self, value):
        self._names = value.split()

    def _delname(self):
        self._names = []

    name = property(
        _getname, _setname, _delname, "List of KEGG identifiers for the Entry."
    )

    # Reactions may be given as a space-separated list of KEGG identifiers
    def _getreaction(self):
        return " ".join(self._reactions)

    def _setreaction(self, value):
        self._reactions = value.split()

    def _delreaction(self):
        self._reactions = []

    reaction = property(
        _getreaction,
        _setreaction,
        _delreaction,
        "List of reaction KEGG IDs for this Entry.",
    )

    # We make sure that the node ID is an integer
    def _getid(self):
        return self._id

    def _setid(self, value):
        self._id = int(value)

    def _delid(self):
        del self._id

    id = property(_getid, _setid, _delid, "The pathway graph node ID for the Entry.")

    @property
    def element(self):
        """Return the Entry as a valid KGML element."""
        # The root is this Entry element
        entry = ET.Element("entry")
        entry.attrib = {
            "id": str(self._id),
            "name": self.name,
            "link": self.link,
            "type": self.type,
        }
        if len(self._reactions):
            entry.attrib["reaction"] = self.reaction
        if len(self.graphics):
            for g in self.graphics:
                entry.append(g.element)
        if len(self.components):
            for c in self.components:
                entry.append(c.element)
        return entry

    @property
    def bounds(self):
        """Coordinate bounds for all Graphics elements in the Entry.

        Return the [(xmin, ymin), (xmax, ymax)] coordinates for the Entry
        Graphics elements.
        """
        xlist, ylist = [], []
        for b in [g.bounds for g in self.graphics]:
            xlist.extend([b[0][0], b[1][0]])
            ylist.extend([b[0][1], b[1][1]])
        return [(min(xlist), min(ylist)), (max(xlist), max(ylist))]

    @property
    def is_reactant(self):
        """Return true if this Entry participates in any reaction in its parent pathway."""
        for rxn in self._pathway.reactions:
            if self._id in rxn.reactant_ids:
                return True
        return False


# Component
class Component:
    """An Entry subelement used to represents a complex node.

    A subelement of the Entry element, used when the Entry is a complex
    node, as described in release KGML v0.7.2
    (http://www.kegg.jp/kegg/xml/docs/)

    The Component acts as a collection (with type 'group', and typically
    its own Graphics subelement), having only an ID.
    """

    def __init__(self, parent):
        """Initialize the class."""
        self._id = None
        self._parent = parent

    # We make sure that the node ID is an integer
    def _getid(self):
        return self._id

    def _setid(self, value):
        self._id = int(value)

    def _delid(self):
        del self._id

    id = property(_getid, _setid, _delid, "The pathway graph node ID for the Entry")

    @property
    def element(self):
        """Return the Component as a valid KGML element."""
        # The root is this Component element
        component = ET.Element("component")
        component.attrib = {"id": str(self._id)}
        return component


# Graphics
class Graphics:
    """An Entry subelement used to represents the visual representation.

    A subelement of Entry, specifying its visual representation, as
    described in release KGML v0.7.2 (http://www.kegg.jp/kegg/xml/docs/)

    Attributes:
     - name         Label for the graphics object
     - x            X-axis position of the object (int)
     - y            Y-axis position of the object (int)
     - coords       polyline coordinates, list of (int, int) tuples
     - type         object shape
     - width        object width (int)
     - height       object height (int)
     - fgcolor      object foreground color (hex RGB)
     - bgcolor      object background color (hex RGB)

    Some attributes are present only for specific graphics types.  For
    example, line types do not (typically) have a width.
    We permit non-DTD attributes and attribute settings, such as

    dash         List of ints, describing an on/off pattern for dashes

    """

    def __init__(self, parent):
        """Initialize the class."""
        self.name = ""
        self._x = None
        self._y = None
        self._coords = None
        self.type = ""
        self._width = None
        self._height = None
        self.fgcolor = ""
        self.bgcolor = ""
        self._parent = parent

    # We make sure that the XY coordinates, width and height are numbers
    def _getx(self):
        return self._x

    def _setx(self, value):
        self._x = float(value)

    def _delx(self):
        del self._x

    x = property(_getx, _setx, _delx, "The X coordinate for the graphics element.")

    def _gety(self):
        return self._y

    def _sety(self, value):
        self._y = float(value)

    def _dely(self):
        del self._y

    y = property(_gety, _sety, _dely, "The Y coordinate for the graphics element.")

    def _getwidth(self):
        return self._width

    def _setwidth(self, value):
        self._width = float(value)

    def _delwidth(self):
        del self._width

    width = property(
        _getwidth, _setwidth, _delwidth, "The width of the graphics element."
    )

    def _getheight(self):
        return self._height

    def _setheight(self, value):
        self._height = float(value)

    def _delheight(self):
        del self._height

    height = property(
        _getheight, _setheight, _delheight, "The height of the graphics element."
    )

    # We make sure that the polyline coordinates are integers, too
    def _getcoords(self):
        return self._coords

    def _setcoords(self, value):
        clist = [int(e) for e in value.split(",")]
        self._coords = [tuple(clist[i : i + 2]) for i in range(0, len(clist), 2)]

    def _delcoords(self):
        del self._coords

    coords = property(
        _getcoords,
        _setcoords,
        _delcoords,
        "Polyline coordinates for the graphics element.",
    )

    # Set default colors
    def _getfgcolor(self):
        return self._fgcolor

    def _setfgcolor(self, value):
        if value == "none":
            self._fgcolor = "#000000"  # this default defined in KGML spec
        else:
            self._fgcolor = value

    def _delfgcolor(self):
        del self._fgcolor

    fgcolor = property(_getfgcolor, _setfgcolor, _delfgcolor, "Foreground color.")

    def _getbgcolor(self):
        return self._bgcolor

    def _setbgcolor(self, value):
        if value == "none":
            self._bgcolor = "#000000"  # this default defined in KGML spec
        else:
            self._bgcolor = value

    def _delbgcolor(self):
        del self._bgcolor

    bgcolor = property(_getbgcolor, _setbgcolor, _delbgcolor, "Background color.")

    @property
    def element(self):
        """Return the Graphics as a valid KGML element."""
        # The root is this Component element
        graphics = ET.Element("graphics")
        if isinstance(self.fgcolor, str):  # Assumes that string is hexstring
            fghex = self.fgcolor
        else:  # Assumes ReportLab Color object
            fghex = "#" + self.fgcolor.hexval()[2:]
        if isinstance(self.bgcolor, str):  # Assumes that string is hexstring
            bghex = self.bgcolor
        else:  # Assumes ReportLab Color object
            bghex = "#" + self.bgcolor.hexval()[2:]
        graphics.attrib = {
            "name": self.name,
            "type": self.type,
            "fgcolor": fghex,
            "bgcolor": bghex,
        }
        for (n, attr) in [
            ("x", "_x"),
            ("y", "_y"),
            ("width", "_width"),
            ("height", "_height"),
        ]:
            if getattr(self, attr) is not None:
                graphics.attrib[n] = str(getattr(self, attr))
        if self.type == "line":  # Need to write polycoords
            graphics.attrib["coords"] = ",".join(
                [str(e) for e in chain.from_iterable(self.coords)]
            )
        return graphics

    @property
    def bounds(self):
        """Coordinate bounds for the Graphics element.

        Return the bounds of the Graphics object as an [(xmin, ymin),
        (xmax, ymax)] tuple.  Coordinates give the centre of the
        circle, rectangle, roundrectangle elements, so we have to
        adjust for the relevant width/height.
        """
        if self.type == "line":
            xlist = [x for x, y in self.coords]
            ylist = [y for x, y in self.coords]
            return [(min(xlist), min(ylist)), (max(xlist), max(ylist))]
        else:
            return [
                (self.x - self.width * 0.5, self.y - self.height * 0.5),
                (self.x + self.width * 0.5, self.y + self.height * 0.5),
            ]

    @property
    def centre(self):
        """Return the centre of the Graphics object as an (x, y) tuple."""
        return (
            0.5 * (self.bounds[0][0] + self.bounds[1][0]),
            0.5 * (self.bounds[0][1] + self.bounds[1][1]),
        )


# Reaction
class Reaction:
    """A specific chemical reaction with substrates and products.

    This describes a specific chemical reaction between one or more
    substrates and one or more products.

    Attributes:
     - id             Pathway graph node ID of the entry
     - names          List of KEGG identifier(s) from the REACTION database
     - type           String: reversible or irreversible
     - substrate      Entry object of the substrate
     - product        Entry object of the product

    """

    def __init__(self):
        """Initialize the class."""
        self._id = None
        self._names = []
        self.type = ""
        self._substrates = set()
        self._products = set()
        self._pathway = None

    def __str__(self):
        """Return an informative human-readable string."""
        outstr = [
            f"Reaction node ID: {self.id}",
            f"Reaction KEGG IDs: {self.name}",
            f"Type: {self.type}",
            f"Substrates: {','.join([s.name for s in self.substrates])}",
            f"Products: {','.join([s.name for s in self.products])}",
        ]
        return "\n".join(outstr) + "\n"

    def add_substrate(self, substrate_id):
        """Add a substrate, identified by its node ID, to the reaction."""
        if self._pathway is not None:
            if int(substrate_id) not in self._pathway.entries:
                raise ValueError(
                    "Couldn't add substrate, no node ID %d in Pathway"
                    % int(substrate_id)
                )
        self._substrates.add(substrate_id)

    def add_product(self, product_id):
        """Add a product, identified by its node ID, to the reaction."""
        if self._pathway is not None:
            if int(product_id) not in self._pathway.entries:
                raise ValueError(
                    "Couldn't add product, no node ID %d in Pathway" % product_id
                )
        self._products.add(int(product_id))

    # The node ID is also the node ID of the Entry that corresponds to the
    # reaction; we get the corresponding Entry when there is an associated
    # Pathway
    def _getid(self):
        return self._id

    def _setid(self, value):
        self._id = int(value)

    def _delid(self):
        del self._id

    id = property(_getid, _setid, _delid, "Node ID for the reaction.")

    # Names may show up as a space-separated list of several KEGG identifiers
    def _getnames(self):
        return " ".join(self._names)

    def _setnames(self, value):
        self._names.extend(value.split())

    def _delnames(self):
        del self.names

    name = property(
        _getnames, _setnames, _delnames, "List of KEGG identifiers for the reaction."
    )

    # products and substrates are read-only properties, returning lists
    # of Entry objects
    @property
    def substrates(self):
        """Return list of substrate Entry elements."""
        return [self._pathway.entries[sid] for sid in self._substrates]

    @property
    def products(self):
        """Return list of product Entry elements."""
        return [self._pathway.entries[pid] for pid in self._products]

    @property
    def entry(self):
        """Return the Entry corresponding to this reaction."""
        return self._pathway.entries[self._id]

    @property
    def reactant_ids(self):
        """Return a list of substrate and product reactant IDs."""
        return self._products.union(self._substrates)

    @property
    def element(self):
        """Return KGML element describing the Reaction."""
        # The root is this Relation element
        reaction = ET.Element("reaction")
        reaction.attrib = {"id": str(self.id), "name": self.name, "type": self.type}
        for s in self._substrates:
            substrate = ET.Element("substrate")
            substrate.attrib["id"] = str(s)
            substrate.attrib["name"] = self._pathway.entries[s].name
            reaction.append(substrate)
        for p in self._products:
            product = ET.Element("product")
            product.attrib["id"] = str(p)
            product.attrib["name"] = self._pathway.entries[p].name
            reaction.append(product)
        return reaction


# Relation
class Relation:
    """A relationship between to products, KOs, or protein and compound.

    This describes a relationship between two products, KOs, or protein
    and compound, as described in release KGML v0.7.2
    (http://www.kegg.jp/kegg/xml/docs/)

    Attributes:
     - entry1 - The first Entry object node ID defining the
       relation (int)
     - entry2 - The second Entry object node ID defining the
       relation (int)
     - type - The relation type
     - subtypes - List of subtypes for the relation, as a list of
       (name, value) tuples

    """

    def __init__(self):
        """Initialize the class."""
        self._entry1 = None
        self._entry2 = None
        self.type = ""
        self.subtypes = []
        self._pathway = None

    def __str__(self):
        """Return a useful human-readable string."""
        outstr = [
            "Relation (subtypes: %d):" % len(self.subtypes),
            "Entry1:",
            str(self.entry1),
            "Entry2:",
            str(self.entry2),
        ]
        for s in self.subtypes:
            outstr.extend([f"Subtype: {s[0]}", str(s[1])])
        return "\n".join(outstr)

    # Properties entry1 and entry2
    def _getentry1(self):
        if self._pathway is not None:
            return self._pathway.entries[self._entry1]
        return self._entry1

    def _setentry1(self, value):
        self._entry1 = int(value)

    def _delentry1(self):
        del self._entry1

    entry1 = property(_getentry1, _setentry1, _delentry1, "Entry1 of the relation.")

    def _getentry2(self):
        if self._pathway is not None:
            return self._pathway.entries[self._entry2]
        return self._entry2

    def _setentry2(self, value):
        self._entry2 = int(value)

    def _delentry2(self):
        del self._entry2

    entry2 = property(_getentry2, _setentry2, _delentry2, "Entry2 of the relation.")

    @property
    def element(self):
        """Return KGML element describing the Relation."""
        # The root is this Relation element
        relation = ET.Element("relation")
        relation.attrib = {
            "entry1": str(self._entry1),
            "entry2": str(self._entry2),
            "type": self.type,
        }
        for (name, value) in self.subtypes:
            subtype = ET.Element("subtype")
            subtype.attrib = {"name": name, "value": str(value)}
            relation.append(subtype)
        return relation
