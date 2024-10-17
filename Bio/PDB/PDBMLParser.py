# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""This module contains a parser for PDBML (PDB XML) files.

PDBML is a representation of PDB data in XML format.

See https://pdbml.wwpdb.org/.
"""

from os import PathLike
from typing import TextIO
from typing import Union
from xml.etree import ElementTree
from xml.etree.ElementTree import Element

import numpy as np

from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder


def _parse_resolution_from(
    tree: ElementTree, namespaces: dict[str, str]
) -> Union[float, None]:
    for candidate in [
        "PDBx:refineCategory/PDBx:refine/PDBx:ls_d_res_high",
        "PDBx:refine_histCategory/PDBx:refine_hist/PDBx:d_res_high",
        "PDBx:em_3d_reconstructionCategory/PDBx:em_3d_reconstruction/PDBx:resolution",
    ]:
        element = tree.find(candidate, namespaces)
        if element:
            text = element.text
            if text:
                return float(text)
    return None


def _parse_header_from(
    tree: ElementTree, namespaces: dict[str, str]
) -> dict[str, Union[str, float]]:
    return {
        "name": tree.find(
            "PDBx:structCategory/PDBx:struct/PDBx:title", namespaces
        ).text,
        "head": tree.find(
            "PDBx:struct_keywordsCategory/PDBx:struct_keywords/PDBx:text", namespaces
        ).text,
        "idcode": tree.find(
            "PDBx:entryCategory/PDBx:entry",
            namespaces,
        ).attrib["id"],
        "deposition_date": tree.find(
            "PDBx:pdbx_database_statusCategory/PDBx:pdbx_database_status/PDBx:recvd_initial_deposition_date",
            namespaces,
        ).text,
        "structure_method": tree.find(
            "PDBx:exptlCategory/PDBx:exptl", namespaces
        ).attrib["method"],
        "resolution": _parse_resolution_from(tree, namespaces),
    }


def _parse_atom_from(element: Element, namespaces: dict[str, str]):
    name = element.find("PDBx:label_atom_id", namespaces).text
    x = float(element.find("PDBx:Cartn_x", namespaces).text)
    y = float(element.find("PDBx:Cartn_y", namespaces).text)
    z = float(element.find("PDBx:Cartn_z", namespaces).text)

    return {
        "name": name,
        "coord": np.array((x, y, z), float),
        "b_factor": float(element.find("PDBx:B_iso_or_equiv", namespaces).text),
        "occupancy": float(element.find("PDBx:occupancy", namespaces).text),
        "altloc": element.find("PDBx:label_alt_id", namespaces).text or " ",
        "fullname": name,
        "serial_number": int(element.attrib["id"]),
        "element": element.find("PDBx:type_symbol", namespaces).text,
    }


def _parse_residue_id_from(
    element: Element, namespaces: dict[str, str]
) -> tuple[str, int, str]:
    assert element.tag == f"{{{namespaces['PDBx']}}}atom_site"
    atom_group = element.find("PDBx:group_PDB", namespaces).text
    component_id = element.find("PDBx:label_comp_id", namespaces).text
    ins_code_element = element.find("PDBx:pdbx_PDB_ins_code", namespaces)

    if atom_group == "HETATM":
        if component_id == "HOH" or component_id == "WAT":
            hetero_field = "W"
        else:
            hetero_field = "H"
    else:
        hetero_field = " "

    sequence_id = int(element.find("PDBx:auth_seq_id", namespaces).text)
    insertion_code = ins_code_element.text if ins_code_element is not None else " "
    return hetero_field, sequence_id, insertion_code


class PDBMLParser:
    """A parser for PDBML (PDB XML) files. See https://pdbml.wwpdb.org/.

    This parser is based on the mmCIF parser also provided in the PDB package in the sense that the structure object
    returned by this parser is equal to the structure returned by the mmCIF parser for any given PDB structure.
    """

    def __init__(self):
        """Initialize a PDBML parser."""
        self.structure_builder = StructureBuilder()

    def get_structure(
        self, source: Union[int, str, bytes, PathLike, TextIO]
    ) -> Structure:
        """Parse and return the PDB structure from XML source.

        :param Union[int, str, bytes, PathLike, TextIO] source: The XML representation of the PDB structure
        :return: the PDB structure
        :rtype: Bio.PDB.Structure.Structure
        """
        namespaces = dict(
            [node for _, node in ElementTree.iterparse(source, ["start-ns"])]
        )

        if hasattr(source, "seek"):
            # This resets the source if source is a file handle.
            source.seek(0)

        tree = ElementTree.parse(source)
        header = _parse_header_from(tree, namespaces)
        self.structure_builder.init_structure(structure_id=header["idcode"])
        self.structure_builder.set_header(header)
        atom_elements = tree.find("PDBx:atom_siteCategory", namespaces)

        builder_model_count = 0
        builder_model_number = None
        builder_chain_id = None
        builder_residue_id = None
        builder_component_id = None
        for element in atom_elements:
            # See https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/atom_site.html
            # for definitions of the different elements.
            model_number = int(element.find("PDBx:pdbx_PDB_model_num", namespaces).text)
            chain_id = element.find("PDBx:auth_asym_id", namespaces).text
            component_id = element.find("PDBx:label_comp_id", namespaces).text
            residue_id = _parse_residue_id_from(element, namespaces)

            if not builder_model_number or model_number != builder_model_number:
                # Ideally, this parser would initialize the model with the model number as the ID.
                # However, the mmCIF parser uses the model's index in the structure (starting from zero) as the ID.
                # We want the structure returned by this parser to be the same as that of the mmCIF parser.
                self.structure_builder.init_model(builder_model_count, model_number)
                builder_model_count += 1
                builder_model_number = model_number
                builder_chain_id = None
                builder_residue_id = None
            if not builder_chain_id or chain_id != builder_chain_id:
                self.structure_builder.init_chain(chain_id)
                builder_chain_id = chain_id
                builder_residue_id = None
            if residue_id != builder_residue_id or component_id != builder_component_id:
                self.structure_builder.init_residue(component_id, *residue_id)
                builder_residue_id = residue_id
                builder_component_id = component_id

            self.structure_builder.init_atom(**_parse_atom_from(element, namespaces))

        return self.structure_builder.get_structure()
