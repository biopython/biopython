# Copyright 2017 Joe Greener. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Write an mmCIF file.

See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for syntax.
"""

import re
from collections import defaultdict

from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBIO import Select, StructureIO

# If certain entries should have a certain order of keys, that is specified here
mmcif_order = {
    "_atom_site": [
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_formal_charge",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    ]
}


_select = Select()


class MMCIFIO(StructureIO):
    """Write a Structure object or a mmCIF dictionary as a mmCIF file.

    Examples
    --------
        >>> from Bio.PDB import MMCIFParser
        >>> from Bio.PDB.mmcifio import MMCIFIO
        >>> parser = MMCIFParser()
        >>> structure = parser.get_structure("1a8o", "PDB/1A8O.cif")
        >>> io=MMCIFIO()
        >>> io.set_structure(structure)
        >>> io.save("bio-pdb-mmcifio-out.cif")
        >>> import os
        >>> os.remove("bio-pdb-mmcifio-out.cif")  # tidy up


    """

    def __init__(self):
        """Initialise."""
        pass

    def set_dict(self, dic):
        """Set the mmCIF dictionary to be written out."""
        self.dic = dic
        # Remove self.structure if it has been set
        if hasattr(self, "structure"):
            delattr(self, "structure")

    def save(self, filepath, select=_select, preserve_atom_numbering=False):
        """Save the structure to a file.

        :param filepath: output file
        :type filepath: string or filehandle

        :param select: selects which entities will be written.
        :type select: object

        Typically select is a subclass of L{Select}, it should
        have the following methods:

         - accept_model(model)
         - accept_chain(chain)
         - accept_residue(residue)
         - accept_atom(atom)

        These methods should return 1 if the entity is to be
        written out, 0 otherwise.
        """
        # Similar to the PDBIO save method, we check if the filepath is a
        # string for a filepath or an open file handle
        if isinstance(filepath, str):
            fp = open(filepath, "w")
            close_file = True
        else:
            fp = filepath
            close_file = False
        # Decide whether to save a Structure object or an mmCIF dictionary
        if hasattr(self, "structure"):
            self._save_structure(fp, select, preserve_atom_numbering)
        elif hasattr(self, "dic"):
            self._save_dict(fp)
        else:
            raise ValueError(
                "Use set_structure or set_dict to set a structure or dictionary to write out"
            )
        if close_file:
            fp.close()

    def _save_dict(self, out_file):
        # Form dictionary where key is first part of mmCIF key and value is list
        # of corresponding second parts
        key_lists = {}
        for key in self.dic:
            if key == "data_":
                data_val = self.dic[key]
            else:
                s = re.split(r"\.", key)
                if len(s) == 2:
                    if s[0] in key_lists:
                        key_lists[s[0]].append(s[1])
                    else:
                        key_lists[s[0]] = [s[1]]
                else:
                    raise ValueError("Invalid key in mmCIF dictionary: " + key)

        # Re-order lists if an order has been specified
        # Not all elements from the specified order are necessarily present
        for key, key_list in key_lists.items():
            if key in mmcif_order:
                inds = []
                for i in key_list:
                    try:
                        inds.append(mmcif_order[key].index(i))
                    # Unrecognised key - add at end
                    except ValueError:
                        inds.append(len(mmcif_order[key]))
                key_lists[key] = [k for _, k in sorted(zip(inds, key_list))]

        # Write out top data_ line
        if data_val:
            out_file.write("data_" + data_val + "\n#\n")

        for key, key_list in key_lists.items():
            # Pick a sample mmCIF value, which can be a list or a single value
            sample_val = self.dic[key + "." + key_list[0]]
            n_vals = len(sample_val)
            # Check the mmCIF dictionary has consistent list sizes
            for i in key_list:
                val = self.dic[key + "." + i]
                if (
                    isinstance(sample_val, list)
                    and (isinstance(val, str) or len(val) != n_vals)
                ) or (isinstance(sample_val, str) and isinstance(val, list)):
                    raise ValueError(
                        "Inconsistent list sizes in mmCIF dictionary: " + key + "." + i
                    )
            # If the value is a single value, write as key-value pairs
            if isinstance(sample_val, str) or (
                isinstance(sample_val, list) and len(sample_val) == 1
            ):
                m = 0
                # Find the maximum key length
                for i in key_list:
                    if len(i) > m:
                        m = len(i)
                for i in key_list:
                    # If the value is a single item list, just take the value
                    if isinstance(sample_val, str):
                        value_no_list = self.dic[key + "." + i]
                    else:
                        value_no_list = self.dic[key + "." + i][0]
                    out_file.write(
                        "{k: <{width}}".format(k=key + "." + i, width=len(key) + m + 4)
                        + self._format_mmcif_col(value_no_list, len(value_no_list))
                        + "\n"
                    )
            # If the value is more than one value, write as keys then a value table
            elif isinstance(sample_val, list):
                out_file.write("loop_\n")
                col_widths = {}
                # Write keys and find max widths for each set of values
                for i in key_list:
                    out_file.write(key + "." + i + "\n")
                    col_widths[i] = 0
                    for val in self.dic[key + "." + i]:
                        len_val = len(val)
                        # If the value requires quoting it will add 2 characters
                        if self._requires_quote(val) and not self._requires_newline(
                            val
                        ):
                            len_val += 2
                        if len_val > col_widths[i]:
                            col_widths[i] = len_val
                # Technically the max of the sum of the column widths is 2048

                # Write the values as rows
                for i in range(n_vals):
                    for col in key_list:
                        out_file.write(
                            self._format_mmcif_col(
                                self.dic[key + "." + col][i], col_widths[col] + 1
                            )
                        )
                    out_file.write("\n")
            else:
                raise ValueError(
                    "Invalid type in mmCIF dictionary: " + str(type(sample_val))
                )
            out_file.write("#\n")

    def _format_mmcif_col(self, val, col_width):
        # Format a mmCIF data value by enclosing with quotes or semicolon lines
        # where appropriate. See
        # https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for
        # syntax.

        # If there is a newline or quotes cannot be contained, use semicolon
        # and newline construct
        if self._requires_newline(val):
            return "\n;" + val + "\n;\n"
        elif self._requires_quote(val):
            # Choose quote character
            if "' " in val:
                return "{v: <{width}}".format(v='"' + val + '"', width=col_width)
            else:
                return "{v: <{width}}".format(v="'" + val + "'", width=col_width)
        # Safe to not quote
        # Numbers must not be quoted
        else:
            return "{v: <{width}}".format(v=val, width=col_width)

    def _requires_newline(self, val):
        # Technically the space can be a tab too
        if "\n" in val or ("' " in val and '" ' in val):
            return True
        else:
            return False

    def _requires_quote(self, val):
        # Technically the words should be case-insensitive
        if (
            " " in val
            or "'" in val
            or '"' in val
            or val[0] in ["_", "#", "$", "[", "]", ";"]
            or val.startswith("data_")
            or val.startswith("save_")
            or val in ["loop_", "stop_", "global_"]
        ):
            return True
        else:
            return False

    def _get_label_asym_id(self, entity_id):
        # Convert a positive integer into a chain ID
        # Goes A to Z, then AA to ZA, AB to ZB etc
        # This is in line with existing mmCIF files
        div = entity_id
        out = ""
        while div > 0:
            mod = (div - 1) % 26
            out += chr(65 + mod)
            div = int((div - mod) / 26)
        return out

    def _save_structure(self, out_file, select, preserve_atom_numbering):
        atom_dict = defaultdict(list)

        for model in self.structure.get_list():
            if not select.accept_model(model):
                continue
            # mmCIF files with a single model have it specified as model 1
            if model.serial_num == 0:
                model_n = "1"
            else:
                model_n = str(model.serial_num)
            # This is used to write label_entity_id and label_asym_id and
            # increments from 1, changing with each molecule
            entity_id = 0
            if not preserve_atom_numbering:
                atom_number = 1
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                chain_id = chain.get_id()
                if chain_id == " ":
                    chain_id = "."
                # This is used to write label_seq_id and increments from 1,
                # remaining blank for hetero residues
                residue_number = 1
                prev_residue_type = ""
                prev_resname = ""
                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        continue
                    hetfield, resseq, icode = residue.get_id()
                    if hetfield == " ":
                        residue_type = "ATOM"
                        label_seq_id = str(residue_number)
                        residue_number += 1
                    else:
                        residue_type = "HETATM"
                        label_seq_id = "."
                    resseq = str(resseq)
                    if icode == " ":
                        icode = "?"
                    resname = residue.get_resname()
                    # Check if the molecule changes within the chain
                    # This will always increment for the first residue in a
                    # chain due to the starting values above
                    if residue_type != prev_residue_type or (
                        residue_type == "HETATM" and resname != prev_resname
                    ):
                        entity_id += 1
                    prev_residue_type = residue_type
                    prev_resname = resname
                    label_asym_id = self._get_label_asym_id(entity_id)
                    for atom in residue.get_unpacked_list():
                        if select.accept_atom(atom):
                            atom_dict["_atom_site.group_PDB"].append(residue_type)
                            if preserve_atom_numbering:
                                atom_number = atom.get_serial_number()
                            atom_dict["_atom_site.id"].append(str(atom_number))
                            if not preserve_atom_numbering:
                                atom_number += 1
                            element = atom.element.strip()
                            if element == "":
                                element = "?"
                            atom_dict["_atom_site.type_symbol"].append(element)
                            atom_dict["_atom_site.label_atom_id"].append(
                                atom.get_name().strip()
                            )
                            altloc = atom.get_altloc()
                            if altloc == " ":
                                altloc = "."
                            atom_dict["_atom_site.label_alt_id"].append(altloc)
                            atom_dict["_atom_site.label_comp_id"].append(
                                resname.strip()
                            )
                            atom_dict["_atom_site.label_asym_id"].append(label_asym_id)
                            # The entity ID should be the same for similar chains
                            # However this is non-trivial to calculate so we write "?"
                            atom_dict["_atom_site.label_entity_id"].append("?")
                            atom_dict["_atom_site.label_seq_id"].append(label_seq_id)
                            atom_dict["_atom_site.pdbx_PDB_ins_code"].append(icode)
                            coord = atom.get_coord()
                            atom_dict["_atom_site.Cartn_x"].append(f"{coord[0]:.3f}")
                            atom_dict["_atom_site.Cartn_y"].append(f"{coord[1]:.3f}")
                            atom_dict["_atom_site.Cartn_z"].append(f"{coord[2]:.3f}")
                            atom_dict["_atom_site.occupancy"].append(
                                str(atom.get_occupancy())
                            )
                            atom_dict["_atom_site.B_iso_or_equiv"].append(
                                str(atom.get_bfactor())
                            )
                            atom_dict["_atom_site.auth_seq_id"].append(resseq)
                            atom_dict["_atom_site.auth_asym_id"].append(chain_id)
                            atom_dict["_atom_site.pdbx_PDB_model_num"].append(model_n)

        # Data block name is the structure ID with special characters removed
        structure_id = self.structure.id
        for c in ["#", "$", "'", '"', "[", "]", " ", "\t", "\n"]:
            structure_id = structure_id.replace(c, "")
        atom_dict["data_"] = structure_id

        # Set the dictionary and write out using the generic dictionary method
        self.dic = atom_dict
        self._save_dict(out_file)
