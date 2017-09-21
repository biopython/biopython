"""Write an mmCIF file.
See https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for syntax.
"""

import re
from collections import defaultdict

from Bio._py3k import basestring
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.PDBIO import Select

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

class MMCIFIO(object):
    """Write a Structure object, a subset of a Structure object or a mmCIF
    dictionary as a mmCIF file.

    Example:
        >>> p=MMCIFParser()
        >>> s=p.get_structure("1fat", "1fat.cif")
        >>> io=MMCIFIO()
        >>> io.set_structure(s)
        >>> io.save("out.cif")

    """

    def __init__(self):
        pass

    def set_structure(self, pdb_object):
        # Check what the user is providing and build a structure appropriately
        # This is duplicated from the PDBIO class
        if pdb_object.level == "S":
            structure = pdb_object
        else:
            sb = StructureBuilder()
            sb.init_structure('pdb')
            sb.init_seg(' ')
            # Build parts as necessary
            if pdb_object.level == "M":
                sb.structure.add(pdb_object)
                self.structure = sb.structure
            else:
                sb.init_model(0)
                if pdb_object.level == "C":
                    sb.structure[0].add(pdb_object)
                else:
                    sb.init_chain('A')
                    if pdb_object.level == "R":
                        try:
                            parent_id = pdb_object.parent.id
                            sb.structure[0]['A'].id = parent_id
                        except Exception:
                            pass
                        sb.structure[0]['A'].add(pdb_object)
                    else:
                        # Atom
                        sb.init_residue('DUM', ' ', 1, ' ')
                        try:
                            parent_id = pdb_object.parent.parent.id
                            sb.structure[0]['A'].id = parent_id
                        except Exception:
                            pass
                        sb.structure[0]['A'].child_list[0].add(pdb_object)

            # Return structure
            structure = sb.structure
        self.structure = structure

    def set_dict(self, dic):
        # Set the mmCIF dictionary to be written out
        self.dic = dic

    def save(self, filepath, select=Select(), preserve_atom_numbering=False):
        if isinstance(filepath, basestring):
            fp = open(filepath, "w")
            close_file = True
        else:
            fp = filepath
            close_file = False
        # Decide whether to save a Structure object or an mmCIF dictionary
        if hasattr(self, 'structure'):
            self._save_structure(fp, select, preserve_atom_numbering)
        elif hasattr(self, 'dic'):
            self._save_dict(fp)
        else:
            raise(ValueError("Use set_structure or set_dict to set a structure or dictionary to write out"))
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
                    raise(ValueError("Invalid key in mmCIF dictionary: "+key))

        # Re-order lists if an order has been specified
        # Not all elements from the specified order are necessarily present
        for key, key_list in key_lists.items():
            if key in mmcif_order:
                inds = [mmcif_order[key].index(i) for i in key_list]
                z = zip(inds, key_list)
                z.sort()
                key_lists[key] = [k for _,k in z]

        # Write out top data_ line
        if data_val:
            out_file.write("data_"+data_val+"\n")
            out_file.write("#\n")

        for key, key_list in key_lists.items():
            # Pick a sample mmCIF value, which can be a list or a single value
            sample_val = self.dic[key+"."+key_list[0]]
            val_type = type(sample_val)
            n_vals = len(sample_val)
            # Check the mmCIF dictionary has consistent list sizes
            for i in key_list:
                val = self.dic[key+"."+i]
                if type(val) != val_type or (val_type == list and len(val) != n_vals):
                    raise(ValueError("Inconsistent list sizes in mmCIF dictionary: "+key+"."+i))
            # If the value is a single value, write as key-value pairs
            if val_type == str:
                m = 0
                # Find the maximum key length
                for i in key_list:
                    if len(i) > m:
                        m = len(i)
                for i in key_list:
                    out_file.write("{k: <{width}}".format(k=key+"."+i, width=len(key)+m+4) + self._format_mmcif_col(self.dic[key+"."+i], len(self.dic[key+"."+i]))+"\n")
            # If the value is a list, write as keys then a value table
            elif val_type == list:
                out_file.write("loop_\n")
                col_widths = {}
                # Write keys and find max widths for each set of values
                for i in key_list:
                    out_file.write(key+"."+i+"\n")
                    col_widths[i] = 0
                    for val in self.dic[key+"."+i]:
                        l = len(val)
                        # If the value requires quoting it will add 2 characters
                        if self._requires_quote(val):
                            l += 2
                        if l > col_widths[i]:
                            col_widths[i] = l
                # Technically the max of the sum of the column widths is 2048

                # Write the values as rows
                for i in range(n_vals):
                    for col in key_list:
                        out_file.write(self._format_mmcif_col(self.dic[key+"."+col][i], col_widths[col]+1))
                    out_file.write("\n")
            else:
                raise(ValueError("Invalid type in mmCIF dictionary: "+str(val_type)))
            out_file.write("#\n")

    def _format_mmcif_col(self, val, col_width):
        # Format a mmCIF data value by enclosing with quotes or semicolon lines
        # where appropriate. See
        # https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for
        # syntax.

        # If there is a newline or quotes cannot be contained, use semicolon
        # and newline construct
        if self._requires_newline(val):
            return "\n;"+val+"\n;\n"
        # Technically these should be case-insensitive
        elif self._requires_quote(val):
            # Choose quote character
            if "' " in val:
                return "{v: <{width}}".format(v="\""+val+"\"", width=col_width)
            else:
                return "{v: <{width}}".format(v="'"+val+"'", width=col_width)
        # Safe to not quote
        # Numbers must not be quoted
        else:
            return "{v: <{width}}".format(v=val, width=col_width)

    def _requires_newline(self, val):
        # Technically the space can be a tab too
        if "\n" in val or ("' " in val and "\" " in val):
            return True
        else:
            return False

    def _requires_quote(self, val):
        if " " in val or "'" in val or "\"" in val or val[0] in ["_", "#", "$", "[", "]", ";"] or val.startswith("data_") or val.startswith("save_") or val in ["loop_", "stop_", "global_"]:
            return True
        else:
            return False

    def _save_structure(self, out_file, select, preserve_atom_numbering):
        atom_dict = defaultdict(list)

        for model in self.structure.get_list():
            if not select.accept_model(model):
                continue
            # mmCIF files with a single model have it specified as model 1
            if model.id == 0:
                model_n = "1"
            else:
                model_n = str(model.id)
            if not preserve_atom_numbering:
                atom_number = 1
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                chain_id = chain.get_id()
                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        continue
                    hetfield, resseq, icode = residue.get_id()
                    if hetfield == " ":
                        residue_type = "ATOM"
                    else:
                        residue_type = "HETATM"
                    resseq = str(resseq)
                    if icode == " ":
                        icode = "?"
                    resname = residue.get_resname()
                    for atom in residue.get_unpacked_list():
                        if select.accept_atom(atom):
                            atom_dict["_atom_site.group_PDB"].append(residue_type)
                            if preserve_atom_numbering:
                                atom_number = atom.get_serial_number()
                            atom_dict["_atom_site.id"].append(str(atom_number))
                            if not preserve_atom_numbering:
                                atom_number += 1
                            if atom.element:
                                atom_dict["_atom_site.type_symbol"].append(atom.element)
                            else:
                                atom_dict["_atom_site.type_symbol"].append("?")
                            atom_dict["_atom_site.label_atom_id"].append(atom.get_name().strip())
                            altloc = atom.get_altloc()
                            if altloc == " ":
                                atom_dict["_atom_site.label_alt_id"].append(".")
                            else:
                                atom_dict["_atom_site.label_alt_id"].append(altloc)
                            atom_dict["_atom_site.label_comp_id"].append(resname)
                            atom_dict["_atom_site.pdbx_PDB_ins_code"].append(icode)
                            coord = atom.get_coord()
                            atom_dict["_atom_site.Cartn_x"].append("%.3f" % coord[0])
                            atom_dict["_atom_site.Cartn_y"].append("%.3f" % coord[1])
                            atom_dict["_atom_site.Cartn_z"].append("%.3f" % coord[2])
                            atom_dict["_atom_site.occupancy"].append(str(atom.get_occupancy()))
                            atom_dict["_atom_site.B_iso_or_equiv"].append(str(atom.get_bfactor()))
                            atom_dict["_atom_site.auth_seq_id"].append(resseq)
                            atom_dict["_atom_site.auth_asym_id"].append(chain_id)
                            atom_dict["_atom_site.pdbx_PDB_model_num"].append(model_n)

        atom_dict["data_"] = self.structure.id
        self.dic = atom_dict
        self._save_dict(out_file)
