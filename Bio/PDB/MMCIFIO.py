"""Write an mmCIF file."""

import re

from Bio.PDB.StructureBuilder import StructureBuilder

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
        self.dic = dic

    def save(self, filepath):
        # Decide whether to save a Structure object or an mmCIF dictionary
        if hasattr(self, 'structure'):
            self._save_structure(filepath)
        else:
            self._save_dict(filepath)

    def _save_dict(self, filepath):
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

        out_file = open(filepath, "w")
        if data_val:
            out_file.write("data_"+data_val+"\n")

        for key, key_list in key_lists.items():
            out_file.write("#\n")
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
                for i in key_list:
                    if len(i) > m:
                        m = len(i)
                for i in key_list:
                    out_file.write("{k: <{width}}".format(k=key+"."+i, width=len(key)+m+4)+ self._format_mmcif_col(self.dic[key+"."+i], len(self.dic[key+"."+i])+2)+"\n")
            # If the value is a list, write as keys then a value table
            elif val_type == list:
                out_file.write("loop_\n")
                col_widths = {}
                # Write keys and find max widths for each set of values
                for i in key_list:
                    out_file.write(key+"."+i+"\n")
                    col_widths[i] = 0
                    for val in self.dic[key+"."+i]:
                        if len(val) > col_widths[i]:
                            col_widths[i] = len(val)
                # Write the values as rows
                for i in range(n_vals):
                    for col in key_list:
                        out_file.write(self._format_mmcif_col(self.dic[key+"."+col][i], col_widths[col]+3))
                    out_file.write("\n")
            else:
                raise(ValueError("Invalid type in mmCIF dictionary: "+str(val_type)))
        out_file.close()

    def _format_mmcif_col(self, val, col_width):
        # Format a mmCIF data value by enclosing with quotes or semicolon lines
        # where appropriate. See
        # https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax for
        # syntax.
        # Actually any whitespace after quote and below, make this a regex?
        # If there is a newline or quotes cannot be contained, use semicolon construct
        if "\n" in val or ("' " in val and "\" " in val):
            return "\n;"+val+"\n;\n"
        # If quoting is required
        elif " " in val or "'" in val or "\"" in val or val[0] in ["_", "#", "$", "'", "\"", "[", "]", ";"] or val.startswith("data_") or val.startswith("save_") or val in ["loop_", "stop_", "global_"]:
            # Choose quote character
            if "' " in val:
                return "{v: <{width}}".format(v="\""+val+"\"", width=col_width)
            else:
                return "{v: <{width}}".format(v="'"+val+"'", width=col_width)
        # Safe to not quote
        else:
            return "{v: <{width}}".format(v=val, width=col_width)

    def _save_structure(self, filepath):
        pass
