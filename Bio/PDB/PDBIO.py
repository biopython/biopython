# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Output of PDB files."""

from Bio.Data.IUPACData import atom_weights # Allowed Elements

_ATOM_FORMAT_STRING="%s%5i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"


class Select(object):
    """
    Default selection (everything) during writing - can be used as base class
    to implement selective output. This selects which entities will be written out.
    """

    def __repr__(self):
        return "<Select all>"

    def accept_model(self, model):
        """
        Overload this to reject models for output.
        """
        return 1

    def accept_chain(self, chain):
        """
        Overload this to reject chains for output.
        """
        return 1

    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        return 1

    def accept_atom(self, atom):
        """
        Overload this to reject atoms for output.
        """
        return 1


class PDBIO(object):
    """
    Write a Structure object (or a subset of a Structure object) as a PDB file.


    Example:
        >>> p=PDBParser()
        >>> s=p.get_structure("1fat", "1fat.pdb")
        >>> io=PDBIO()
        >>> io.set_structure(s)
        >>> io.save("out.pdb")
    """
    def __init__(self, use_model_flag=0):
        """
        @param use_model_flag: if 1, force use of the MODEL record in output.
        @type use_model_flag: int
        """
        self.use_model_flag=use_model_flag
    
    # private mathods

    def _get_atom_line(self, atom, hetfield, segid, atom_number, resname, 
        resseq, icode, chain_id, charge="  "):
        """Returns an ATOM PDB string (PRIVATE)."""
        if hetfield!=" ":
            record_type="HETATM"
        else:
            record_type="ATOM  "
        if atom.element:
            element = atom.element.strip().upper()
            if element.capitalize() not in atom_weights:
                raise ValueError("Unrecognised element %r" % atom.element)
            element = element.rjust(2)
        else:
            element = "  "
        name=atom.get_fullname()
        altloc=atom.get_altloc()
        x, y, z=atom.get_coord()
        bfactor=atom.get_bfactor()
        occupancy=atom.get_occupancy()
        args=(record_type, atom_number, name, altloc, resname, chain_id,
            resseq, icode, x, y, z, occupancy, bfactor, segid,
            element, charge)
        return _ATOM_FORMAT_STRING % args

    # Public methods

    def set_structure(self, structure):
        self.structure=structure

    def save(self, file, select=Select(), write_end=0):
        """
        @param file: output file
        @type file: string or filehandle 

        @param select: selects which entities will be written.
        @type select: 
            select hould have the following methods:
                - accept_model(model)
                - accept_chain(chain)
                - accept_residue(residue)
                - accept_atom(atom)
            These methods should return 1 if the entity
            is to be written out, 0 otherwise.

            Typically select is a subclass of L{Select}.
        """
        get_atom_line=self._get_atom_line
        if isinstance(file, basestring):
            fp=open(file, "w")
            close_file=1
        else:
            # filehandle, I hope :-)
            fp=file
            close_file=0
        # multiple models?
        if len(self.structure)>1 or self.use_model_flag:
            model_flag=1
        else:
            model_flag=0
        for model in self.structure.get_list():
            if not select.accept_model(model):
                continue
            # necessary for ENDMDL 
            # do not write ENDMDL if no residues were written
            # for this model
            model_residues_written=0
            atom_number=1
            if model_flag:
                fp.write("MODEL      %s\n" % model.serial_num)
            for chain in model.get_list():
                if not select.accept_chain(chain):
                    continue
                chain_id=chain.get_id()
                # necessary for TER 
                # do not write TER if no residues were written
                # for this chain
                chain_residues_written=0
                for residue in chain.get_unpacked_list():
                    if not select.accept_residue(residue):
                        continue
                    hetfield, resseq, icode=residue.get_id()
                    resname=residue.get_resname()  
                    segid=residue.get_segid()
                    for atom in residue.get_unpacked_list():
                        if select.accept_atom(atom):
                            chain_residues_written=1
                            model_residues_written=1
                            s=get_atom_line(atom, hetfield, segid, atom_number, resname,
                                resseq, icode, chain_id)
                            fp.write(s)
                            atom_number=atom_number+1
                if chain_residues_written:
                    fp.write("TER\n")
            if model_flag and model_residues_written:
                fp.write("ENDMDL\n")
            if write_end:
                fp.write('END\n')
        if close_file:
            fp.close()

if __name__=="__main__":
    
    from Bio.PDB.PDBParser import PDBParser

    import sys

    p=PDBParser(PERMISSIVE=True)

    s=p.get_structure("test", sys.argv[1])

    io=PDBIO()
    io.set_structure(s)
    io.save("out1.pdb")

    fp=open("out2.pdb", "w")
    s1=p.get_structure("test1", sys.argv[1])
    s2=p.get_structure("test2", sys.argv[2])
    io=PDBIO(1)
    io.set_structure(s1)
    io.save(fp)
    io.set_structure(s2)
    io.save(fp, write_end=1)
    fp.close()




