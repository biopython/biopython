# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

"""Parser for PDB files."""

# For using with statement in Python 2.5 or Jython
from __future__ import with_statement

import warnings

import numpy

from Bio.File import as_handle

from Bio.PDB.PDBExceptions import \
        PDBConstructionException, PDBConstructionWarning
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.parse_pdb_header import _parse_pdb_header_list


# If PDB spec says "COLUMNS 18-20" this means line[17:20]


class PDBParser(object):
    """
    Parse a PDB file and return a Structure object.
    """

    def __init__(self, PERMISSIVE=True, get_header=False,
                 structure_builder=None, QUIET=False):
        """
        The PDB parser call a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instanciated by the
        PDBParser object itself, but if the user provides his own StructureBuilder
        object, the latter is used instead.

        Arguments:
        
        o PERMISSIVE - Evaluated as a Boolean. If false, exceptions in
        constructing the SMCRA data structure are fatal. If true (DEFAULT),
        the exceptions are caught, but some residues or atoms will be missing.
        THESE EXCEPTIONS ARE DUE TO PROBLEMS IN THE PDB FILE!.

        o structure_builder - an optional user implemented StructureBuilder class. 

        o QUIET - Evaluated as a Boolean. If true, warnings issued in constructing
        the SMCRA data will be supressed. If false (DEFAULT), they will be shown.
        These warnings might be indicative of problems in the PDB file!        
        """
        if structure_builder!=None:
            self.structure_builder=structure_builder
        else:
            self.structure_builder=StructureBuilder()
        self.header=None
        self.trailer=None
        self.line_counter=0
        self.PERMISSIVE=bool(PERMISSIVE)
        self.QUIET=bool(QUIET)

    # Public methods

    def get_structure(self, id, file):
        """Return the structure.

        Arguments:
        o id - string, the id that will be used for the structure
        o file - name of the PDB file OR an open filehandle
        """

        if self.QUIET:
            warning_list = warnings.filters[:]
            warnings.filterwarnings('ignore', category=PDBConstructionWarning)
            
        self.header=None
        self.trailer=None
        # Make a StructureBuilder instance (pass id of structure as parameter)
        self.structure_builder.init_structure(id)

        with as_handle(file) as handle:
            self._parse(handle.readlines())

        self.structure_builder.set_header(self.header)
        # Return the Structure instance
        structure = self.structure_builder.get_structure()
        
        if self.QUIET:
            warnings.filters = warning_list
        
        return structure

    def get_header(self):
        "Return the header."
        return self.header

    def get_trailer(self):
        "Return the trailer."
        return self.trailer

    # Private methods
    
    def _parse(self, header_coords_trailer):
        "Parse the PDB file."
        # Extract the header; return the rest of the file
        self.header, coords_trailer=self._get_header(header_coords_trailer)
        # Parse the atomic data; return the PDB file trailer
        self.trailer=self._parse_coordinates(coords_trailer)
    
    def _get_header(self, header_coords_trailer):
        "Get the header of the PDB file, return the rest."
        structure_builder=self.structure_builder
        i = 0
        for i in range(0, len(header_coords_trailer)):
            structure_builder.set_line_counter(i+1)
            line=header_coords_trailer[i]
            record_type=line[0:6] 
            if(record_type=='ATOM  ' or record_type=='HETATM' or record_type=='MODEL '):
                break
        header=header_coords_trailer[0:i]
        # Return the rest of the coords+trailer for further processing
        self.line_counter=i
        coords_trailer=header_coords_trailer[i:]
        header_dict=_parse_pdb_header_list(header)
        return header_dict, coords_trailer
    
    def _parse_coordinates(self, coords_trailer):
        "Parse the atomic data in the PDB file."
        local_line_counter=0
        structure_builder=self.structure_builder
        current_model_id=0
        # Flag we have an open model
        model_open=0
        current_chain_id=None
        current_segid=None
        current_residue_id=None
        current_resname=None
        for i in range(0, len(coords_trailer)):
            line=coords_trailer[i]
            record_type=line[0:6]
            global_line_counter=self.line_counter+local_line_counter+1
            structure_builder.set_line_counter(global_line_counter)
            if(record_type=='ATOM  ' or record_type=='HETATM'):
                # Initialize the Model - there was no explicit MODEL record
                if not model_open:
                    structure_builder.init_model(current_model_id)
                    current_model_id+=1
                    model_open=1
                fullname=line[12:16]
                # get rid of whitespace in atom names
                split_list=fullname.split()
                if len(split_list)!=1:
                    # atom name has internal spaces, e.g. " N B ", so
                    # we do not strip spaces
                    name=fullname
                else:
                    # atom name is like " CA ", so we can strip spaces
                    name=split_list[0]
                altloc=line[16:17]
                resname=line[17:20]
                chainid=line[21:22]
                try:
                    serial_number=int(line[6:11])
                except:
                    serial_number=0
                resseq=int(line[22:26].split()[0])   # sequence identifier   
                icode=line[26:27]           # insertion code
                if record_type=='HETATM':       # hetero atom flag
                    if resname=="HOH" or resname=="WAT":
                        hetero_flag="W"
                    else:
                        hetero_flag="H"
                else:
                    hetero_flag=" "
                residue_id=(hetero_flag, resseq, icode)
                # atomic coordinates
                try:
                    x=float(line[30:38]) 
                    y=float(line[38:46]) 
                    z=float(line[46:54])
                except:
                    #Should we allow parsing to continue in permissive mode?
                    #If so what coordindates should we default to?  Easier to abort!
                    raise PDBConstructionException(\
                        "Invalid or missing coordinate(s) at line %i." \
                        % global_line_counter)
                coord=numpy.array((x, y, z), 'f')
                # occupancy & B factor
                try:
                    occupancy=float(line[54:60])
                except:
                    self._handle_PDB_exception("Invalid or missing occupancy",
                                               global_line_counter)
                    occupancy = 0.0 #Is one or zero a good default?
                try:
                    bfactor=float(line[60:66])
                except:
                    self._handle_PDB_exception("Invalid or missing B factor",
                                               global_line_counter)
                    bfactor = 0.0 #The PDB use a default of zero if the data is missing
                segid=line[72:76]
                element=line[76:78].strip()
                if current_segid!=segid:
                    current_segid=segid
                    structure_builder.init_seg(current_segid)
                if current_chain_id!=chainid:
                    current_chain_id=chainid
                    structure_builder.init_chain(current_chain_id)
                    current_residue_id=residue_id
                    current_resname=resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException, message:
                        self._handle_PDB_exception(message, global_line_counter)
                elif current_residue_id!=residue_id or current_resname!=resname:
                    current_residue_id=residue_id
                    current_resname=resname
                    try:
                        structure_builder.init_residue(resname, hetero_flag, resseq, icode)
                    except PDBConstructionException, message:
                        self._handle_PDB_exception(message, global_line_counter) 
                # init atom
                try:
                    structure_builder.init_atom(name, coord, bfactor, occupancy, altloc,
                                                fullname, serial_number, element)
                except PDBConstructionException, message:
                    self._handle_PDB_exception(message, global_line_counter)
            elif(record_type=='ANISOU'):
                anisou=map(float, (line[28:35], line[35:42], line[43:49], line[49:56], line[56:63], line[63:70]))
                # U's are scaled by 10^4 
                anisou_array=(numpy.array(anisou, 'f')/10000.0).astype('f')
                structure_builder.set_anisou(anisou_array)
            elif(record_type=='MODEL '):
                try:
                    serial_num=int(line[10:14])
                except:
                    self._handle_PDB_exception("Invalid or missing model serial number",
                                               global_line_counter)
                    serial_num=0
                structure_builder.init_model(current_model_id,serial_num)
                current_model_id+=1
                model_open=1
                current_chain_id=None
                current_residue_id=None
            elif(record_type=='END   ' or record_type=='CONECT'):
                # End of atomic data, return the trailer
                self.line_counter=self.line_counter+local_line_counter
                return coords_trailer[local_line_counter:]
            elif(record_type=='ENDMDL'):
                model_open=0
                current_chain_id=None
                current_residue_id=None
            elif(record_type=='SIGUIJ'):
                # standard deviation of anisotropic B factor
                siguij=map(float, (line[28:35], line[35:42], line[42:49], line[49:56], line[56:63], line[63:70]))
                # U sigma's are scaled by 10^4
                siguij_array=(numpy.array(siguij, 'f')/10000.0).astype('f')   
                structure_builder.set_siguij(siguij_array)
            elif(record_type=='SIGATM'):
                # standard deviation of atomic positions
                sigatm=map(float, (line[30:38], line[38:45], line[46:54], line[54:60], line[60:66]))
                sigatm_array=numpy.array(sigatm, 'f')
                structure_builder.set_sigatm(sigatm_array)
            local_line_counter=local_line_counter+1
        # EOF (does not end in END or CONECT)
        self.line_counter=self.line_counter+local_line_counter
        return []

    def _handle_PDB_exception(self, message, line_counter):
        """
        This method catches an exception that occurs in the StructureBuilder
        object (if PERMISSIVE), or raises it again, this time adding the 
        PDB line number to the error message.
        """
        message="%s at line %i." % (message, line_counter)
        if self.PERMISSIVE:
            # just print a warning - some residues/atoms may be missing
            warnings.warn("PDBConstructionException: %s\n"
                          "Exception ignored.\n"
                          "Some atoms or residues may be missing in the data structure."
                          % message, PDBConstructionWarning)
        else:
            # exceptions are fatal - raise again with new message (including line nr)
            raise PDBConstructionException(message)


if __name__=="__main__":

    import sys

    p=PDBParser(PERMISSIVE=True)

    filename = sys.argv[1]
    s=p.get_structure("scr", filename)

    for m in s:
        p=m.get_parent()
        assert(p is s)
        for c in m:
            p=c.get_parent()
            assert(p is m)
            for r in c:
                print r
                p=r.get_parent()
                assert(p is c)
                for a in r:
                    p=a.get_parent()
                    if not p is r:
                        print p, r
                    
                
                
        
