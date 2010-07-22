# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# NACCESS interface adapted from Bio/PDB/DSSP.py

import os, sys, tempfile
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.AbstractPropertyMap import AbstractResiduePropertyMap, AbstractAtomPropertyMap

"""Interface for the program NACCESS.

See: http://wolf.bms.umist.ac.uk/naccess/

errors likely to occur with the binary:
default values are often due to low default settings in accall.pars
- e.g. max cubes error: change in accall.pars and recompile binary

use naccess -y, naccess -h or naccess -w to include HETATM records
"""

def run_naccess(model, pdb_file, probe_size = None, z_slice = None, \
                naccess = 'naccess', temp_path = '/tmp/'):
    
    # make temp directory; chdir to temp directory, 
    # as NACCESS writes to current working directory
    tmp_path = tempfile.mktemp(dir = temp_path)
    os.mkdir(tmp_path)
    old_dir = os.getcwd()
    os.chdir(tmp_path)
    
    # file name must end with '.pdb' to work with NACCESS
    # -> create temp file of existing pdb
    #    or write model to temp file
    tmp_pdb_file = tempfile.mktemp('.pdb', dir = tmp_path)
    if pdb_file:
        os.system('cp %s %s' % (pdb_file, tmp_pdb_file))
    else:
        writer = PDBIO()
        writer.set_structure(model.get_parent())
        writer.save(tmp_pdb_file)

    # create the command line and run
    # catch standard out & err
    command = '%s %s ' % (naccess, tmp_pdb_file)
    if probe_size:
        command += '-p %s ' % probe_size
    if z_slice:
        command += '-z %s ' % z_slice
    in_, out, err = os.popen3(command)
    in_.close()
    stdout = out.readlines()
    out.close()
    stderr = err.readlines()
    err.close()

    # get the output, then delete the temp directory
    rsa_file = tmp_pdb_file[:-4] + '.rsa'
    rf = open(rsa_file)
    rsa_data = rf.readlines()
    rf.close()
    asa_file = tmp_pdb_file[:-4] + '.asa'
    af = open(asa_file)
    asa_data = af.readlines()
    af.close()
    os.chdir(old_dir)
    os.system('rm -rf %s >& /dev/null' % tmp_path)
    return rsa_data, asa_data

def process_rsa_data(rsa_data):
    # process the .rsa output file: residue level SASA data
    naccess_rel_dict = {}
    for line in rsa_data:
        if line.startswith('RES'):
            res_name = line[4:7]
            chain_id = line[8]
            resseq = int(line[9:13])
            icode = line[13]
            res_id = (' ', resseq, icode)
            naccess_rel_dict[(chain_id, res_id)] = { \
                'res_name': res_name,
                'all_atoms_abs': float(line[16:22]),
                'all_atoms_rel': float(line[23:28]),
                'side_chain_abs': float(line[29:35]),
                'side_chain_rel': float(line[36:41]),
                'main_chain_abs': float(line[42:48]),
                'main_chain_rel': float(line[49:54]),
                'non_polar_abs': float(line[55:61]),
                'non_polar_rel': float(line[62:67]),
                'all_polar_abs': float(line[68:74]),
                'all_polar_rel': float(line[75:80]) } 
    return naccess_rel_dict

def process_asa_data(rsa_data):
    # process the .asa output file: atomic level SASA data
    naccess_atom_dict = {}
    for line in rsa_data:
        atom_serial = line[6:11]
        full_atom_id = line[12:16]
        atom_id = full_atom_id.strip()
        altloc = line[16]
        resname = line[17:20]
        chainid = line[21]
        resseq = int(line[22:26])
        icode = line[26]
        res_id = (' ', resseq, icode)
        id = (chainid, res_id, atom_id)
        asa = line[54:62]               # solvent accessibility in Angstrom^2
        vdw = line[62:68]               # van der waal radius
        naccess_atom_dict[id] = asa
    return naccess_atom_dict


class NACCESS(AbstractResiduePropertyMap):
    
    def __init__(self, model, pdb_file = None,
                 naccess_binary = 'naccess', tmp_directory = '/tmp'):
        res_data, atm_data = run_naccess(model, pdb_file, naccess = naccess_binary,
                                         temp_path = tmp_directory)
        naccess_dict = process_rsa_data(res_data)
        res_list = []
        property_dict={}
        property_keys=[]
        property_list=[]
        # Now create a dictionary that maps Residue objects to accessibility
        for chain in model:
            chain_id=chain.get_id()
            for res in chain:
                res_id=res.get_id()
                if (chain_id, res_id) in naccess_dict:
                    item = naccess_dict[(chain_id, res_id)]
                    res_name = item['res_name']
                    assert (res_name == res.get_resname())
                    property_dict[(chain_id, res_id)] = item
                    property_keys.append((chain_id, res_id))
                    property_list.append((res, item))
                    res.xtra["EXP_NACCESS"]=item
                else:
                    pass
        AbstractResiduePropertyMap.__init__(self, property_dict, property_keys, 
                property_list)

class NACCESS_atomic(AbstractAtomPropertyMap):

    def __init__(self, model, pdb_file = None,
                 naccess_binary = 'naccess', tmp_directory = '/tmp'):
        res_data, atm_data = run_naccess(model, pdb_file, naccess = naccess_binary,
                                         temp_path = tmp_directory)
        self.naccess_atom_dict = process_asa_data(atm_data)
        atom_list = []
        property_dict={}
        property_keys=[]
        property_list=[]
        # Now create a dictionary that maps Atom objects to accessibility
        for chain in model:
            chain_id = chain.get_id()
            for residue in chain:
                res_id = residue.get_id()
                for atom in residue:
                    atom_id = atom.get_id()
                    full_id=(chain_id, res_id, atom_id)
                    if full_id in self.naccess_atom_dict:
                        asa = self.naccess_atom_dict[full_id]
                        property_dict[full_id]=asa
                        property_keys.append((full_id))
                        property_list.append((atom, asa))
                        atom.xtra['EXP_NACCESS']=asa
        AbstractAtomPropertyMap.__init__(self, property_dict, property_keys, 
                property_list)


if __name__=="__main__":
    
    import sys
    from Bio.PDB import PDBParser
    
    p=PDBParser()
    s=p.get_structure('X', sys.argv[1])
    model=s[0]

    n = NACCESS(model, sys.argv[1])
    for e in n.get_iterator():
        print e

