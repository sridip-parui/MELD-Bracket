import parmed as pmd ## best friends!!
from parmed.modeller.residue import ResidueTemplate
from parmed.amber.offlib import AmberOFFLibrary
from numpy import mean
import numpy as np
from numpy import linalg as LA


class MoleculeBase(object):
    '''
    Base class for curating related lib, frcmod files as well as preparing tleap
    input string and performing translational moves.
    '''
    
    def __init__(self):
        self._translation_vector = np.zeros(3)
        self._disulfide_list = []
        self._frcmod_files = []
        self._lib_files = []
        
    def set_translation_vector(self,translation_vector):
        self._translation_vector = np.array(translation_vector)
        
    def add_disulfide(self,res_index_i,res_index_j):
        self._disulfide_list.append((res_index_i, res_index_j))
        
    def add_frcmod_file(self, fname):
        self._frcmod_files.append(fname)
    
    def add_lib_file(self, fname):
        self._lib_files.append(fname)

    ################## helper functions ########################
    def _gen_translation_string(self,mol_id):
        x=self._translation_vector[0]
        y=self._translation_vector[1]
        z=self._translation_vector[2]
        return f"""translate {mol_id} {{ {x} {y} {z} }}"""
    
    def _gen_disulfide_string(self, mol_id):
        disulfide_strings = []
        for i, j in self._disulfide_list:
            d = f"bond {mol_id}.{i}.SG {mol_id}.{j}.SG"
            disulfide_strings.append(d)
        return disulfide_strings
    
    def _gen_read_frcmod_string(self):
        frcmod_string = []
        for p in self._frcmod_files:
            frcmod_string.append(f"loadAmberParams {p}")
        return frcmod_string
    
    def _gen_read_lib_string(self):
        lib_string = []
        for p in self._lib_files:
            lib_string.append(f"loadoff {p}")
        return lib_string

    
class MoleculeFromMol2(MoleculeBase):
    
    def __init__(self, mol2_path,mol_id=None,unit_name=None):
        super(MoleculeFromMol2, self).__init__()
        self._mol2_path = mol2_path
        self._mol_id = mol_id
        self._unit_name = unit_name
        
    def prepare_for_tleap(self, mol_id,unit_name):
        '''
        This function update the unit_name in old lib and mol2 file to the new one\
        (unit_name). It relies on parmed to do the heavy lifting.
        '''
        if len(self._lib_files) != 1:
            raise Exception("Did you properly add lib and frcmod file?")
            
        updated_mol2_path = self._update_(self._mol2_path,unit_name=unit_name)
        self._mol2_path = updated_mol2_path 

        assert len(self._lib_files) == 1
        updated_offlib_path = self._update_(self._lib_files[0],unit_name=unit_name)
        self._lib_files = [updated_offlib_path]
        
    def _update_(self,libORmol2File,unit_name=None):
        '''Update the name of the residue template in libORmol2File'''
        libORmol2 = pmd.load_file(libORmol2File)
        if isinstance(libORmol2,dict): # It is a OFFLibrary file
            if len(libORmol2) == 1:
                _,residueTemplate = libORmol2.popitem()
            else:
                raise Exception(f"Warning: More than one template in {libORmol2File}.")
        elif isinstance(libORmol2,ResidueTemplate): # It is a mol2 file
            residueTemplate = libORmol2
        else:
            raise Exception(f"I don't know how to handle {libORmol2File}.")
        
        residueTemplate.name = unit_name    # update the name of the template
        mol_id, fl_type = libORmol2File.split(".")
        newFileName = f"{mol_id}_{unit_name}.{fl_type}"
        residueTemplate.save(newFileName)
        return newFileName
    
    def prepare_translation_vector(
        self,ref_mol_path=None,
        anchor_xyz=None,offset=None
        ):
        '''
        This function takes the ref_mol and compute the translational vector 
        that will move the curr_mol(self) from it by offset (in $\AA$).
        '''
        mol_xyz = pmd.load_file(self._mol2_path).coordinates
        mol_com_xyz = mean(mol_xyz,axis=0)

        if (ref_mol_path is not None and 
            os.path.isfile(ref_mol_path) and 
            offset is not None):
            # check whether a valide file or a offset vector is provided.
            ref_mol_xyz = pmd.load_file(ref_mol_path).coordinates
            ref_com_xyz = mean(ref_mol_xyz,axis=0)
            com_v = (mol_com_xyz - ref_com_xyz)
            com_v /= LA.norm(com_v)
            self._translation_vector = com_v*offset
        
        elif type(anchor_xyz).__module__ == np.__name__:
            # check if a valide numpy array vector is given.
            self._translation_vector = anchor_xyz - mol_com_xyz
        else:
            raise Exception("I don't know how to handle the input.")
            
    def generate_tleap_input(self, mol_id,unit_name):
        leap_cmds = []
        ##leap_cmds.append("source leaprc.gaff")
        leap_cmds.extend(self._gen_read_frcmod_string())
        leap_cmds.extend(self._gen_read_lib_string())
        leap_cmds.append(f"mol_{mol_id} = loadmol2 {mol_id}_{unit_name}.mol2")
            #.format(mol_id=mol_id,unit_name=unit_name))
        leap_cmds.extend(self._gen_disulfide_string("mol_{}".format(mol_id)))
        leap_cmds.append(self._gen_translation_string("mol_{}".format(mol_id)))
        return leap_cmds 

    
class MoleculeFromPdb(MoleculeBase):
    
        def __init__(self, pdb_path,mol_id=None):
            super(MoleculeFromPdb, self).__init__()
            self._mol_id = mol_id
            self._unit_name = None
            self._pdb_path = pdb_path
            
        def prepare_for_tleap(self,mol_id,unit_name):
            '''
            Just a place holder.
            '''
            pass 
            
        def generate_tleap_input(self, mol_id,unit_name=None):
            leap_cmds = []
            ##leap_cmds.append("source leaprc.gaff")
            leap_cmds.extend(self._gen_read_frcmod_string())
            leap_cmds.extend(self._gen_read_lib_string())
            #leap_cmds.append("mol_{mol_id} = loadpdb {mol_id}.pdb".format(mol_id=mol_id))
            leap_cmds.append("mol_{mol_id} = loadpdb {rec_fl}".format(mol_id=mol_id, rec_fl=self._pdb_path))
            leap_cmds.extend(self._gen_disulfide_string("mol_{}".format(mol_id)))
            leap_cmds.append(self._gen_translation_string("mol_{}".format(mol_id)))
            return leap_cmds


import subprocess


class systemBuilder(object):
    
    def __init__(
        self, forcefield=["ff14sb","gaff2"],explicit_solvent=True, 
        solvent_forcefield="tip3p", solvent_distance=8.,
        explicit_ions=True, p_ion="Na+", n_ion="Cl-"
        ):
            
        self._forcefield = []
        self._set_forcefield(forcefield)

        self._explicit_solvent = explicit_solvent
        self._explicit_ions = explicit_ions
        if self._explicit_solvent:
            self._set_solvent_forcefield(solvent_forcefield)
            self._solvent_dist = solvent_distance

            if self._explicit_ions:
                self._set_positive_ion_type(p_ion)
                self._set_negative_ion_type(n_ion)
                
    def build_system_from_molecules(
        self, molecules, patchers=None,
        leap_header_cmds=None, output_prefix=None
        ):  
        if patchers is None:
            patchers = []
        if leap_header_cmds is None:
            leap_header_cmds = []
        if isinstance(leap_header_cmds, str):
            leap_header_cmds = [leap_header_cmds]
        if output_prefix is None:
            raise Exception("Please provide a proper output prefix.")

        leap_cmds = []
        mol_ids = []
        leap_cmds.extend(self._generate_leap_header())
        leap_cmds.extend(leap_header_cmds)
        for mol in molecules:
            mol_id = mol._mol_id
            mol_ids.append(mol_id)
            mol.prepare_for_tleap(mol_id,mol._unit_name)
            leap_cmds.extend(mol.generate_tleap_input(mol_id,mol._unit_name))
        if self._explicit_solvent:
            leap_cmds.extend(self._generate_solvent(mol_ids))
            leap_cmds.extend(self._generate_leap_footer(["solute"],output_prefix))
        else:
            leap_cmds.extend(self._generate_leap_footer(mol_ids))
        
        with open("leap.parm.in","w") as leap_prefix:
            for cmd in leap_cmds:
                if cmd.startswith("solvateBox"):
                    break
                leap_prefix.write(cmd+"\n")
        
        with open("tleap.in", "w") as tleap_file:
            tleap_string = "\n".join(leap_cmds)
            tleap_file.write(tleap_string)
        subprocess.check_call("tleap -f tleap.in > tleap.out", shell=True)
        
    
    def _set_forcefield(self, forcefield):
        ff_dict = {
            "ff12sb": "leaprc.ff12SB",
            "ff14sb": "leaprc.protein.ff14SB",
            "ff14sbside": "leaprc.protein.ff14SBonlysc",
            "gaff2": "leaprc.gaff2",
            "gaff2_mod": "leaprc.gaff2_FClBrS"
            }
        for ff in forcefield:
            try:
                self._forcefield.append(ff_dict[ff])
            except KeyError:
                raise RuntimeError(f"Unknown forcefield: {ff}")

    #def _set_gb_radii(self, gb_radii):
    #    allowed = ["mbondi2", "mbondi3"]
    #    if gb_radii not in allowed:
    #        raise RuntimeError("Unknown gb_radii: {gb_radii}".format(gb_radii=gb_radii))
    #    else:
    #        self._gb_radii = gb_radii

    def _set_solvent_forcefield(self, solvent_forcefield):
        ff_dict = {
            "spce": "leaprc.water.spce",
            "spceb": "leaprc.water.spceb",
            "opc":  "leaprc.water.opc",
            "tip3p": "leaprc.water.tip3p",
            "tip4pew": "leaprc.water.tip4pew",
            }
        
        box_dict = {
            "spce": "SPCBOX",
            "spceb": "SPCBOX",
            "opc":  "OPCBOX",
            "tip3p": "TIP3PBOX",
            "tip4pew": "TIP4PEWBOX",
            }
        try:
            self._solvent_forcefield = ff_dict[solvent_forcefield]
            self._solvent_box = box_dict[solvent_forcefield]
        except KeyError:
            raise RuntimeError(f"Unknown solvent_model: {solvent_model}")
            
    def _set_positive_ion_type(self, ion_type):
        allowed = ["Na+", "K+", "Li+", "Rb+", "Cs+", "Mg+"]
        if ion_type not in allowed:
            raise RuntimeError(f"Unknown ion_type: {ion_type}")
        else:
            self._p_ion = ion_type

    def _set_negative_ion_type(self, ion_type):
        allowed = ["Cl-", "I-", "Br-", "F-"]
        if ion_type not in allowed:
            raise RuntimeError(f"Unknown ion_type: {ion_type}")
        else:
            self._n_ion = ion_type

    def _generate_leap_header(self):
        leap_cmds = []

        for ff in self._forcefield:
            leap_cmds.append(f"source {ff}")
        if self._explicit_solvent:
            leap_cmds.append(f"source {self._solvent_forcefield}")
        return leap_cmds
    
    def _generate_solvent(self, mol_ids):
        leap_cmds = []
        list_of_mol_ids = ""
        for mol_id in mol_ids:
            list_of_mol_ids += f"mol_{mol_id} "
        
        leap_cmds.append(f"solute = combine {{ {list_of_mol_ids} }}")
         
        if self._explicit_ions:
            leap_cmds.append(f"addIons solute {self._p_ion} 0")
            leap_cmds.append(f"addIons solute {self._n_ion} 0")
        
        leap_cmds.append("check solute")    # check integrity of the solute 
        leap_cmds.append(f"solvateBox solute {self._solvent_box} {self._solvent_dist}")
        
        return leap_cmds

    def _generate_leap_footer(self, mol_ids,output_prefix):
        leap_cmds = []
        list_of_mol_ids = ""
        for mol_id in mol_ids:
            list_of_mol_ids += f"{mol_id} "

        leap_cmds.append(f"sys = combine {{ {list_of_mol_ids} }}")
        leap_cmds.append(f"saveAmberParm sys {output_prefix}.prmtop {output_prefix}.inpcrd")

        leap_cmds.append("quit")
        return leap_cmds

import os
import shutil


class FileBase(object):
    
    def __init__(self,fileSrc=None,destDir=None):
        self._fileSrc_ = fileSrc
        self._destDir_ = destDir
        self._checkFile_()    # make sure we get a valid file.
    
    def copyFile(self):
        shutil.copy(self._fileSrc_,self._destDir_)
    
    def _checkFile_(self):
        if (self._fileSrc_ is None or
           not os.path.isfile(self._fileSrc_)
           ):
            raise Exception("Please provide a valid file.")
        
        if (self._destDir_ is None):
            raise Exception("Please provid a destination.")
        else: # we get a destination
            if not os.path.isdir(self._destDir_): # dest dir not exist
                try:
                    os.makedirs(self._destDir_)   # try to create one
                except OSError as e:          # failed to create it.
                    print(f"Failed to create due to {e.strerror} on {e.filename}")
    
import parmed as pmd ## best friends!!
import numpy as np
import mdtraj as md
import math


class SampleSurfacePoint:
    '''
    This function sample points on a suface defined by radius and a macromolecule 
    PDB file. It then returns a list of points
    '''
        
    def __init__(self,protPDB=None,radius=None,surfPointSep=None,NPoints=None):
        assert protPDB is not None and radius is not None and surfPointSep is not None \
        and NPoints is not None
        
        self._center = pmd.load_file(protPDB).coordinates.mean(axis=0)
        self._r = radius
        self._NPoints = NPoints
        self._separation = surfPointSep
        
    def _getPointFastest(self):
        """
        Idea adopts from here: https://karthikkaranth.me/blog/generating-random-points-in-a-sphere/
        """

        u,v=np.random.random(2)
        theta = u * 2.0 * math.pi
        phi = math.acos(2.0 * v - 1.0)
        sinTheta = math.sin(theta)
        cosTheta = math.cos(theta)
        sinPhi = math.sin(phi)
        cosPhi = math.cos(phi)
        x = self._r * sinPhi * cosTheta
        y = self._r * sinPhi * sinTheta
        z = self._r * cosPhi
        xyz = np.array([x,y,z]) + self._center
        return xyz

    def getPointAwayFromEachOther(self):
    
        candidates = []
        if len(candidates) == 0:
            candidates.append(self._getPointFastest())
        while len(candidates) < self._NPoints:
            tmpPoint = self._getPointFastest()
            xyz = np.array(candidates)
            if (np.sqrt(np.square(xyz - tmpPoint).sum(axis=1)) > self._separation).all():
                candidates.append(tmpPoint)
        return np.array(candidates)
    
    def getPointDict(self,mol_name_list):
        xyz_dict = {}
        xyz_mat = self.getPointAwayFromEachOther()
        assert xyz_mat.shape[0] == len(mol_name_list)-1,"Too many mols."
        for i in range(len(mol_name_list)-1):
            xyz_dict[i] = xyz_mat[i,:]
        #for i in range(len(mol_name_list)):
        #    xyz_dict[mol_name_list[i]] = xyz_mat[i,:]
        return xyz_dict
        
def count_wat(prmtop):
    top = pmd.load_file(prmtop)
    watIdx = 0
    nonWatIonIdx = 0
    # this for loop is very expensive, I need to optimize it.
    for resid in top.residues:
        if resid.name != "Cl-" and resid.name != "Na+":
            nonWatIonIdx += 1
    
        if resid.name == "WAT":
            watIdx += 1
    nonWatIonIdx -= watIdx
    return watIdx, nonWatIonIdx



class solventEquilizer:
    
    def __init__(
        self,targetNumOfSolvent=None,buffer=None,
        top=None,crd=None,
        tol=None,soluteres=None,
            leapin="leap.parm.in",molname="solute",loadpdb="no",boxMode=1
    ):
        self._targetNumOfSolvent = targetNumOfSolvent
        self._buffer = buffer
        self._top = top
        self._crd = crd
        self._tol = tol
        self._numOfSolute = soluteres
        self._leapin = leapin
        self._molname = molname
        self._loadpdb = loadpdb
        self._boxMode = boxMode # default rectangular box
        self._excutable = "/home/frank/Script/Snippet/Solvate.sh"
        
    def run(self):
        with open("equilize.in","w") as fl:
            fl.write(f"target {self._targetNumOfSolvent}\n")
            fl.write(f"buffer {self._buffer}\n")
            fl.write(f"top    {self._top}\n")
            fl.write(f"crd    {self._crd}\n")
            fl.write(f"loadpdb {self._loadpdb}\n")
            fl.write(f"leapin {self._leapin}\n")
            fl.write(f"tol    {self._tol}\n")
            fl.write(f"mode {self._boxMode}\n")
            fl.write(f"soluteres {self._numOfSolute}\n")
            fl.write(f"molname {self._molname}\n")
            
        subprocess.check_call(f"{self._excutable} equilize.in > equilize.out", 
                              shell=True
                             )
