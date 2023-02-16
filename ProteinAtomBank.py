import os
import mdtraj as md

class AtomBankForProtein:
    def __init__(self):
        self._recPDBFile = "rec_amber.pdb"

    @classmethod
    def CA(self):
        return ["CA"]

    @classmethod
    def backbone(self):
        return ["CA","C","N"]

    @classmethod
    def heavyAtoms(self,residByName):
        atoms = {
                 "ALA":['N','C','O','CA','CB'],
                 "VAL":['N','C','O','CA','CB','CG1','CG2'],
                 "LEU":['N','C','O','CA','CB','CG','CD1','CD2'],
                 "ILE":['N','C','O','CA','CB','CG1','CG2','CD1'],
                 "PHE":['N','C','O','CA','CB','CG','CD1','CE1','CZ','CE2','CD2'],
                 "TRP":['N','C','O','CA','CB','CG','CD1','NE1','CE2','CZ2','CH2','CZ3','CE3','CD2'],
                 "MET":['N','C','O','CA','CB','CG','SD','CE'],
                 "PRO":['N','C','O','CD','CG','CB','CA'],
                 "ASP":['N','C','O','CA','CB','CG','OD1','OD2'],
                 "GLU":['N','C','O','CA','CB','CG','CD','OE1','OE2'],
                 "LYS":['N','C','O','CA','CB','CG','CD','CE','NZ'],
                 "ARG":['N','C','O','CA','CB','CG','CD','NE','CZ','NH1','NH2'],
                 "HIS":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
                 "HID":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
                 "HIE":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
                 "HIP":['N','C','O','CA','CB','CG','ND1','CE1','NE2','CD2'],
                 "GLY":['N','C','O','CA'],
                 "SER":['N','C','O','CA','CB','OG'],
                 "THR":['N','C','O','CA','CB','CG2','OG1'],
                 "CYS":['N','C','O','CA','CB','SG'],
                 "CYX":['N','C','O','CA','CB','SG'],
                 "TYR":['N','C','O','CA','CB','CG','CD1','CE1','CZ','OH','CE2','CD2'],
                 "ASN":['N','C','O','CA','CB','CG','OD1','ND2'],
                 "GLN":['N','C','O','CA','CB','CG','CD','OE1','NE2'],
        }
        assert residByName in atoms, "Please specify a valid protein residue name."
        return atoms[residByName]
