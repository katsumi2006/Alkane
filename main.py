from rdkit import Chem
from rdkit.Chem import AllChem
import pyvista as pv
import pubchempy as pcp
from math import sqrt 
import numpy as np
plotter = pv.Plotter()

CPK_map = {  #more to be added
    "bond": "gray",
    "H": "white",
    "C": "black",
    "N": "blue",
    "O": "red",
    "F": "green",
    "Cl": "green",
    "Br": "darkred",
    "I": "darkviolet",
    "P": "orange"
} 
Radii = {
    "H": 0.2,
    "C": 0.3,
    "bond": 0.08,
    "N": 0.3,
    "O": 0.3,
    "F": 0.3,
    "Cl": 0.3,
    "Br": 0.3,
    "I": 0.3,
    "P": 0.3
}

offset = 0.10

choice = int(input("IUPAC (1), SMILES(2): "))

iupac_name = ""
smiles = "" 

if choice == 1:
    iupac_name = str(input("IUPAC Name:"))
    compounds = pcp.get_compounds(iupac_name, 'name')
    try:
        smiles = compounds[0].connectivity_smiles
    except:
        print("Compound or Molecule wasn't found!")
        exit()

elif choice == 2:
    smiles = str(input("SMILES: "))
    try:
        compounds = pcp.get_compounds(smiles,'smiles')
    except:
        print("Compound or Molecule wasn't found!")
        exit()
    else:
        iupac_name = compounds[0].iupac_name

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
Chem.Kekulize(mol, clearAromaticFlags=True)

AllChem.EmbedMolecule(mol, AllChem.ETKDG()) 
AllChem.UFFOptimizeMolecule(mol) # minimize the energy level

conf = mol.GetConformer()

for atom in mol.GetAtoms():
    pos = np.array(conf.GetAtomPosition(atom.GetIdx()))
    sphere = pv.Sphere(radius=Radii[atom.GetSymbol()], center=pos)
    plotter.add_mesh(sphere, color=CPK_map[atom.GetSymbol()])
    

for bond in mol.GetBonds():
    typ = bond.GetBondType()
    Begin = bond.GetBeginAtom()
    End = bond.GetEndAtom()
    A = np.array(conf.GetAtomPosition(Begin.GetIdx()))
    B = np.array(conf.GetAtomPosition(End.GetIdx()))

    drctn = B - A
    hght = np.linalg.norm(drctn)
    unit_drctn = drctn / hght
    cntr = (A + B) * 0.5

    if abs(drctn[0]) <= abs(drctn[1]) and abs(drctn[0]) <= abs(drctn[2]):
        u = np.array([1, 0, 0])
    elif abs(drctn[1]) <= abs(drctn[0]) and abs(drctn[1]) <= abs(drctn[2]):
        u = np.array([0, 1, 0])
    else:
        u = np.array([0, 0, 1])
    
    p = np.cross(drctn, u)
    unit_p = p / np.linalg.norm(p)
    
    offsetV = unit_p * offset
    
    if typ == Chem.rdchem.BondType.SINGLE:
        cylinder = pv.Cylinder(center=cntr, direction=unit_drctn, radius=Radii["bond"], height=hght)
        plotter.add_mesh(cylinder, color=CPK_map["bond"])
    if typ == Chem.rdchem.BondType.DOUBLE:
        cylinder1 = pv.Cylinder(center=cntr + offsetV, direction=unit_drctn, radius=Radii["bond"], height=hght)
        plotter.add_mesh(cylinder1, color=CPK_map["bond"])

        cylinder2 = pv.Cylinder(center=cntr - offsetV, direction=unit_drctn, radius=Radii["bond"], height=hght)
        plotter.add_mesh(cylinder2, color=CPK_map["bond"])
    if typ == Chem.rdchem.BondType.TRIPLE:
        cylinder1 = pv.Cylinder(center=cntr, direction=unit_drctn, radius=Radii["bond"], height=hght)
        plotter.add_mesh(cylinder1, color=CPK_map["bond"])
        cylinder2 = pv.Cylinder(center=cntr + offsetV * 2, direction=unit_drctn, radius=Radii["bond"], height=hght)
        plotter.add_mesh(cylinder2, color=CPK_map["bond"])
        cylinder3 = pv.Cylinder(center=cntr - offsetV * 2, direction=unit_drctn, radius=Radii["bond"], height=hght)
        plotter.add_mesh(cylinder3, color=CPK_map["bond"])

plotter.add_text(iupac_name, font_size=10, position="upper_right")
plotter.show()
