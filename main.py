from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import pyvista as pv

c_color = "black"
c_radius = 0.3
h_color = "blue"
h_radius = 0.2

#n = int(input("n: "))
#smiles = "C"*n
smiles = input("SMILES:")
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(mol)
conf = mol.GetConformer()

class Atom:
    def __init__(self, symbol, idx, position, color, radius):
        self.symbol = symbol
        self.idx = idx
        self.position = position
        self.color = color
        self.radius = radius
        self.bonds = []
    def add_bond(self, other_atom):
        self.bonds.append(other_atom)

atoms = []


idx_to_atom = {}

for atom in mol.GetAtoms():
    idx = atom.GetIdx()
    pos = conf.GetAtomPosition(idx)
    symbol = atom.GetSymbol()
    color = c_color if symbol == "C" else h_color
    radius = c_radius if symbol == "C" else h_radius

    atom_obj = Atom(symbol, idx, (pos.x, pos.y, pos.z), color, radius)
    atoms.append(atom_obj)
    idx_to_atom[idx] = atom_obj

# Add bonds
for bond in mol.GetBonds():
    a1 = idx_to_atom[bond.GetBeginAtomIdx()]
    a2 = idx_to_atom[bond.GetEndAtomIdx()]
    a1.add_bond(a2)
    a2.add_bond(a1)  #both sides have the bond

plotter = pv.Plotter()
for atom in atoms:
    
    sphere = pv.Sphere(radius=atom.radius, center=atom.position)
    plotter.add_mesh(sphere, color=atom.color)


    for bonded_atom in atom.bonds:
        if atom.idx < bonded_atom.idx:
            start = np.array(atom.position)
            end = np.array(bonded_atom.position)
            line = pv.Line(start, end)
            plotter.add_mesh(line, color="grey", line_width=5)


plotter.show()
