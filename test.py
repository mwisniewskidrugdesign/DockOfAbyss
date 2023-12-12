import pandas as pd
import sys
import os

pdbqt_proteins = os.listdir('/home2/sfglab/mwisniewski/PhD/data/lp_pdbbin/protein/pdbqt')
mol_proteins = os.listdir('/home2/sfglab/mwisniewski/PhD/data/lp_pdbbin/protein/mol2')

pdb_ligands = os.listdir('/home2/sfglab/mwisniewski/PhD/data/lp_pdbbin/ligand/pdb')
pdbqt_ligands  = os.listdir('/home2/sfglab/mwisniewski/PhD/data/lp_pdbbin/ligand/pdbqt')

pdb_native_ligands = os.listdir('/home2/sfglab/mwisniewski/PhD/data/lp_pdbbin/native_ligand/pdb')
pdbqt_native_ligands  = os.listdir('/home2/sfglab/mwisniewski/PhD/data/lp_pdbbin/native_ligand/pdbqt')

df = pd.read_csv('/home2/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv')
mask = df['CL1'] == True
df = df[mask]
df = df['pdbid'].tolist()

pdbqt_proteins_to_fix = []
mol_proteins_to_fix = []
pdb_ligands_to_fix = []
pdbqt_ligands_to_fix = []
pdb_native_ligands_to_fix = []
pdbqt_native_ligands_to_fix = []

for molecule in df:

    if molecule+_'protein.pdbqt' not in pdbqt_proteins:
        pdbqt_proteins_to_fix.append(molecule)

    if molecule+'_protein.mol2' not in mol_proteins:
        mol_proteins_to_fix.append(molecule)

    if molecule+'_ligand.pdb' not in pdb_ligands:
        pdb_ligands_to_fix.append(molecule)

    if molecule + '_ligand.pdbqt' not in pdbqt_ligands:
        pdbqt_ligands_to_fix.append(molecule)

    if molecule+'_ligand.pdb' not in pdb_native_ligands:
        pdb_native_ligands_to_fix.append(molecule)

    if molecule+'_ligand.pdbqt' not in pdbqt_native_ligands:
        pdbqt_native_ligands_to_fix.append(molecule)

with open('pdbqt_proteins', 'w') as f:
    for line in pdbqt_proteins_to_fix:
        f.write(f"{line}\n")

with open('mol_proteins', 'w') as f:
    for line in mol_proteins_to_fix:
        f.write(f"{line}\n")

with open('pdb_ligands', 'w') as f:
    for line in pdb_ligands_to_fix:
        f.write(f"{line}\n")

with open('pdbqt_ligands', 'w') as f:
    for line in pdbqt_ligands_to_fix:
        f.write(f"{line}\n")

with open('pdb_native_ligands', 'w') as f:
    for line in pdb_native_ligands_to_fix:
        f.write(f"{line}\n")

with open('pdbqt_native_ligands', 'w') as f:
    for line in pdbqt_native_ligands_to_fix:
        f.write(f"{line}\n")
