# Dock of Abyss

Dataset library:
```commandline
├── /proteins
│ ├── /pdb                  #pdb protein files
│ ├── /pdbqt                #pdbqt protein files for smina
│ └── /fasta                #sequence protein files for clustering
├── /ligands
│ ├── /sdf                  #sdf ligand files
│ ├── /mol2                 #mol2 ligand files
│ ├── /pdb                  #pdb ligand files
│ ├── /pdbqt                #pdbqt ligand files for smina
│ └── /smi                  #from generating minimized ligands
├── /native_ligands     #from PDB complex binding pocket
│ ├── /sdf                  #sdf native_ligand files
│ ├── /mol2                 #mol2 native_ligand files
│ ├── /pdb                  #pdb native_ligand files
│ ├── /pdbqt                #pdbqt native_ligand files for smina
│ └── /smi                  #for generating minimized ligands
├── /logs               #log files for pipeline operations
├── /docs               #documents
├── /docking_scores
│ ├── /smina
│ │ ├── /pdbqt          #smina files with modes conformations
│ │ ├── /log            #smina files with modes affinity and rmsd
│ │ └── /atom_terms     #smina files with sf components per atom
```
## Molecular Docking Software
### Autodock Smina
Autodock Smina was used with standard scoring function weights.

The receptors were prepared with **MGL-tools Python package script**: *prepare_receptor4.py* and ligands with **openbabel**.

**The command line:**
```commandline
smina -r <receptor.pdbqt> -l <ligand.pdbqt> --autobox_ligand <native_ligand.pdbqt> --autobox_add 8 --exhaustiveness 64 --num_modes <no_modes> -o <output_file.pdbqt> --atom_terms <atom_terms_output_file.txt> --log output_file.log --atom_term_data --cpu 12 --min_rmsd_filter 0 --energy_range 1000000
```

### RxDock
RxDock command:
```commandline
reference ligand method for mapping cavity
scoring function:
S_total = S_inter + S_intra + S_site + S_restrain
S_inter = W_inter_vdW * S_inter_vdW + W_inter_polar * S_inter_polar ...
standard scoring function file = RbtInterIdxSF.prm 
https://rxdock.gitlab.io/documentation/devel/html/reference-guide/scoring-functions.html

check diffrent scoring functions:

RbtInterGridSF.prm      #As above, but vdW term uses a precalculated grid
RbtSolvIdxSF.prm        #Intermolecular scoring function definition (desolvation scoring function, SF5)
calcgrid_vdw1.prm       #vdW term only (ECUT = 1), for calculating vdW grid (used by rbcalcgrid)
calcgrid_vdw5.prm       #vdW term only (ECUT = 5), for calculating vdW grid (used by rbcalcgrid)
Tripos52_vdw.prm        #vdW term parameter file
Tripos52_dihedrals.prm  #Dihedral term parameter file
solvation_asp.prm       #Desolvation term parameter file


SystemPreparation

```
