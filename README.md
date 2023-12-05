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
Smina command:
```commandline
smina -r <receptor.pdbqt> -l <ligand.pdbqt> --autobox_ligand <native_ligand.pdbqt> --autobox_add 8 --exhaustiveness 64 --num_modes <no_modes> -o <output_file.pdbqt> --atom_terms <atom_terms_output_file.txt> --log output_file.log --atom_term_data --cpu 12 --min_rmsd_filter 0 --energy_range 1000000

```

RxDock command:
