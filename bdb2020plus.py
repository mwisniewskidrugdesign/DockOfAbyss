from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

import settings


pd.set_option('display.max_columns', None)

generate_library_step = True
convert_to_pdbqt_step=True
smina_docking_step=True

def pipeline(datadir,bdb2020plus_datadir,bdb2020plus_df):

    if generate_library_step:
        '''Generate the library for bdb2020plus operations'''
        generator.generate_libraray(datadir)
        workspace = generator.GetDataset(datadir,bdb2020plus_df)                                        #create workspace
        workspace.bdb2020plus(bdb2020plus_datadir)                                          #copy files to workspace

    if convert_to_pdbqt_step:
        for index,row in bdb2020plus_df.iterrows():                                                     #iteruje po dataframe
            print(str(index)+'.'+row['pdbid'])
            library = generator.Converter(row['pdbid'],datadir)                             #Klasa konwertera
            library.to_pdbqt(protein=True,ligand=True,native_ligand=True)                    #konwertuje do pdbqt

    if smina_docking_step:
        bdb2020plus_df_prep = datasets.DatasetPreparation(bdb2020plus_df)
        molecules = bdb2020plus_df_prep.get_molecules('pdbid')           #Generate list of molecules
        print(molecules)

        bdb2020plus_docking = docking.Smina(datadir)
        bdb2020plus_docking.smina_dirs()
        smina_matrix = bdb2020plus_docking.create_smina_matrix(molecules[:],molecules[:],50)

        for molecule_idx,molecule in enumerate(molecules[:]):
            print('Docking '+molecule+' to '+molecule+'.')
            bdb2020plus_docking.smina_files(molecule,molecule,molecule)
            bdb2020plus_docking.smina_docking(50)
            bdb2020plus_docking.read_experimental_affinity(bdb2020plus_df,molecule,molecule)
            bdb2020plus_docking.read_scoring_function()
            bdb2020plus_docking.read_atom_term_function(50)

            bdb2020plus_docking.fill_smina_matrix(molecule_idx,molecule_idx)

        bdb2020plus_docking.save_matrix('smina_matrix')


