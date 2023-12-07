from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

pd.set_option('display.max_columns', None)

generate_library_step = True
convert_step=True                #ADD IF !!!!!!!!!!!!!!!!!!!!!!!!!! SUCH US DOCKING PROGRAMS LIST
docking_step=True
docking_programs=['smina']

def pipeline(datadir,bdb2020plus_datadir,bdb2020plus_df,no_modes):

    if generate_library_step:
        '''Generate the library for bdb2020plus operations'''
        generator.generate_libraray(datadir)
        workspace = generator.GetDataset(datadir,bdb2020plus_datadir,bdb2020plus_df)                                    #create workspace
        workspace.bdb2020plus()                                                  #copy files to workspace

    if convert_step:
        for index,row in bdb2020plus_df.iterrows():                                                 #   Converter loop for every PDB structure from BDB2020+ set
            print(str(index)+'.'+row['pdbid'])                                                      ##  Print PDB structure code
            library = generator.Converter(row['pdbid'],datadir)                                     ##  Converter Class: library
            library.pdb_to_pdbqt(protein=True,ligand=True,native_ligand=True)                       ##  Convert files from pdb to pdbqt (define which ones: protein, ligand, native_ligand)
            library.pdb_to_mol(protein=True,ligand=True,native_ligand=True)                         ##  Convert files from pdb to mol (define which ones)

    if docking_step:
        bdb2020plus_df_prep = datasets.DatasetPreparation(bdb2020plus_df)                           #   Generate Molecule list Class - is it neccessery in this case?
        molecules = bdb2020plus_df_prep.get_molecules('pdbid')                                      ##  Generating molecules list from PDB structure codes

        if 'smina' in docking_programs:
            bdb2020plus_smina_docking = docking.Smina(datadir)                                      # SMINA Docking Class
            bdb2020plus_smina_docking.smina_dirs()                                                              ## Generate output dirs for SMINA docking
            smina_matrix = bdb2020plus_smina_docking.create_smina_matrix(molecules[:2],molecules[:2],no_modes)  ## Generate empty SMINA outputs matrix

        # if 'rxdock' in docking_programs:
        #     bdb2020plus_rx_docking = docking.RxDock(datadir)
        #     bdb2020plus_rx_docking.rxdock_dirs()
        #     rxdock_matrix = None

        for molecule_idx,molecule in enumerate(molecules[:2]):                                      ## Docking Loop for molecules from list generated earlier
            print('Docking '+molecule+' to '+molecule+'. With: \n',docking_programs)                ## Print PDB structure code
            if 'smina' in docking_programs:
                smina_docking_error_number = 0
                while True:                                                                         ## Loop - necessery to generate exactly 100 modes, sometimes with random seed it's generating smaller number of conforms
                    try:
                        bdb2020plus_smina_docking.smina_files(molecule,molecule,molecule)           ## set variables for specific molecule files
                        bdb2020plus_smina_docking.smina_docking(no_modes)                           ## smina docking function
                        bdb2020plus_smina_docking.read_experimental_affinity(bdb2020plus_df,molecule,molecule)  ## reading experimental affinity data for specific molecule from dataframe
                        bdb2020plus_smina_docking.read_scoring_function()                                       ## reading scoring function predicted binding affinity from output
                        bdb2020plus_smina_docking.read_atom_term_function(no_modes)                             ## reading atom terms sf's components from output
                        bdb2020plus_smina_docking.fill_smina_matrix(molecule_idx,molecule_idx)                  ## fill smina matrix with output datas                      ### MAYBE inserts read_* functions into this one
                    except:
                        smina_docking_error_number+=1
                        print('Smina proposed less modes than expected for docking '+molecule+' to '+molecule+'. '+str(smina_docking_error_number)+'st time.')
                        continue
                    break
            # if 'rxdock' in docking_programs:
            #     rx_docking_error_number = 0
            #     bdb2020plus_rx_docking.rxdock_files(molecule,molecule,molecule)
            #     bdb2020plus_rx_docking.rxdock_system_preparation()
            #     #while True:
            #         #try:
            #             #print(elo)
            #         #except:
            #             #rx_docking_error_number+=1
            #             #print('RxDock proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + rx_docking_error_number + 'st time.')

        if 'smina' in docking_programs:
            bdb2020plus_smina_docking.save_matrix('smina_matrix')                                       ## save SMINA matrix


