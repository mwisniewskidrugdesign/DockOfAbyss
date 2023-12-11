from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

pd.set_option('display.max_columns', None)

generate_library_step = True
convert_step=True                #ADD IF !!!!!!!!!!!!!!!!!!!!!!!!!! SUCH US DOCKING PROGRAMS LIST
docking_step=True
docking_programs=['rxdock']

def diagonal_pipeline(datadir,rawdir,df,no_modes):

    if generate_library_step:
        '''Generate the library for bdb2020plus operations'''
        generator.generate_libraray(datadir)
        workspace = generator.GetDataset(datadir,rawdir,df, pdb_id_column='pdbid')                                    #create workspace
        workspace.bdb2020plus()                                                  #copy files to workspace

    if convert_step:
        for index,row in df.iterrows():                                                 #   Converter loop for every PDB structure from BDB2020+ set
            print(str(index)+'.'+row['pdbid'])                                                      ##  Print PDB structure code
            library = generator.Converter(row['pdbid'],datadir)                                     ##  Converter Class: library
            library.pdb_to_pdbqt(protein=True,ligand=True,native_ligand=True)                       ##  Convert files from pdb to pdbqt (define which ones: protein, ligand, native_ligand)
            library.pdb_to_mol(protein=True,ligand=True,native_ligand=True)                         ##  Convert files from pdb to mol (define which ones)

    if docking_step:
        df_prep = datasets.DatasetPreparation(df)                           #   Generate Molecule list Class - is it neccessery in this case?
        molecules = df_prep.get_molecules('pdbid')                                      ##  Generating molecules list from PDB structure codes

        if 'smina' in docking_programs:
            smina_docking = docking.Smina(datadir)                                      # SMINA Docking Class
            smina_docking.smina_dirs()                                                              ## Generate output dirs for SMINA docking
            smina_matrix = smina_docking.create_smina_matrix(molecules[:1],molecules[:1],no_modes)  ## Generate empty SMINA outputs matrix

        if 'rxdock' in docking_programs:
            rx_docking = docking.RxDock(datadir)
            rx_docking.rxdock_dirs()
            rxdock_matrix = rx_docking.create_rxdock_matrix(molecules[:1],molecules[:1],no_modes)

        for molecule_idx,molecule in enumerate(molecules[:1]):                                      ## Docking Loop for molecules from list generated earlier
            print('Docking '+molecule+' to '+molecule+'. With: \n',docking_programs)                ## Print PDB structure code
            if 'smina' in docking_programs:
                smina_docking_error_number = 0
                while True:                                                                         ## Loop - necessery to generate exactly 100 modes, sometimes with random seed it's generating smaller number of conforms
                    try:
                        smina_docking.smina_files(molecule,molecule,molecule)           ## set variables for specific molecule files
                        smina_docking.smina_docking(no_modes)                           ## smina docking function
                        smina_docking.read_experimental_affinity(df,molecule,molecule)  ## reading experimental affinity data for specific molecule from dataframe
                        smina_docking.read_scoring_function()                                       ## reading scoring function predicted binding affinity from output
                        smina_docking.read_atom_term_function(no_modes)                             ## reading atom terms sf's components from output
                        smina_docking.fill_smina_matrix(molecule_idx,molecule_idx)                  ## fill smina matrix with output datas                      ### MAYBE inserts read_* functions into this one
                    except:
                        smina_docking_error_number+=1
                        print('Smina proposed less modes than expected for docking '+molecule+' to '+molecule+'. '+str(smina_docking_error_number)+'st time.')
                        continue
                    break
            if 'rxdock' in docking_programs:
                rx_docking_error_number = 0
                rx_docking.rxdock_files(molecule,molecule,molecule)
                rx_docking.rxdock_system_preparation()
                while True:
                    try:
                        rx_docking.rxdock_files(molecule,molecule,molecule)
                        #rx_docking.rxdock_docking(no_modes)
                        rx_docking.read_output()
                    except:
                        rx_docking_error_number+=1
                        print('RxDock proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + rx_docking_error_number + 'st time.')
                        continue
                    break

        if 'smina' in docking_programs:
            smina_docking.save_matrix('smina_matrix')                                       ## save SMINA matrix


