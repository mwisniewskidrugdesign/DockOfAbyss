from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

pd.set_option('display.max_columns', None)

generate_library_step = True
convert_step=True                #ADD IF !!!!!!!!!!!!!!!!!!!!!!!!!! SUCH US DOCKING PROGRAMS LIST
docking_step=True
docking_programs=['smina','rxdock']

def pipeline(datadir,bdb2020plus_datadir,bdb2020plus_df,no_modes):

    if generate_library_step:
        '''Generate the library for bdb2020plus operations'''
        generator.generate_libraray(datadir)
        workspace = generator.GetDataset(datadir,bdb2020plus_df)                                        #create workspace
        workspace.bdb2020plus(bdb2020plus_datadir)                                          #copy files to workspace

    if convert_step:
        for index,row in bdb2020plus_df.iterrows():                                                     #iteruje po dataframe
            print(str(index)+'.'+row['pdbid'])
            library = generator.Converter(row['pdbid'],datadir)                             #Klasa konwertera
            library.pdb_to_pdbqt(protein=True,ligand=True,native_ligand=True)                    #konwertuje do pdbqt
            library.pdb_to_mol(protein=True,ligand=True,native_ligand=True)

    if docking_step:
        bdb2020plus_df_prep = datasets.DatasetPreparation(bdb2020plus_df)
        molecules = bdb2020plus_df_prep.get_molecules('pdbid')           #Generate list of molecules
        print(molecules)

        if 'smina' in docking_programs:
            bdb2020plus_smina_docking = docking.Smina(datadir)
            bdb2020plus_smina_docking.smina_dirs()
            smina_matrix = bdb2020plus_smina_docking.create_smina_matrix(molecules[:],molecules[:],no_modes)

        if 'rxdock' in docking_programs:
            bdb2020plus_rx_docking = docking.RxDock(datadir)
            bdb2020plus_rx_docking.rxdock_dirs()
            rxdock_matrix = None

        for molecule_idx,molecule in enumerate(molecules[:]):
            print('Docking '+molecule+' to '+molecule+'. With: \n',docking_programs)
            if 'smina' in docking_programs:
                smina_docking_error_number = 0
                while True:
                    print(molecule)
                    try:
                        bdb2020plus_smina_docking.smina_files(molecule,molecule,molecule)
                        bdb2020plus_smina_docking.smina_docking(no_modes)
                        bdb2020plus_smina_docking.read_experimental_affinity(bdb2020plus_df,molecule,molecule)
                        bdb2020plus_smina_docking.read_scoring_function()
                        bdb2020plus_smina_docking.read_atom_term_function(no_modes)
                        bdb2020plus_smina_docking.fill_smina_matrix(molecule_idx,molecule_idx)
                    except:
                        smina_docking_error_number+=1
                        print('Smina proposed less modes than expected for docking '+molecule+' to '+molecule+'. '+str(smina_docking_error_number)+'st time.')
                        continue
                    break
            if 'rxdock' in docking_programs:
                rx_docking_error_number = 0
                #while True:
                    #try:
                        #print(elo)
                    #except:
                        #rx_docking_error_number+=1
                        #print('RxDock proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + rx_docking_error_number + 'st time.')

        if 'smina' in docking_programs:
            bdb2020plus_smina_docking.save_matrix('smina_matrix')


