from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

pd.set_option('display.max_columns', None)

generate_library_step = True
convert_step=True                #ADD IF !!!!!!!!!!!!!!!!!!!!!!!!!! SUCH US DOCKING PROGRAMS LIST
docking_step=False
docking_programs=['smina','rxdock']

def diagonal_pipeline(datadir,rawdir,df,no_modes,pdb_id_column,batch_start,batch_end):
  #prep DF step for Clear 1 or Clear 2 !!!!
  mask = df['CL1'] == True
  df=df[mask]
  print(df)
  if generate_library_step:
    '''Generate the library for LP_PDBBIND operations'''
    generator.generate_libraray(datadir)
    workspace = generator.GetDataset(datadir,rawdir,df,pdb_id_column)                                    #create workspace
    workspace.lp_pdbbind()                                                    #copy files to workspace

  if convert_step:
    df=df[batch_start:batch_end] #to edit
    for index,row in df.iterrows():
      print(str(index) + '.' + row['pdbid'])
      library = generator.Converter(row['pdbid'],datadir)
      #library.pdb_to_pdbqt(protein=True,pocket=True,ligand=False,native_ligand=False)
      #library.pdb_to_mol(protein=True, pocket=True, ligand=False, native_ligand=False)
      library.mol_to_pdb(protein=False,pocket=False,ligand=True,native_ligand=True)
      library.pdb_to_pdbqt(protein=False,pocket=False,ligand=True,native_ligand=True)

      ### add for other

  if docking_step:
    df_prep = datasets.DatasetPreparation(df)  # Generate Molecule list Class - is it neccessery in this case?
    molecules = df_prep.get_molecules('pdbid')[batch_start:batch_end]  ##  Generating molecules list from PDB structure codes

    if 'smina' in docking_programs:
      smina_docking = docking.Smina(datadir)  # SMINA Docking Class
      smina_docking.smina_dirs()  ## Generate output dirs for SMINA docking
      smina_matrix = smina_docking.create_smina_matrix(molecules[:], molecules[:],
                                                       no_modes)  ## Generate empty SMINA outputs matrix

    if 'rxdock' in docking_programs:
      rx_docking = docking.RxDock(datadir)
      rx_docking.rxdock_dirs()
      rxdock_matrix = rx_docking.create_rxdock_matrix(molecules[:], molecules[:], no_modes)

    for molecule_idx, molecule in enumerate(molecules[:]):  ## Docking Loop for molecules from list generated earlier
      print('Docking ' + molecule + ' to ' + molecule + '. With: \n', docking_programs)  ## Print PDB structure code
      if 'smina' in docking_programs:
        smina_docking_error_number = 0
        while True:  ## Loop - necessery to generate exactly 100 modes, sometimes with random seed it's generating smaller number of conforms
          try:
            smina_docking.smina_files(molecule, molecule, molecule)  ## set variables for specific molecule files
            smina_docking.smina_docking(no_modes)  ## smina docking function
            smina_docking.read_experimental_affinity(df, molecule,
                                                     molecule)  ## reading experimental affinity data for specific molecule from dataframe
            smina_docking.read_scoring_function()  ## reading scoring function predicted binding affinity from output
            smina_docking.read_atom_term_function(no_modes)  ## reading atom terms sf's components from output
            smina_docking.fill_smina_matrix(molecule_idx,
                                            molecule_idx)  ## fill smina matrix with output datas                      ### MAYBE inserts read_* functions into this one
          except:
            smina_docking_error_number += 1
            print('Smina proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + str(
              smina_docking_error_number) + 'st time.')
            continue
          break
      if 'rxdock' in docking_programs:
        rx_docking_error_number = 0
        rx_docking.rxdock_files(molecule, molecule, molecule)
        rx_docking.rxdock_system_preparation()
        while True:
          try:
            rx_docking.rxdock_files(molecule, molecule, molecule)
            rx_docking.rxdock_docking(no_modes)
            rx_docking.read_output(molecule_idx, molecule_idx)
          except:
            rx_docking_error_number += 1
            print(
              'RxDock proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + rx_docking_error_number + 'st time.')
            continue
          break

    if 'smina' in docking_programs:
      smina_docking.save_matrix('smina_matrix')  ## save SMINA matrix
    if 'rxdock' in docking_programs:
      rx_docking.save_matrix('rxdock_matrix')
