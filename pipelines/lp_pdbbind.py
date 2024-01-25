import sys

import numpy as np

import settings
from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

pd.set_option('display.max_columns', None)

def diagonal_pipeline(datadir: str, rawdir: str,df: pd.DataFrame,no_modes: int,pdb_id_column: str,batch_start,batch_end,docking_programs=[],steps=[],matrix_type=None):

  generate_library_step = False
  convert_step = False
  docking_step = False
  matrix_step = False

  if 'generate_library' in steps:
    generate_library_step = True
  if 'convert' in steps:
    convert_step = True
  if 'docking' in steps:
    docking_step = True
  if 'matrix' in steps:
    matrix_step = True

  mask = df['CL1'] == True
  df=df[mask]

  if generate_library_step:
    '''Generate the workspace for LP_PDBBIND operations'''
    df = df[batch_start:batch_end]
    generator.generate_libraray(datadir)
    workspace = generator.GetDataset(datadir,
                                     rawdir,
                                     df,
                                     pdb_id_column)
    workspace.lp_pdbbind()

  if convert_step:

    '''Convert all nessesscery files for LP_PDBBIND operations'''

    settings.init()

    df = df[batch_start:batch_end]
    print(settings.to_pdbqt['ligand'])
    for index,row in df.iterrows():

      print(str(index) + '.' + row['pdbid'])

      library = generator.Converter(datadir,row['pdbid'])

      library.pdb_to_pdbqt(protein=settings.to_pdbqt['protein'],
                           pocket=settings.to_pdbqt['pocket'])
      library.pdb_to_mol(protein=settings.to_mol['protein'],
                         pocket=settings.to_mol['pocket'])
      library.mol_to_pdb(ligand=settings.to_pdb['ligand'],
                         native_ligand=settings.to_pdb['native_ligand'])
      library.pdb_to_pdbqt(ligand=settings.to_pdbqt['ligand'],
                           native_ligand=settings.to_pdbqt['native_ligand'])
      library.mol_to_sdf(ligand=settings.to_sdf['ligand'],
                         native_ligand=settings.to_sdf['native_ligand'])

  if docking_step:
    print('Docking Step')
    docking_df = df[batch_start:batch_end]
    df_prep = datasets.DatasetPreparation(docking_df)  # Generate Molecule list Class - is it neccessery in this case?
    molecules = df_prep.get_molecules('pdbid')  ##  Generating molecules list from PDB structure codes

    if 'smina' in docking_programs:

      smina_docking = docking.Smina(datadir)  # SMINA Docking Class
      smina_docking.smina_dirs()  ## Generate output dirs for SMINA docking

    if 'rxdock' in docking_programs:

      rx_docking = docking.RxDock(datadir)
      rx_docking.rxdock_dirs()

    if 'gnina' in docking_programs:

      gnina_docking = docking.Gnina(datadir)
      gnina_docking.gnina_dirs()

    if 'diffdock' in docking_programs:
      diffdock_docking = docking.DiffDock(datadir)
      diffdock_docking.diffdock_dirs()

    for molecule_idx, molecule in enumerate(molecules[:]):  ## Docking Loop for molecules from list generated earlier

      print('Docking ' + molecule + ' to ' + molecule + '. With: \n', docking_programs)  ## Print PDB structure code

      #  Docking Sub-Module

      if 'smina' in docking_programs:

        smina_docking_error_number = 0
        smina_docking.smina_files(molecule, molecule, molecule)

        while True:  ## Loop - necessery to generate exactly 100 modes, sometimes with random seed it's generating smaller number of conforms

          smina_checker = smina_docking.files_checker()

          if smina_checker == True:
            print('\tFiles Checker: Brakuje Plikow lub Powstało ponizej 50 konformacji')
            results = smina_docking.smina_docking(no_modes)  ## smina docking function
            smina_modes_checker = smina_docking.modes_checker()

            if smina_modes_checker == False:
              print('\tModes Checker: Dokowanie zakończone sukcesem')
              break
            else:
              if 'Parse error' in results.stderr:
                print('\tDocking Error')
                break
              smina_docking_error_number +=1
              print('Smina proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + str(smina_docking_error_number) + 'st time.')
              continue
          else:
            break
      if 'rxdock' in docking_programs:

        rx_docking_error_number = 0
        rx_docking.rxdock_files(molecule, molecule, molecule)
        rx_docking.rxdock_system_preparation()
        rxdock_checker = rx_docking.files_checker()
        print('RxDock Checker: ',rxdock_checker)

        while True:
          if rxdock_checker == True:
            try:
              rx_docking.rxdock_docking(no_modes)
            except:
              rx_docking_error_number += 1
              print('RxDock proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + str(rx_docking_error_number) + 'st time.')
              continue
          break

        #rx_docking.read_output(molecule_idx, molecule_idx)
      if 'gnina' in docking_programs:
        gnina_docking_error_number = 0
        gnina_docking.gnina_files(molecule, molecule, molecule)
        gnina_checker = gnina_docking.gnina_files_checker()
        while True:  ## Loop - necessery to generate exactly 50 modes, sometimes with random seed it's generating smaller number of conforms
          if gnina_checker == True:
            try:
              gnina_docking.gnina_docking(no_modes)
            except:
              gnina_docking_error_number += 1
              print('Smina proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + str(gnina_docking_error_number) + 'st time.')
              continue
          break
      if 'diffdock' in docking_programs:
        diffdock_docking_error_number = 0
        diffdock_docking.diffdock_files(molecule, molecule, molecule)
        diffdock_checker = diffdock_docking.diffdock_files_checker()
        while True:  ## Loop - necessery to generate exactly 50 modes, sometimes with random seed it's generating smaller number of conforms
          if diffdock_checker == True:
            try:
              diffdock_docking.diffdock_docking(no_modes)
            except:
              diffdock_docking_error_number += 1
              print('Smina proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + str(diffdock_docking_error_number) + 'st time.')
              continue
          break

  if matrix_step:

    df_prep = datasets.DatasetPreparation(df)  # Generate Molecule list Class - is it neccessery in this case?
    molecules = df_prep.get_molecules('pdbid')  ##  Generating molecules list from PDB structure codes

    if 'smina' in docking_programs:

      smina_matrix = docking.Smina(datadir)
      smina_matrix.create_smina_matrix(molecules,molecules,no_modes)

      for molecule_idx, molecule in enumerate(molecules):
        print(molecule)
        smina_matrix.smina_files(molecule, molecule, molecule)
        checker = smina_matrix.files_checker()
        if checker == False:
          if matrix_type == 'scoring':
            smina_matrix.read_experimental_affinity(df, molecule,molecule)  ## reading experimental affinity data for specific molecule from dataframe
            smina_matrix.read_scoring_function()  ## reading scoring function predicted binding affinity from output
            smina_matrix.read_atom_term_function(no_modes)  ## reading atom terms sf's components from output
            smina_matrix.fill_scoring_matrix(molecule_idx,molecule_idx)  ## fill smina matrix with output datas
          elif matrix_type == 'rmsd':
            smina_matrix.read_rmsd()
            smina_matrix.fill_rmsd_matrix(molecule_idx,molecule_idx)
        else:
          if matrix_type == 'scoring':
            smina_matrix.read_experimental_affinity(df,molecule,molecule)
            smina_matrix.predicted_binding_affinity = [0.] * no_modes
            smina_matrix.sf_components = [[0. for i in range(no_modes)]] * 5
            smina_matrix.fill_scoring_matrix(molecule_idx,molecule_idx)
          elif matrix_type=='rmsd':
            smina_matrix.lb_rmsds=[np.NaN] * no_modes
            smina_matrix.ub_rmsds=[np.NaN] * no_modes
            smina_matrix.fill_rmsd_matrix(molecule_idx,molecule_idx)
      if matrix_type == 'scoring':
        smina_matrix.save_matrix('smina_scoring_matrix')
      elif matrix_type =='rmsd':
        smina_matrix.save_matrix('smina_rmsd_matrix')
    if 'rxdock' in docking_programs:
      rxdock_matrix = docking.RxDock(datadir)
      rxdock_matrix.create_rxdock_matrix(molecules,molecules,no_modes)

      for molecule_idx, molecule in enumerate(molecules):
        print(molecule)
        rxdock_matrix.rxdock_files(molecule, molecule, molecule)
        checker = rxdock_matrix.matrix_file_checker()

        if checker == True:
          if matrix_type == 'scoring':
            rxdock_matrix.read_scoring_output(molecule_idx,molecule_idx)
          elif matrix_type == 'rmsd':
            rmsd_checker = rxdock_matrix.rmsd_file_checker()
            if rmsd_checker == False:
              rxdock_matrix.rxdock_rmsd()
            rxdock_matrix.read_rmsd_output(molecule_idx,molecule_idx)
        elif checker == False:
          if matrix_type == 'scoring':
            for mode_idx in range(no_modes):
              rxdock_matrix.mode_values = tf.convert_to_tensor(np.array([np.NaN] * 300))
              rxdock_matrix.fill_rxdock_matrix(molecule_idx,molecule_idx,mode_idx)
          elif matrix_type == 'rmsd':
            rxdock_matrix.rmsd=np.NaN
            for mode_idx in range(no_modes):
              rxdock_matrix.fill_rxdock_matrix(molecule_idx,molecule,mode_idx)
      if matrix_type == 'scoring':
        rxdock_matrix.save_matrix('rxdock_scoring_matrix')
      elif matrix_type == 'rmsd':
        rxdock_matrix.save_matrix('rxdock_rmsd_matrix')
    if 'gnina' in docking_programs:
      sys.exit()
