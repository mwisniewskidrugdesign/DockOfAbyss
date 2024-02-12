import os.path
import sys
import settings
import numpy as np
import pandas as pd
from utils.files import smina_output_files_checker, rxdock_output_files_checker, gnina_output_files_checker, diffdock_output_files_checke
from protein_ligand import datasets, generator, docking, matrix
from typing import List

pd.set_option('display.max_columns', None)

def diagonal_pipeline(pdb_id_column: str, #  data frame column with specified pdb_id column
                      batch_start: int,      # if u want to run your pipeline in batches, there u should specify the start index based on dataframe index
                      batch_end: int,        # if u want to run your pipeline in batches, there u should specify the end index based on dataframe index
                      ):
  '''Diagonal pipeline'''

  settings.init()

  df = pd.read_csv(settings.raw_dataframe)
  df = df[df['CL1'] == True]
  df.to_csv(settings.datadir+'/docs/results/dataframe_of_molecules_to_docking.csv')
  df = df[batch_start:batch_end]

  whole_dataset = datasets.DatasetPreparation(df)  # Generate Molecule list Class - is it neccessery in this case?
  molecules = whole_dataset.get_molecules('pdbid')  # Generating molecules list from PDB structure codes

  if 'generate_library' in settings.steps:

    '''Generate the workspace from raw data directory'''
    generator.generate_libraray()
    workspace = generator.GetDataset(pdb_id_column)
    workspace.lp_pdbbind()
  if 'convert' in settings.steps:

    '''Convert all specified in your batch and dataframe molecule files from your data directory'''

    for index,row in df.iterrows():

      print('[CONVERTING]',str(index) / len(df))


      library = generator.Converter(row['pdbid'])

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
  if 'docking' in settings.steps:

    '''Process docking for your dataset'''

    if 'smina' in settings.docking_programs:

      smina_docking = docking.Smina()  #  SMINA Docking Class
      smina_docking.smina_dirs()              #  Generate output dirs for SMINA docking
    if 'rxdock' in settings.docking_programs:

      rx_docking = docking.RxDock()
      rx_docking.rxdock_dirs()
    if 'gnina' in settings.docking_programs:

      gnina_docking = docking.Gnina()
      gnina_docking.gnina_dirs()
    if 'diffdock' in settings.docking_programs:
      diffdock_docking = docking.DiffDock()
      diffdock_docking.diffdock_dirs()


    for molecule_idx, molecule in enumerate(molecules):  #  Docking loop for molecules from generated earlier list

      print(molecule_idx,'/',len(molecules))

      #  Docking Sub-Module

      if 'smina' in settings.docking_programs:

        smina_docking_error_number = 0
        #  Generating input and output files variables
        smina_docking.smina_files(molecule, molecule, molecule)

        #  Loop - necessery to generate exactly <n> modes, sometimes with random seed it's generating smaller number of conforms
        while True:

          #  checking whether output files didn't exists or if output file exists check them
          smina_checker, _ = smina_output_files_checker(smina_docking.log_output_file, smina_docking.pdbqt_output_file, smina_docking.atom_terms_output_file)
          if smina_checker == True:

            #  Smina docking function
            results = smina_docking.smina_docking()

            #  Checking whether output generates proper x number of modes
            _, smina_modes_checker = smina_output_files_checker(smina_docking.log_output_file, smina_docking.pdbqt_output_file, smina_docking.atom_terms_output_file)
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
            print('Output files exist')
            break
      if 'rxdock' in settings.docking_programs:

        rx_docking_error_number = 0
        rx_rmsd_error_number=0
        rx_docking.rxdock_files(molecule, molecule, molecule)
        rx_docking.rxdock_system_preparation()
        rx_cavity_checker,_,_ = rxdock_output_files_checker(cavity_file=rx_docking.cavity1_grd_file)

        while True:
          if rx_cavity_checker == True:
            try:
              _,output_checker,_ = rxdock_output_files_checker(output_file=rx_docking.rx_output+'.sd')
              if output_checker==False:
                rx_docking.rxdock_docking()
            except:
              rx_docking_error_number += 1
              print('RxDock proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + str(rx_docking_error_number) + 'st time.')
              continue
            try:
              _,output_checker,rmsd_checker = rxdock_output_files_checker(output_file=rx_docking.rx_output+'.sd',rmsd_file=rx_docking.rx_output+'_rmsd.sdf')
              if output_checker == True and rmsd_checker == False:
                rx_docking.rxdock_rmsd()
            except:
              rx_rmsd_error_number +=1
              print('RxDock RMSD error')
              continue
          break
      if 'gnina' in settings.docking_programs:

        gnina_docking_error_number = 0
        gnina_docking.gnina_files(molecule, molecule, molecule)

        while True: #  Loop - necessery to generate exactly 50 modes, sometimes with random seed it's generating smaller number of conforms

          gnina_output_checker, gnina_modes_checker, _ = gnina_output_files_checker(log_output_file=gnina_docking.log_output_file,
                                                                                    sdf_gz_output_file=gnina_docking.sdf_gz_output_file,
                                                                                    atom_terms_output_file=gnina_docking.atom_terms_output_file)
          if gnina_output_checker == False or gnina_modes_checker == False:
            gnina_docking.gnina_docking()

            _, gnina_modes_checker, _ = gnina_output_files_checker(log_output_file=gnina_docking.log_output_file)
            if gnina_modes_checker == False:

              gnina_docking_error_number += 1
              print('Gnina proposed less modes than expected for docking ' + molecule + ' to ' + molecule + '. ' + str(gnina_docking_error_number) + 'st time.')
              continue

          gnina_output_checker, _, gnina_rmsd_checker = gnina_output_files_checker(log_output_file=gnina_docking.log_output_file,
                                                                      sdf_gz_output_file=gnina_docking.sdf_gz_output_file,
                                                                      atom_terms_output_file=gnina_docking.atom_terms_output_file,
                                                                      rmsd_output_file=gnina_docking.rmsd_output_file)

          if gnina_output_checker == True and gnina_rmsd_checker==False:
            gnina_docking.gnina_rmsd_calc()
            _, _, gnina_rmsd_checker = gnina_output_files_checker(log_output_file=gnina_docking.log_output_file,
                                                                  sdf_gz_output_file=gnina_docking.sdf_gz_output_file,
                                                                  atom_terms_output_file=gnina_docking.atom_terms_output_file,
                                                                  rmsd_output_file=gnina_docking.rmsd_output_file)

            if gnina_rmsd_checker == False:
              continue
          break
      if 'diffdock' in settings.docking_programs:
        csv_checker = diffdock_output_files_checker()
        print(checker)
        if csv_checker == False:
          for molecule_idx, molecule in enumerate(molecules):
            print(molecule)
            diffdock_docking.diffdock_files(molecule,molecule)
            diffdock_docking.update_diffdock_dataframe()
          diffdock_docking.save_diffdock_dataframe()

    if 'diffdock' in settings.docking_programs:
      diffdock_docking.diffdock_docking()

  if 'matrix' in settings.steps:

    matrix_preparation = generator.Documents()
    docking_output_boolean_dataframe = matrix_preparation.get_diagonal_unification(molecules,molecules)
    docking_output_boolean_dataframe.to_csv(settings.datadir+'/docs/docking_output_booleans.csv')

    docking_output_all_true_dataframe = docking_output_boolean_dataframe[~docking_output_boolean_dataframe[settings.matrix_programs].all(axis=1)]
    docking_output_all_true_dataframe.to_csv(settings.datadir+'/docs/docking_output_all_true.csv')

    matrix_preparation = datasets.DatasetPreparation(docking_output_all_true_dataframe)
    matrix_molecules = matrix_preparation.get_molecules('proteins')
    with open(settings.datadir+'/docs/matrix_molecules.txt', 'w') as file:
      for molecule in matrix_molecules:
        file.write(str(molecule) + '\n')

    if 'smina' in settings.docking_programs:

      smina = matrix.SminaResultsMatrix(matrix_molecules,matrix_molecules,7,df)
      smina_matrix = smina.create_matrix()

      for mol_idx, matrix_molecule in enumerate(matrix_molecules):
        smina.protein = matrix_molecule
        smina.ligand = matrix_molecule
        experimental_affinity = smina.read_experimental_affinity()
        rmsds = smina.read_RMSD()
        predicted_affinities = smina.read_predicted_binding_affinity()
        scoring_components = smina.read_atom_terms()

        for mod_idx in range(settings.number_of_models):

          smina.fill_matrix(mol_idx,mol_idx,mod_idx)

        smina.save_matrix('smina_matrix')
    if 'rxdock' in settings.docking_programs:
      rxdock = matrix.RxDockResultsMatrix(matrix_molecules,matrix_molecules,21,df)
      rxdock_matrix = rxdock.create_matrix()

      for mol_idx, matrix_molecule in enumerate(matrix_molecules):
        rxdock.protein = matrix_molecule
        rxdock.ligand = matrix_molecule
        experimental_affinity = rxdock.read_experimental_affinity()
        rmsds,predicted_affinities,scoring_components = rxdock.read_data()

        for mod_idx in range(settings.number_of_models):

          rxdock.fill_matrix(mol_idx,mol_idx,mod_idx)

        rxdock.save_matrix('rxdock_matrix')
    if 'gnina' in settings.docking_programs:
      print('elo')
