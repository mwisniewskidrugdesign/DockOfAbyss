from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

pd.set_option('display.max_columns', None)

generate_library_step = True
convert_step=True                #ADD IF !!!!!!!!!!!!!!!!!!!!!!!!!! SUCH US DOCKING PROGRAMS LIST
docking_step=True
docking_programs=['smina']

def pipeline(datadir,lp_pdbbind_datadir,lp_pdbbind_df,no_mode):
  #prep DF step for Clear 1 or Clear 2 !!!!
  if generate_library_step:
    '''Generate the library for LP_PDBBIND operations'''
    generator.generate_libraray(datadir)
    workspace = generator.GetDataset(datadir,lp_pdbbind_datadir,lp_pdbbind_df)                                    #create workspace
    workspace.lp_pdbbind()                                                    #copy files to workspace

  if convert_step:
    for index,row in lp_pdbbind_df.iterrows():
      print(str(index) + '.' + row['pdbid'])
      library = generator.Converter(row['pdbid'],datadir)
      library.pdb_to_pdbqt(protein=True,pocket=True,ligand=True,native_ligand=True)
      ### add for other

  if docking_step:
    df_prep = datasets.DatasetPreparation(lp_pdbbind_df)
    molecules = df_prep.get_molecules('pdbid')



    ##### TO DO
    if 'smina' in docking_programs:
      smina_docking = docking.Smina(datadir)
    if 'smina' in docking_programs:
      smina_docking_error_number = 0
      while True:
        try:

