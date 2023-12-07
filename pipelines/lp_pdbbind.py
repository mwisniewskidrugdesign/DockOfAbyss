from protein_ligand import generator
from protein_ligand import docking
from protein_ligand import datasets
import pandas as pd

pd.set_option('display.max_columns', None)

generate_library_step = True
convert_step=True                #ADD IF !!!!!!!!!!!!!!!!!!!!!!!!!! SUCH US DOCKING PROGRAMS LIST
docking_step=True

def pipeline(datadir,lp_pdbbind_datadir,lp_pdbbind_df,no_mode):
  if generate_library_step:
    '''Generate the library for LP_PDBBIND operations'''
    generator.generate_libraray(datadir)
    workspace = generator.GetDataset(datadir,lp_pdbbind_datadir,lp_pdbbind_df)                                    #create workspace
    workspace.lp_pdbbind()                                                    #copy files to workspace
      
