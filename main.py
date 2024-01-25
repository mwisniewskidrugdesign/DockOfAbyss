import settings
from pipelines import bdb2020plus, lp_pdbbind
from verify import docking_outputs
import sys
import pandas as pd

pd.set_option('display.max_columns', None)      # wyświetlanie df
settings.init()                                 # inicjalizacja zmiennych odnośnie ścieżek i innych w przyszłości

### SYS argv ###

#
batch_start = int(sys.argv[1])
if batch_start==0:
    batch_start=None
batch_end = int(sys.argv[2])
if batch_end==0:
    batch_end=None

# programs = list(sys.argv[4].split(','))
# steps = list(sys.argv[5].split(','))

###########

### DIRS ###

if settings.station == 'inka':
    datadir='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/bdb2020plus'
    bdb2020plus_datadir = ''

elif settings.station == 'eden':
    datadir='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/lp_pdbbind'
    bdb2020plus_datadir=''
    lp_pdbbind_dir='/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/raw/PDBbind'

############

### DATAFRAMES ###

if settings.station == 'inka':
    bdb2020plus_df = pd.read_csv('')
    lp_pdbbind_df = pd.read_csv('')

elif settings.station == 'eden':
    bdb2020plus_df = pd.read_csv('/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/BDB2020+.csv')
    lp_pdbbind_df = pd.read_csv('/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv')

##################

### PIPELINES ###

bdb2020plus_pipeline = False
mpro_pipeline = False
egfr_pipeline = False
lp_pdbbind_pipeline = False
verification_pipeline = False

if settings.pipeline == 'bdb2020':
    bdb2020plus_pipeline = True

elif settings.pipeline == 'lp_pdbbind':
    lp_pdbbind_pipeline = True

elif settings.pipeline == 'verification':
    verification_pipeline = True

#############

### PIPELINES ###

if bdb2020plus_pipeline:
    bdb2020plus.diagonal_pipeline(datadir,
                                  bdb2020plus_datadir,
                                  bdb2020plus_df,
                                  settings.number_of_models)

if lp_pdbbind_pipeline:
    lp_pdbbind.diagonal_pipeline(datadir,
                                 lp_pdbbind_dir,
                                 lp_pdbbind_df,
                                 settings.number_of_models,
                                 'pdbid',
                                 batch_start,batch_end,
                                 settings.docking_programs,
                                 settings.steps,
                                 settings.matrix_type)

if verification_pipeline:
    docking_outputs.non_docked(datadir,
                               lp_pdbbind_df,
                               type='diag')

print('~~~~Fin~~~~~')
