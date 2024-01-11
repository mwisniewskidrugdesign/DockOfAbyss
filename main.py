import settings
from pipelines import bdb2020plus, lp_pdbbind
from verify import docking_outputs
import sys
import pandas as pd

pd.set_option('display.max_columns', None)      # wyświetlanie df
settings.init()                                 # inicjalizacja zmiennych odnośnie ścieżek i innych w przyszłości

### SYS argv ###

pipeline = sys.argv[1]

batch_start = int(sys.argv[2])
if batch_start==0:
    batch_start=None

batch_end = int(sys.argv[3])
if batch_end==0:
    batch_end=None

programs = list(sys.argv[4].split(','))
steps = list(sys.argv[5].split(','))

###########

### DIRS ###

if settings.station == 'inka':
    datadir='/mnt/raid/mwisniewski/PhD/data/bdb2020plus'
    bdb2020plus_datadir = '/mnt/raid/mwisniewski/PhD/data/LP-PDBBind/dataset/BDB2020+/dataset'

elif settings.station == 'eden':
    datadir='/home2/sfglab/mwisniewski/PhD/data/lp_pdbbind'
    bdb2020plus_datadir='/home2/sfglab/mwisniewski/PhD/data/raw/bdb2020plus'
    lp_pdbbind_dir='/home2/sfglab/mwisniewski/PhD/data/raw/PDBbind'

############

### DATAFRAMES ###

if settings.station == 'inka':
    bdb2020plus_df = pd.read_csv('/mnt/raid/mwisniewski/PhD/data/LP-PDBBind/dataset/BDB2020+/BDB2020+.csv')
    lp_pdbbind_df = pd.read_csv('/home2/sfglab/mwisniewski/PhD/data/others/PDBbind')

elif settings.station == 'eden':
    bdb2020plus_df = pd.read_csv('/home2/sfglab/mwisniewski/PhD/data/dataframes/BDB2020+.csv')
    lp_pdbbind_df = pd.read_csv('/home2/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv')

##################

### PIPELINES ###

bdb2020plus_pipeline = False
mpro_pipeline = False
egfr_pipeline = False
lp_pdbbind_pipeline = False
verification_pipeline = False

if pipeline == 'bdb2020':
    bdb2020plus_pipeline = True
elif pipeline == 'lp_pdbbind':
    lp_pdbbind_pipeline = True
elif pipeline == 'verification':
    verification_pipeline = True

#############

### PIPELINES ###

if bdb2020plus_pipeline:
    bdb2020plus.diagonal_pipeline(datadir, bdb2020plus_datadir, bdb2020plus_df,50)

if lp_pdbbind_pipeline:
    lp_pdbbind.diagonal_pipeline(datadir,lp_pdbbind_dir,lp_pdbbind_df,50,'pdbid',batch_start,batch_end,programs,steps) #batch_start,batch_end

if verification_pipeline:
    docking_outputs.non_docked(datadir,lp_pdbbind_df,type='diag')
