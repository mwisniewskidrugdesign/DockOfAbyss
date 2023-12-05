import settings
from pipelines import bdb2020plus

import pandas as pd


pd.set_option('display.max_columns', None)      # wyświetlanie df
settings.init()                                 # inicjalizacja zmiennych odnośnie ścieżek i innych w przyszłości

### DIRS ###
# na to osobny plik z configiem
datadir='/mnt/raid/mwisniewski/PhD/data/bdb2020plus'
bdb2020plus_datadir = '/mnt/raid/mwisniewski/PhD/data/LP-PDBBind/dataset/BDB2020+/dataset'
############

### DATAFRAMES ###
# na to osobny plik z configiem
bdb2020plus_df = pd.read_csv('/mnt/raid/mwisniewski/PhD/data/LP-PDBBind/dataset/BDB2020+/BDB2020+.csv')
##################

### PIPELINES ###
# na to zrobię $1 argument w prompt
bdb2020plus_pipeline = True
mpro_pipeline = False
egfr_pipeline = False
lp_pdbbind_pipeline = False

#############

### PIPELINES ###

if bdb2020plus_pipeline:
    bdb2020plus.pipeline(datadir, bdb2020plus_datadir, bdb2020plus_df,50)


