import settings
from pipelines import bdb2020plus, lp_pdbbind
from verify import docking_outputs
import sys
import pandas as pd

pd.set_option('display.max_columns', None)      # Display
settings.init()                                 # inicjalizacja zmiennych odnośnie ścieżek i innych w przyszłości

### SYS argv ###

#
batch_start = int(sys.argv[1])
if batch_start==0:
    batch_start=None
batch_end = int(sys.argv[2])
if batch_end==0:
    batch_end=None

##################

### PIPELINES ###

if settings.pipeline == 'bdb2020':
    bdb2020plus.diagonal_pipeline(settings.datadir,
                                  settings.raw_datadir,
                                  settings.raw_dataframe,
                                  settings.number_of_models)

if settings.pipeline == 'lp_pdbbind':
    lp_pdbbind.diagonal_pipeline('pdbid',batch_start,batch_end)


print('~~~~Fin~~~~~')
