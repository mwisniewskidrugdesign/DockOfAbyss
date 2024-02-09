import settings
from pipelines import bdb2020plus, lp_pdbbind
from verify.docking_outputs import diag_smina_verification,diag_gnina_verification,diag_rxdock_verification
from protein_ligand.datasets import DatasetPreparation
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

elif settings.pipeline == 'lp_pdbbind':
    lp_pdbbind.diagonal_pipeline('pdbid',batch_start,batch_end)

else:
    df = pd.read_csv(settings.raw_dataframe)
    df = df[df['CL1'] == True]

    docking_verification = DatasetPreparation(df)
    molecules = docking_verification.get_molecules('pdbid')

    print(molecules)

    diag_smina_verification(molecules)
    diag_rxdock_verification(molecules)
    diag_gnina_verification(molecules)

print('~~~~Fin~~~~~')
