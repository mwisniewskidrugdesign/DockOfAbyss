# Read the file and extract RMSD values
file_path = "test/3vdb_output.sdf"  # Replace with the actual path to your file
import pandas as pd
pd.set_option('display.max_columns', None)      # Display

# List of keys to extract from the file
# List of keys to extract from the file
keys_to_extract = [
    '<SCORE>', '<SCORE.INTER>', '<SCORE.INTER.CONST>', '<SCORE.INTER.POLAR>',
    '<SCORE.INTER.REPUL>', '<SCORE.INTER.ROT>', '<SCORE.INTER.VDW>',
    '<SCORE.INTER.norm>', '<SCORE.INTRA>', '<SCORE.INTRA.DIHEDRAL>',
    '<SCORE.INTRA.DIHEDRAL.0>', '<SCORE.INTRA.POLAR>', '<SCORE.INTRA.POLAR.0>',
    '<SCORE.INTRA.REPUL>', '<SCORE.INTRA.REPUL.0>', '<SCORE.INTRA.VDW>',
    '<SCORE.INTRA.VDW.0>', '<SCORE.INTRA.norm>', '<SCORE.RESTR>',
    '<SCORE.RESTR.CAVITY>', '<SCORE.RESTR.norm>', '<SCORE.SYSTEM>',
    '<SCORE.SYSTEM.CONST>', '<SCORE.SYSTEM.DIHEDRAL>', '<SCORE.SYSTEM.norm>',
    '<SCORE.HEAVY>', '<SCORE.norm>'
]

# Dictionary to store extracted values
with open(file_path, 'r') as output_file:
    list_of_modes = output_file.read().split('$$$$')[:-1]
    modes=[]
    for mode_index, mode in enumerate(list_of_modes):
        one_model_values = ['SCORE']+mode.split('>  <SCORE>')[-1].split('\n')[1:]
        predicted_binding_affinity = one_model_values[0]
        one_model_values = one_model_values[:]
        keys=one_model_values[::3]
        values=one_model_values[1::3]
        #keys, values = zip(one_model_values[::3],one_model_values[1::3])
        #values = one_model_values[1::3]
        keys = [key.replace('>  <','').replace('>','') for key in keys]
        together = dict(zip(keys,values))
        modes.append(together)

    unique = set(key for mode in modes for key in mode)
    fill_dicts = {key: 0 for key in unique}
    for mode in modes:
        for key in unique:
            mode.setdefault(key,0)
    df = pd.DataFrame(modes)
    print(df[df['SCORE']=='-32.3094']['SCORE.SYSTEM'])
    print(df)

    rmsds=df['RMSD'].tolist()
    df_components = df.drop(['RMSD', 'SCORE'], axis=1)

    scoring_components_per_complex = [row.tolist() for _, row in df_components.iterrows()]
    print(len(scoring_components_per_complex))

