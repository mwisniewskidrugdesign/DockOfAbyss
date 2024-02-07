import pandas as pd

columns = ['complex_name', 'protein_path', 'ligand_description', 'protein_sequence']
diffdock_df = pd.DataFrame(columns=columns)
print(diffdock_df)