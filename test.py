import os

path='LP_PDBBind.csv'
if os.path.exists(path):
    with open(path, 'r') as fp:
        x = len(fp.readlines())

checker = not os.path.exists('LP_PDBBind.csv') or x != 75
print(checker)