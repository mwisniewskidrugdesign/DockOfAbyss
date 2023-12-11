import numpy as np
import tensorflow as tf

# smina_matrix = np.load('smina_matrix.npy',allow_pickle=True)
# rxdock_matrix = np.load('rxdock_matrix.npy',allow_pickle=True)
#
# print(smina_matrix[:,:,0]) #dla 115
#
# smina_diag = np.diag(smina_matrix[:,:,0])
# print('\n\nSMINA DIAGONAL\n\n')

import os
import sys
import glob
from subprocess import Popen, PIPE
def check_mgl_library():
    if os.path.exists('/Library/MGLTools'):
        MGL_LIB = glob.glob('/Library/MGLTools/*')
        pythonpath = os.path.join(MGL_LIB[0], 'bin', 'pythonsh')
        utilities24 = os.path.join(MGL_LIB[0], 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24')
        prep_receptor_script = os.path.join(utilities24, 'prepare_receptor4.py')
        return (pythonpath, prep_receptor_script)
    else:
        raise Exception('The prepare_receptor4.py script cannot befound  because MGLTools is not installed.')
def execute_receptor_prep(receptor_pdb, output_pdbqt):
    pythonpath, prep_script = check_mgl_library()
    args = (pythonpath, prep_script, receptor_pdb, output_pdbqt)
    cmd = '%s %s -r %s -o %s -A checkhydrogens' % args
    p = Popen(cmd, shell=True, stdout=PIPE)
    stdout, stderr = p.communicate()
    print 'Done!'

if __name__ == '__main__':
    try:
        receptor_pdb, output_pdbqt = sys.argv[1:]
        execute_receptor_prep(receptor_pdb, output_pdbqt)
    except ValueError:
        print '\nUsage: prepare_receptor.py <receptor_pdb> <output_pdbqt>\n'