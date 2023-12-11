import numpy as np
import tensorflow as tf

smina_matrix = np.load('smina_matrix.npy',allow_pickle=True)
rxdock_matrix = np.load('rxdock_matrix.npy',allow_pickle=True)

print(smina_matrix[:,:,0]) #dla 115

smina_diag = np.diag(smina_matrix[:,:,0])
print('\n\nSMINA DIAGONAL\n\n')
