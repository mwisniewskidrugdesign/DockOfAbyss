import subprocess
import os
import sys

os.environ['RBT_HOME'] = '/mnt/raid/mwisniewski/PhD/data'
print(os.environ['RBT_HOME'])
inp = os.environ['RBT_HOME']+'/bdb2020plus/ligand/sdf/6TRX_ligand.sdf'
output =os.environ['RBT_HOME']+'/bdb2020plus/docs/6TRX-6TRX'
system = os.environ['RBT_HOME']+'/bdb2020plus/docs/temp/6TRX-6TRX.prm'
command= 'rbdock -i %s -o %s -r %s -p dock.prm -n 100' % (inp, output, system)

result = subprocess.run([command],shell=True, capture_output=True,text=True)
print(result)
print(result.stdout)
print('error: ')
print(result.stderr)