import pandas
import pandas as pd

import settings
import os
from typing import List

def smina_verification(sminadir: str,molecule_dockings: List) -> List:

    smina_problems = []

    pdbqt_files = os.listdir(sminadir+'/pdbqt')
    print(pdbqt_files)
    log_files = os.listdir(sminadir+'/logs')
    atom_term_files = os.listdir(sminadir+'/atom_terms')

    molecule_dockings = molecule_dockings

    for molecule_docking in molecule_dockings:
        if (molecule_docking+'.log' not in log_files
                or molecule_docking+'.pdbqt' not in pdbqt_files
                or molecule_docking+'_atom_terms.txt' not in atom_term_files):
            smina_problems.append(molecule_docking)

    return smina_problems

def rxdock_verification(rxdockdir: str, molecule_dockings: List) -> List:

    rxdock_problems=[]

    rxdockdir = os.listdir(rxdockdir)

    for molecule_docking in molecule_dockings:
        if molecule_docking+'.sd' not in rxdockdir:
            rxdock_problems.append(molecule_docking)

    return rxdock_problems

def gnina_verification(gninadir: str,molecule_dockings: List) -> List:

    gnina_problems = []

    sdf_gz_files = os.listdir(gninadir+'/sdf_gz')
    log_files = os.listdir(gninadir+'/logs')
    atom_term_files = os.listdir(gninadir+'/atom_terms')

    molecule_dockings = molecule_dockings

    for molecule_docking in molecule_dockings:
        if (molecule_docking+'.log' not in log_files
                or molecule_docking+'.sdf.gz' not in sdf_gz_files
                or molecule_docking+'_atom_terms.txt' not in atom_term_files):
            gnina_problems.append(molecule_docking)

    return gnina_problems
def save_problems(datadir: str, smina: List, rxdock: List, gnina: List):

    smina_doc = datadir+'/docs/smina_problems.txt'
    rxdock_doc = datadir+'/docs/rxdock_problems.txt'
    gnina_doc = datadir+'/docs/gnina_problems.txt'

    with open(smina_doc, 'w') as output:
        for row in smina:
            output.write(str(row) + '\n')

    with open(rxdock_doc, 'w') as output:
        for row in rxdock:
            output.write(str(row) + '\n')

    with open(gnina_doc, 'w') as output:
        for row in gnina:
            output.write(str(row) + '\n')
def non_docked(datadir: str, dataframe: pd.DataFrame, type='diag') -> List:

    sminadir = datadir + '/docking_scores/smina'
    rxdockdir = datadir + '/docking_scores/rxdock'
    gninadir = datadir + '/docking_scores/gnina'

    molecules = dataframe['pdbid'].to_list()

    if type == 'diag':
        molecule_dockings = [str(molecule)+'_'+str(molecule) for molecule in molecules]
        print('There was: ',str(len(molecule_dockings)),' dockings conducted.')
    elif type == 'all':
        molecule_dockings=[]
        for protein in molecules:
            for ligand in molecules:
                molecule_dockings.append(str(protein)+'_'+str(ligand))
        print('There was: ',str(len(molecule_dockings)),' dockings conducted.')
    print(molecule_dockings)
    smina_problems = smina_verification(sminadir,molecule_dockings)
    rxdock_problems = rxdock_verification(rxdockdir,molecule_dockings)
    gnina_problems = gnina_verification(gninadir,molecule_dockings)

    save_problems(datadir, smina_problems,rxdock_problems,gnina_problems)

