import pandas as pd
import settings
import os
from typing import List

settings.init()

def diag_smina_verification(molecules) -> List:

    settings.init()

    pdbqt_dir = settings.datadir+'dockings_scores/smina/pdbqt'
    log_dir = settings.datadir+'dockings_scores/smina/logs'
    atom_terms_dir = settings.datadir+'dockings_scores/smina/atom_terms'

    pdbqt_booleans = []
    log_booleans = []
    atom_terms_booleans = []
    for molecule in molecules:
        complex = molecule+'_'+molecule
        pdbqt_file = pdbqt_dir + '/' + complex + '.pdbqt'
        log_file = log_dir + '/' + complex + '.log'
        atom_terms_file = atom_terms_dir+'/'+complex+'_atom_terms.txt'

        pdbqt_boolean = True if os.path.exists(pdbqt_file) and os.path.getsize(pdbqt_file) > 0
        log_boolean = True if os.path.exists(log_file) and os.path.getsize(log_file) > 0
        atom_terms_boolean = True if os.path.exists(atom_terms_file) and os.path.getsize(atom_terms_file) > 0


        pdbqt_booleans.append(pdbqt_boolean)
        log_booleans.append(log_boolean)
        atom_terms_booleans.append(atom_terms_boolean)

    dict={'complex_id':molecules,
          'smina_pdbqt_boolean':pdbqt_booleans,
          'smina_log_boolean':log_booleans,
          'smina_atom_terms_boolean': atom_terms_booleans,
          }
    }
    df = pd.DataFrame(dict)
    df.to_csv(settings.datadir+'/docs/smina_booleans.csv')
    return df
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
        for index,row in enumerate(smina):
            output.write(str(index+1)+'.'+str(row) + '\n')

    with open(rxdock_doc, 'w') as output:
        for index,row in enumerate(rxdock):
            output.write(str(index+1)+'.'+str(row) + '\n')

    with open(gnina_doc, 'w') as output:
        for index,row in enumerate(gnina):
            output.write(str(index+1)+'.'+str(row) + '\n')
def non_docked(datadir: str, dataframe: pd.DataFrame, type='diag') -> List:
    mask = dataframe['CL1'] == True
    dataframe = dataframe[mask]
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
    smina_problems = smina_verification(sminadir,molecule_dockings)
    rxdock_problems = rxdock_verification(rxdockdir,molecule_dockings)
    gnina_problems = gnina_verification(gninadir,molecule_dockings)

    save_problems(datadir, smina_problems,rxdock_problems,gnina_problems)
