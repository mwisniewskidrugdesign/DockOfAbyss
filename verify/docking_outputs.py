import pandas as pd
import settings
import os
from typing import List

settings.init()

def diag_smina_verification(molecules: List) -> pd.DataFrame:

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

        pdbqt_boolean = True if (os.path.exists(pdbqt_file) and os.path.getsize(pdbqt_file) > 0) else False
        log_boolean = True if (os.path.exists(log_file) and os.path.getsize(log_file) > 0) else False
        atom_terms_boolean = True if (os.path.exists(atom_terms_file) and os.path.getsize(atom_terms_file) > 0) else False


        pdbqt_booleans.append(pdbqt_boolean)
        log_booleans.append(log_boolean)
        atom_terms_booleans.append(atom_terms_boolean)

    dict_boolean={'complex_id':molecules,
                  'smina_pdbqt_boolean':pdbqt_booleans,
                  'smina_log_boolean':log_booleans,
                  'smina_atom_terms_boolean': atom_terms_booleans
                  }
    df = pd.DataFrame(dict_boolean)
    df.to_csv(settings.datadir+'/docs/results/smina_booleans.csv')
    return df
def diag_rxdock_verification(molecules: List) -> pd.DataFrame:

    settings.init()

    rxdock_dir = settings.datadir+'dockings_scores/smina/pdbqt'

    rxdock_sd_booleans = []
    rxdock_sdf_booleans = []

    for molecule in molecules:
        complex = molecule+'_'+molecule
        sdf_file = rxdock_dir + '/' + complex + '_rmsd.sdf'
        sd_file = rxdock_dir + '/' + complex + '.sd'

        sdf_booleans = True if (os.path.exists(sdf_file) and os.path.getsize(sdf_file) > 0) else False
        sd_booleans = True if (os.path.exists(sd_file) and os.path.getsize(sd_file) > 0) else False


        rxdock_sdf_booleans.append(sdf_booleans)
        rxdock_sd_booleans.append(sd_booleans)

    dict_boolean={'complex_id':molecules,
                  'rxdock_sd_boolean':rxdock_sd_booleans,
                  'rxdock_sdf_boolean':rxdock_sdf_booleans
                  }
    df = pd.DataFrame(dict_boolean)
    df.to_csv(settings.datadir+'/docs/results/rxdock_booleans.csv')
    return df

def diag_gnina_verification(molecules: List) -> pd.DataFrame:
    settings.init()

    sdf_gz_dir = settings.datadir + 'dockings_scores/gnina/sdf_gz'
    rmsd_dir = settings.datadir + 'docking_scores/gnina/rmsd'
    log_dir = settings.datadir + 'dockings_scores/gnina/logs'
    atom_terms_dir = settings.datadir + 'dockings_scores/gnina/atom_terms'

    sdf_gz_booleans = []
    rmsd_booleans = []
    log_booleans = []
    atom_terms_booleans = []

    for molecule in molecules:

        complex = molecule + '_' + molecule
        sdf_gz_file = sdf_gz_dir + '/' + complex + '.pdbqt'
        log_file = log_dir + '/' + complex + '.log'
        atom_terms_file = atom_terms_dir + '/' + complex + '_atom_terms.txt'
        rmsd_file = rmsd_dir + '/' + complex + '_rmsd.txt'

        sdf_gz_boolean = True if (os.path.exists(sdf_gz_file) and os.path.getsize(sdf_gz_file) > 0) else False
        log_boolean = True if (os.path.exists(log_file) and os.path.getsize(log_file) > 0) else False
        atom_terms_boolean = True if (os.path.exists(atom_terms_file) and os.path.getsize(atom_terms_file) > 0) else False
        rmsd_boolean = True if (os.path.exists(rmsd_file) and os.path.getsize(rmsd_file) > 0) else False

        sdf_gz_booleans.append(sdf_gz_boolean)
        log_booleans.append(log_boolean)
        atom_terms_booleans.append(atom_terms_boolean)
        rmsd_booleans.append(rmsd_boolean)

    dict_boolean = {'complex_id': molecules,
                    'gnina_sdf_gz_boolean': pdbqt_booleans,
                    'gnina_log_boolean': log_booleans,
                    'gnina_atom_terms_boolean': atom_terms_booleans,
                    'gnina_rmsd_boolean':rmsd_booleans
                    }
    df = pd.DataFrame(dict_boolean)
    df.to_csv(settings.datadir + '/docs/results/gnina_booleans.csv')
    return df

