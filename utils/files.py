import os
import pandas as pd
import numpy as np
import sys
import settings


def smina_output_files_checker(log_file,pdbqt_file,atom_terms_file):
    '''
    Checking whether SMINA output files didnt exists or else check them.
    '''

    settings.init()
    if os.path.exists(log_file):
        with open(log_file, 'r') as fp:
            x = len(fp.readlines())
        checker = x != settings.number_of_models + 25 or not (
                os.path.exists(log_file) and os.path.exists(pdbqt_file) and os.path.exists(atom_terms_file))
        modes_checker = x != settings.number_of_models + 25 or not os.path.exists(log_file)
    else:
        checker = True
        modes_checker = True
    return checker, modes_checker
def rxdock_output_files_checker(cavity_file='',output_file='',rmsd_file=''):
    '''
    Checking whether RxDock output files didnt exists or else check them.
    '''

    settings.init()
    cavity_checker = os.path.exists(cavity_file)
    output_checker = os.path.exists(output_file)
    rmsd_checker = os.path.exists(rmsd_file)

    return cavity_checker, output_checker, rmsd_checker
def gnina_output_files_checker(log_output_file='',sdf_gz_output_file='',atom_terms_output_file='',rmsd_output_file=''):
    '''
    Checking whether GNINA output files didnt exists or else check them.
    '''

    settings.init()

    n_of_generated_modes = None
    try:
        with open(log_output_file,'r') as log_file:
            starter = '-----+------------+------------+----------'
            lines = log_file.readlines()
            if any(starter in line for line in lines):
                start_index = next(i for i, line in enumerate(lines) if starter in line)
                n_of_generated_modes = len(lines) - start_index - 1
    except:
        n_of_generated_modes = None

    output_checker = os.path.exists(log_output_file) and os.path.exists(sdf_gz_output_file) and os.path.exists(atom_terms_output_file)
    modes_checker = n_of_generated_modes == settings.number_of_models
    rmsd_checker = os.path.exists(rmsd_output_file)

    return output_checker, modes_checker, rmsd_checker
def diffdock_output_files_checker(output_dir='',csv_file=''):
    '''
    Checking whether DiffDock files did not exists or else check them.
    '''
    settings.init()
    if os.path.exists(output_dir):
        diffdock_files = [file for file in os.listdir(output_dir)
                          if os.path.isfile(os.path.join(output_dir, file)) and os.path.getsize(os.path.join(output_dir, file)) > 0]
        directory_checker = len(diffdock_files) >= 50
    else:
        directory_checker = False
    csv_checker = os.path.exists(csv_file)

    return directory_checker, csv_checker
