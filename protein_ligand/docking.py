import os
import subprocess
import numpy as np
import pandas as pd
import tensorflow as tf
from typing import List
import settings

class Smina:
    def __init__(self):
        #  dirs
        settings.init()
        self.smina_dir = settings.datadir + '/docking_scores/smina'
        self.pdbqt_smina_dir = self.smina_dir + '/pdbqt'
        self.atom_terms_smina_dir = self.smina_dir + '/atom_terms'
        self.logs_smina_dir = self.smina_dir + '/logs'
        #  input files
        self.protein_file=''
        self.ligand_file=''
        self.native_ligand_file=''
        #  output files
        self.pdbqt_output_file = ''
        self.atom_terms_output_file =''
        self.log_output_file = ''
    def smina_dirs(self):
        '''Create Directories for SMINA outputs'''

        if not os.path.exists(self.smina_dir):
            makedir = subprocess.run(['mkdir ' + self.smina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.pdbqt_smina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.atom_terms_smina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.logs_smina_dir], shell=True, capture_output=True, text=True)
    def smina_files(self,protein: str,ligand: str,native_ligand: str):
        '''Specify Input and Output files for SMINA'''

        self.protein_file = settings.datadir+'/protein/pdbqt/'+protein+'_protein.pdbqt'
        self.ligand_file = settings.datadir+'/ligand/pdbqt/'+ligand+'_ligand.pdbqt'
        self.native_ligand_file = settings.datadir+'/native_ligand/pdbqt/'+native_ligand+'_ligand.pdbqt'

        self.pdbqt_output_file = settings.pdbqt_smina_dir + '/' + protein + '_' + ligand + '.pdbqt'
        self.atom_terms_output_file = settings.atom_terms_smina_dir + '/' + protein + '_' + ligand + '_atom_terms.txt'
        self.log_output_file = settings.logs_smina_dir + '/' + protein + '_' + ligand + '.log'
    def smina_docking(self):
        '''Autodock Smina's docking function'''

        smina_command=[settings.smina_tools_dir,
                       '-r',self.protein_file,
                       '-l',self.ligand_file,
                       '--autobox_ligand',self.native_ligand_file,
                       '--autobox_add','8',
                       '--exhaustiveness','32',
                       '--num_modes',str(settings.number_of_models),
                       '-o',self.pdbqt_output_file,
                       '--atom_terms',self.atom_terms_output_file,
                       '--log',self.log_output_file,
                       '--atom_term_data',
                       '--cpu','3',
                       '--min_rmsd_filter','0',
                       '--energy_range','10000']
        docking = subprocess.run(smina_command, shell=False, capture_output=True, text=True)
        print(docking.stderr)
        print(docking.stdout)
class RxDock:
    def __init__(self,system_file='/mnt/evafs/groups/sfglab/mwisniewski/PhD/DockOfAbyss/configs/rxdock_config.prm'):
        settings.init()
        self.rxdock_dir = settings.datadir + '/docking_scores/rxdock'

        self.protein = None
        self.ligand = None
        self.native_ligand = None

        self.protein_file=''
        self.ligand_file=''
        self.native_ligand_file=''

        self.dock_prm_file = '/mnt/evafs/groups/sfglab/mwisniewski/anaconda3/envs/abyss/share/rxdock-2013.1.1_148c5bd1-1/data/scripts/dock.prm'
        self.system_file = system_file          #rx_dock raw system filepath
        self.system_prepared_file = None        #rx_dock prepared system filepath
        self.rx_output = None                   #rx_dock output filepath
    def rxdock_dirs(self):
        if not os.path.exists(self.rxdock_dir):
            makedir = subprocess.run(['mkdir ' + self.rxdock_dir], shell=True, capture_output=True, text=True)
    def rxdock_files(self,protein,ligand,native_ligand):

        self.protein = protein
        self.ligand = ligand
        self.native_ligand = native_ligand

        self.protein_file = settings.datadir + '/protein/mol2/' + protein + '_protein.mol2'
        self.ligand_file = settings.datadir + '/ligand/sdf/' + ligand + '_ligand.sdf'
        self.native_ligand_file = settings.datadir + '/native_ligand/sdf/' + native_ligand + '_ligand.sdf'
        self.system_prepared_file = settings.datadir+'/docs/temp/'+self.protein+'_'+self.ligand+'.prm'
        self.cavity1_grd_file = settings.datadir+'/docs/temp/'+self.protein+'_'+self.ligand+'_cav1.grd'
        self.rx_output = self.rxdock_dir+ '/' + self.protein + '_' + self.ligand
    def rxdock_system_preparation(self):
        with open(self.system_file,'r') as file:
            system_filedata = file.read()

        system_filedata = system_filedata.replace('{title}', self.protein+'_'+self.ligand)                          #define title
        system_filedata = system_filedata.replace('{receptor_file}', self.protein_file)                   #define receptor.mol2 filepath
        system_filedata = system_filedata.replace('{native_ligand_file}', self.native_ligand_file)

        with open(self.system_prepared_file, 'w') as file:
            file.write(system_filedata)
    def rxdock_docking(self):
        os.environ['RBT_HOME'] = settings.datadir
        command = ['rbcavity','-W','-d','-r',self.system_prepared_file]
        result = subprocess.run(command, shell=False, capture_output=True, text=True)
        print('cavity:')
        print(result.stderr)
        print(result.stdout)
        command = ['rbdock','-i',self.ligand_file,'-o',self.rx_output,'-r',self.system_prepared_file,'-p',self.dock_prm_file,'-n',str(settings.number_of_models)]
        result = subprocess.run(command, shell=False, capture_output=True, text=True)
        print('rbdock:')
        print(result.stderr)
        print(result.stdout)
    def rxdock_rmsd(self):
        os.environ['RBT_HOME'] = settings.datadir
        rmsd_output = self.rx_output+'_rmsd.sdf'
        docking_output=self.rx_output+'.sd'
        command = ['sdrmsd','-o',rmsd_output,self.native_ligand_file,docking_output]
        result = subprocess.run(command,shell=False,capture_output=True,text=True)
        print('sdrmsd:')
        print(result.stderr)
        print(result.stdout)
class Gnina:
    def __init__(self):
        self.gnina_dir = settings.datadir+'/docking_scores/gnina'
        self.sdf_gz_gnina_dir = self.gnina_dir+'/sdf_gz'
        self.atom_terms_gnina_dir = self.gnina_dir +'/atom_terms'
        self.logs_gnina_dir = self.gnina_dir+'/logs'
        self.rmsd_gnina_dir = self.gnina_dir+'/rmsd'
        self.protein_file=''
        self.ligand_file=''
        self.native_ligand_file=''
        self.sdf_gz_output_file = ''
        self.atom_terms_output_file =''
        self.log_output_file = ''
        self.rmsd_output_file=''
    def gnina_dirs(self):
        if not os.path.exists(self.gnina_dir):
            makedir = subprocess.run(['mkdir ' + self.gnina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.sdf_gz_gnina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.atom_terms_gnina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.logs_gnina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.rmsd_gnina_dir], shell=True, capture_output=True, text=True)
    def gnina_files(self,protein,ligand,native_ligand):

        self.protein_file = settings.datadir+'/protein/pdb/'+protein+'_protein.pdb'
        self.ligand_file = settings.datadir+'/ligand/sdf/'+ligand+'_ligand.sdf'
        self.native_ligand_file = settings.datadir+'/native_ligand/sdf/'+native_ligand+'_ligand.sdf'

        self.sdf_gz_output_file = self.sdf_gz_gnina_dir + '/' + protein + '_' + ligand + '.sdf.gz'
        self.atom_terms_output_file = self.atom_terms_gnina_dir+'/'+protein+'_'+ligand+'_atom_terms.txt'
        self.log_output_file = self.logs_gnina_dir + '/' + protein + '_' + ligand + '.log'
        self.rmsd_output_file = self.rmsd_gnina_dir+'/'+protein+'_'+ligand+'_rmsd.txt'
    def gnina_docking(self):
        '''GNINA docking'''
        print('Gnina docking')
        gnina_command=['singularity','run','--nv','--bind','/mnt',settings.gnina_container,'gnina',
                       '-r',self.protein_file,
                       '-l',self.ligand_file,
                       '--autobox_ligand',self.native_ligand_file,
                       '--autobox_add','8',
                       '--exhaustiveness','32',
                       '--num_modes',str(settings.number_of_models),
                       '-o',self.sdf_gz_output_file,
                       '--atom_terms',self.atom_terms_output_file,
                       '--log',self.log_output_file,
                       '--atom_term_data',
                       '--cpu','3',
                       '--min_rmsd_filter','0']

        docking = subprocess.run(gnina_command, shell=False, capture_output=True, text=True)

        print(docking.stdout)
        print(docking.stderr)
    def gnina_rmsd_calc(self):
        '''GNINA output RMSD calculation'''
        print('Gnina output RMSD calculation')
        obrms_command=['singularity','run','--nv','--bind','/mnt',settings.gnina_container,'obrms','--firstonly',self.native_ligand_file,self.sdf_gz_output_file]
        rmsd_calculation = subprocess.run(obrms_command, shell=False, capture_output=True, text=True)
        with open(self.rmsd_output_file,'w') as rmsd_file:
            rmsd_file.write(rmsd_calculation.stdout)
        print(rmsd_calculation.stdout)
        print(rmsd_calculation.stderr)
class DiffDock:
    def __init__(self):
        settings.init()

        self.diffdock_dir = settings.datadir + '/docking_scores/diffdock'
        self.diffdock_results_dir = self.diffdock_dir+'/results'
        self.columns=['complex_name','protein_path','ligand_description','protein_sequence']
        self.diffdock_df = pd.DataFrame(columns=self.columns)
        self.csv_diffdock_file = self.diffdock_dir+'/protein_ligand.csv'
        self.complex = None
    def diffdock_dirs(self):
        '''Create Directories for DiffDock outputs'''

        if not os.path.exists(self.diffdock_dir):
            makedir = subprocess.run(['mkdir ' + self.diffdock_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir' + self.diffdock_results_dir], shell=True, capture_output=True, text=True)
    def diffdock_files(self,protein,ligand):
        '''Specify Input and Output files for DiffDock'''

        self.complex = protein+'_'+ligand
        self.protein_file = settings.datadir+'/protein/pdb/'+protein+'_protein.pdb'
        self.ligand_file = settings.datadir+'/ligand/sdf/'+ligand+'_ligand.sdf'
    def update_diffdock_dataframe(self):
        '''Update the dataframe. I have to make it as an update instead of create because of two types of pipelines: diagonal and all vs all.
        Maybe it will be changed in the future.'''

        complex_data = [self.complex,self.protein_file,self.ligand_file,None]
        new_line = pd.DataFrame([complex_data],columns=self.columns)

        self.diffdock_df = pd.concat([self.diffdock_df,new_line],ignore_index=True)
    def save_diffdock_dataframe(self):
        '''Save the dataframe.'''

        self.diffdock_df.to_csv(self.csv_diffdock_file,index=False)
    def files_checker(self):
        csv_checker = os.path.exists(self.diffdock_dir+'/protein_ligand.csv')
        return csv_checker
    def diffdock_docking(self):
        '''This method running basic DiffDock script for predicting the binding of complex'''
        diffdock_command=[
            'python', settings.diffdock_inference_script_path,
            '--protein_ligand_csv',self.diffdock_dir+'/protein_ligand.csv',
            '--out_dir',self.diffdock_results_dir,
            '--inference_steps','20',
            '--samples_per_complex','40',
            '--batch_size','10',
            '--actual_steps','18',
            '--no_final_step_noise'
        ]

        docking = subprocess.run(diffdock_command, shell=False, capture_output=True, text=True)
        print(docking.stderr)
        print(docking.stdout)

