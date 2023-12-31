import os
import subprocess
import numpy as np
import tensorflow as tf

import settings

class Smina:
    def __init__(self,datadir):
        self.datadir = datadir
        self.smina_dir = datadir + '/docking_scores/smina'
        self.pdbqt_smina_dir = self.smina_dir + '/pdbqt'
        self.atom_terms_smina_dir = self.smina_dir + '/atom_terms'
        self.logs_smina_dir = self.smina_dir + '/logs'
        self.protein_file=''
        self.ligand_file=''
        self.native_ligand_file=''
        self.pdbqt_output_file = ''
        self.atom_terms_output_file =''
        self.log_output_file = ''
        self.matrix = None
        self.experimental_affinity = None
        self.predicted_binding_affinity = None
        self.sf_components=None
    def smina_dirs(self):
        if not os.path.exists(self.smina_dir):
            makedir = subprocess.run(['mkdir ' + self.smina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.pdbqt_smina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.atom_terms_smina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.logs_smina_dir], shell=True, capture_output=True, text=True)
    def smina_files(self,protein,ligand,native_ligand):

        self.protein_file = self.datadir+'/protein/pdbqt/'+protein+'_protein.pdbqt'
        self.ligand_file = self.datadir+'/ligand/pdbqt/'+ligand+'_ligand.pdbqt'
        self.native_ligand_file = self.datadir+'/native_ligand/pdbqt/'+native_ligand+'_ligand.pdbqt'

        self.pdbqt_output_file = self.pdbqt_smina_dir + '/' + protein + '_' + ligand + '.pdbqt'
        self.atom_terms_output_file = self.atom_terms_smina_dir + '/' + protein + '_' + ligand + '_atom_terms.txt'
        self.log_output_file = self.logs_smina_dir + '/' + protein + '_' + ligand + '.log'
    def files_checker(self):
        if os.path.exists(self.log_output_file):
            with open(self.log_output_file, 'r') as fp:
                x = len(fp.readlines())
                print(x)
        checker = not os.path.exists(self.log_output_file) or x != 75
        return checker
    def smina_docking(self,no_modes):
        smina_command=[settings.smina_tools_dir,'-r',self.protein_file,'-l',self.ligand_file,'--autobox_ligand',self.native_ligand_file,'--autobox_add','8','--exhaustiveness','32','--num_modes',str(no_modes),'-o',self.pdbqt_output_file,'--atom_terms',self.atom_terms_output_file,'--log',self.log_output_file,'--atom_term_data','--cpu','3','--min_rmsd_filter','0','--energy_range','10000']
        docking = subprocess.run(smina_command, shell=False, capture_output=True, text=True)
    def read_scoring_function(self):

        pdbqt_output = open(self.pdbqt_output_file,'r')
        self.predicted_binding_affinity = [line.replace('REMARK minimizedAffinity ', '').replace('\n', '') for line in pdbqt_output.readlines() if 'REMARK minimizedAffinity' in line]
        return self.predicted_binding_affinity
    def read_atom_term_function(self,no_modes):
        atom_term_output = open(self.atom_terms_output_file,'r')

        gauss = [list() for i in range(no_modes)]
        gauss2 = [list() for i in range(no_modes)]
        repulsion = [list() for i in range(no_modes)]
        hydrophobic = [list() for i in range(no_modes)]
        non_dir_h_bond = [list() for i in range(no_modes)]

        atom_term_output = atom_term_output.read().split('END\n')[:-1]
        atom_term_output = [i.split('\n') for i in atom_term_output]

        for idx,i in enumerate(atom_term_output):
            i = i[1:-1]

            for jdx,j in enumerate(i):
                j=j.split('> ')[-1].split(' ')
                gauss[idx].append(float(j[0]))
                gauss2[idx].append(float(j[1]))
                repulsion[idx].append(float(j[2]))
                hydrophobic[idx].append(float(j[3]))
                non_dir_h_bond[idx].append(float(j[4]))

            gauss[idx] = float(sum(gauss[idx]))
            gauss2[idx] = float(sum(gauss2[idx]))
            repulsion[idx] = float(sum(repulsion[idx]))
            hydrophobic[idx] = float(sum(hydrophobic[idx]))
            non_dir_h_bond[idx] = float(sum(non_dir_h_bond[idx]))
        self.sf_components = [gauss,gauss2,repulsion,hydrophobic,non_dir_h_bond]
        return self.sf_components
    def read_experimental_affinity(self,df,protein,ligand,affinity_column='value'):
        if protein == ligand:
            self.experimental_affinity = df[df['pdbid']==protein][affinity_column].values[0]
    def create_smina_matrix(self,proteins,ligands,no_modes):

        values = ['Predicted Binding Affinity','Gauss1','Gauss2','Repulsion','Hydrophobic','Non_dir_h_bond']

        no_proteins = len(proteins)
        print(no_proteins)
        no_ligands = len(ligands)
        print(no_ligands)
        no_modes = no_modes
        print(no_modes)
        self.matrix = np.empty((no_proteins,no_ligands,no_modes),dtype=object)
        print('elo')
        return self.matrix
    def fill_smina_matrix(self,pidx,lidx):

        for mode_idx in range(len(self.matrix[0][0])):
            if pidx == lidx:
                mode_values = [float(self.predicted_binding_affinity[mode_idx])]
            else:
                mode_values = [float(self.predicted_binding_affinity[mode_idx])]
            for i in range(5):
                mode_values.append(float(self.sf_components[i][mode_idx]))
            mode_values = np.array(mode_values)
            mode_values = tf.convert_to_tensor(mode_values)
            self.matrix[pidx,lidx,mode_idx] = mode_values
    def save_matrix(self,output):
        output = self.datadir+'/docs/'+output
        np.save(output, self.matrix)
class RxDock:
    def __init__(self,datadir,system_file='/home2/sfglab/mwisniewski/PhD/DockOfAbyss/configs/rxdock_config.prm'):

        self.datadir = datadir
        self.rxdock_dir = datadir + '/docking_scores/rxdock'

        self.protein = None
        self.ligand = None
        self.native_ligand = None

        self.protein_file=''
        self.ligand_file=''
        self.native_ligand_file=''

        self.dock_prm_file = '/home2/sfglab/mwisniewski/anaconda3/envs/dock_of_abyss/share/rxdock-2013.1.1_148c5bd1-1/data/scripts/dock.prm'
        self.system_file = system_file          #rx_dock raw system filepath
        self.system_prepared_file = None        #rx_dock prepared system filepath
        self.rx_output = None                   #rx_dock output filepath

        self.matrix = None
        self.mode_values = None
        self.values = None
        self.experimental_affinity = None
    def rxdock_dirs(self):
        if not os.path.exists(self.rxdock_dir):
            makedir = subprocess.run(['mkdir ' + self.rxdock_dir], shell=True, capture_output=True, text=True)
    def rxdock_files(self,protein,ligand,native_ligand):

        self.protein = protein
        self.ligand = ligand
        self.native_ligand = native_ligand

        self.protein_file = self.datadir + '/protein/mol2/' + protein + '_protein.mol2'
        self.ligand_file = self.datadir + '/ligand/sdf/' + ligand + '_ligand.sdf'
        self.native_ligand_file = self.datadir + '/native_ligand/sdf/' + native_ligand + '_ligand.sdf'
        self.system_prepared_file = self.datadir+'/docs/temp/'+self.protein+'-'+self.ligand+'.prm'
        self.rx_output = self.rxdock_dir + '/' + self.protein + '-' + self.ligand
    def rxdock_system_preparation(self):
        with open(self.system_file,'r') as file:
            system_filedata = file.read()

        system_filedata = system_filedata.replace('{title}', self.protein+'-'+self.ligand)                          #define title
        system_filedata = system_filedata.replace('{receptor_file}', self.protein_file)                   #define receptor.mol2 filepath
        system_filedata = system_filedata.replace('{native_ligand_file}', self.native_ligand_file)

        with open(self.system_prepared_file, 'w') as file:
            file.write(system_filedata)
    def files_checker(self):
        check_file=self.datadir+'/docs/temp/'+self.protein+'-'+self.ligand+'_cav1.grd'
        if not os.path.exists(check_file):
            return True
        else:
            return False

    def rxdock_docking(self,no_modes):
        os.environ['RBT_HOME'] = self.datadir
        command = ['rbcavity','-was','-d','-r',self.system_prepared_file]
        result = subprocess.run(command, shell=False, capture_output=True, text=True)
        command = ['rbdock','-i',self.ligand_file,'-o',self.rx_output,'-r',self.system_prepared_file,'-p',self.dock_prm_file,'-n',str(no_modes)]
        result = subprocess.run(command, shell=False, capture_output=True, text=True)
    def create_rxdock_matrix(self,proteins,ligands,no_modes):

        self.values = ['<SCORE>', '<SCORE.INTER>>', '<SCORE.INTER.CONST>', '<SCORE.INTER.POLAR>',
                  '<SCORE.INTER.REPUL>', '<SCORE.INTER.ROT>', '<SCORE.INTER.VDW>', '<SCORE.INTER.norm>', '<SCORE.INTRA>',
                  '<SCORE.INTRA.DIHEDRAL>', '<SCORE.INTRA.DIHEDRAL.0>', '<SCORE.INTRA.POLAR>', '<SCORE.INTRA.POLAR.0>',
                  '<SCORE.INTRA.REPUL>', '<SCORE.INTRA.REPUL.0>', '<SCORE.INTRA.VDW>', '<SCORE.INTRA.VDW.0>',
                  '<SCORE.INTRA.norm>', '<SCORE.RESTR>', '<SCORE.RESTR.CAVITY>', '<SCORE.RESTR.norm>', '<SCORE.SYSTEM>',
                  '<SCORE.SYSTEM.CONST>', '<SCORE.SYSTEM.DIHEDRAL>', '<SCORE.SYSTEM.norm>', '<SCORE.HEAVY>', '<SCORE.norm>']

        no_proteins = len(proteins)
        no_ligands = len(ligands)
        no_modes = no_modes

        self.matrix = np.empty((no_proteins,no_ligands,no_modes),dtype=object)
        return self.matrix
    def fill_rxdock_matrix(self,pidx,lidx,mode_index):
        self.matrix[pidx,lidx,mode_index] = self.mode_values
    def read_experimental_affinity(self,df,protein,ligand,affinity_column='pKa'):
        if protein == ligand:
            self.experimental_affinity = df[df['pdbid']==protein][affinity_column].values[0]
    def read_output(self, protein_index, ligand_index):

        self.rx_output = self.rx_output + '.sd'
        output_file = open(self.rx_output,'r')
        list_of_modes = output_file.read().split('$$$$')[:-1]

        for mode_index, mode in enumerate(list_of_modes[:]):
            self.mode_values = mode.split('>  <SCORE>')[-1].split('\n')[1:-2:3]
            self.mode_values = np.array(self.mode_values)
            self.mode_values = tf.convert_to_tensor(self.mode_values)
            self.fill_rxdock_matrix(protein_index,ligand_index,mode_index)
    def save_matrix(self,output):
        output = self.datadir+'/docs/'+output
        np.save(output, self.matrix)
class Gnina:
    def __init__(self,datadir):
        self.datadir = datadir
        self.gnina_dir = datadir + '/docking_scores/gnina'
        self.sdf_gz_gnina_dir = self.gnina_dir + '/sdf_gz'
        self.atom_terms_gnina_dir = self.gnina_dir + '/atom_terms'
        self.logs_gnina_dir = self.gnina_dir+ '/logs'
        self.protein_file=''
        self.ligand_file=''
        self.native_ligand_file=''
        self.sdf_gz_output_file = ''
        self.atom_terms_output_file =''
        self.log_output_file = ''
        self.matrix = None
        self.experimental_affinity = None
        self.predicted_binding_affinity = None
        self.sf_components=None
    def gnina_dirs(self):
        if not os.path.exists(self.gnina_dir):
            makedir = subprocess.run(['mkdir ' + self.gnina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.sdf_gz_gnina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.atom_terms_gnina_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.logs_gnina_dir], shell=True, capture_output=True, text=True)
    def gnina_files(self,protein,ligand,native_ligand):

        self.protein_file = self.datadir+'/protein/pdb/'+protein+'_protein.pdb'
        self.ligand_file = self.datadir+'/ligand/sdf/'+ligand+'_ligand.sdf'
        self.native_ligand_file = self.datadir+'/native_ligand/sdf/'+native_ligand+'_ligand.sdf'

        self.sdf_gz_output_file = self.sdf_gz_gnina_dir + '/' + protein + '_' + ligand + '.sdf.gz'
        self.atom_terms_output_file = self.atom_terms_gnina_dir + '/' + protein + '_' + ligand + '_atom_terms.txt'
        self.log_output_file = self.logs_gnina_dir + '/' + protein + '_' + ligand + '.log'
    def gnina_files_checker(self):
        checker = not os.path.exists(self.log_output_file)
        return checker
    def gnina_docking(self,no_modes):
        gnina_command=['singularity','run','--nv',settings.gnina_container,'gnina','-r',self.protein_file,'-l',self.ligand_file,'--autobox_ligand',self.native_ligand_file,'--autobox_add','8','--exhaustiveness','32','--num_modes',str(no_modes),'-o',self.sdf_gz_output_file,'--atom_terms',self.atom_terms_output_file,'--log',self.log_output_file,'--atom_term_data','--cpu','3','--min_rmsd_filter','0']
        print('running gnina')
        docking = subprocess.run(gnina_command, shell=False, capture_output=True, text=True)
        print(docking.stdout)
        print(docking.stderr)
