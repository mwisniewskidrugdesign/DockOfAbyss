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
    def smina_docking(self,no_modes):
        smina_command = settings.smina_tools_dir+' -r ' + self.protein_file + ' -l ' + self.ligand_file + ' --autobox_ligand ' + self.native_ligand_file + ' --autobox_add 8 --exhaustiveness 8 --num_modes '+str(no_modes)+' -o ' + self.pdbqt_output_file + ' --atom_terms ' + self.atom_terms_output_file + ' --log ' + self.log_output_file + ' --atom_term_data --cpu 12 --min_rmsd_filter 0 --energy_range 1000000'
        docking = subprocess.run([smina_command], shell=True, capture_output=True, text=True)
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
        #for i in atom_term_output:
            #print(i)
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
    def read_experimental_affinity(self,df,protein,ligand):
        if protein == ligand:
            self.experimental_affinity = df[df['pdbid']==protein]['pKa'].values[0]
    def create_smina_matrix(self,proteins,ligands,no_modes):

        values = ['pKa','Predicted Binding Affinity','Gauss1','Gauss2','Repulsion','Hydrophobic','Non_dir_h_bond']

        no_proteins = len(proteins)
        no_ligands = len(ligands)
        no_modes = no_modes

        self.matrix = np.empty((no_proteins,no_ligands,no_modes),dtype=object)
        return self.matrix
    def fill_smina_matrix(self,pidx,lidx):
        for mode_idx in range(len(self.matrix[0][0])):
            if pidx == lidx:
                mode_values = [(self.experimental_affinity),float(self.predicted_binding_affinity[mode_idx])]
            else:
                mode_values = [np.NaN,float(self.predicted_binding_affinity[mode_idx])]
            for i in range(5):
                mode_values.append(float(self.sf_components[i][mode_idx]))
            mode_values = np.array(mode_values)
            mode_values = tf.convert_to_tensor(mode_values)
            self.matrix[pidx,lidx,mode_idx] = mode_values
    def save_matrix(self,output):
        output = self.datadir+'/docs/'+output
        np.save(output, self.matrix)

class RxDock:
    def __init__(self,datadir):
        self.datadir = datadir
        self.rxdock_dir = datadir + '/docking_scores/rxdock'
        self.pdbqt_rxdock_dir = self.rxdock_dir + '/pdbqt'
        self.atom_terms_rxdock_dir = self.rxdock_dir + '/atom_terms'
        self.logs_rxdock_dir = self.rxdock_dir + '/logs'
    def rxdock_dirs(self):
        if not os.path.exists(self.rxdock_dir):
            makedir = subprocess.run(['mkdir ' + self.rxdock_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.pdbqt_rxdock_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.atom_terms_rxdock_dir], shell=True, capture_output=True, text=True)
            makedir = subprocess.run(['mkdir ' + self.logs_rxdock_dir], shell=True, capture_output=True, text=True)
    def rxdock_files(self,protein,ligand,native_ligand):
        return None
    def rxdock_system_preparation(self):
        return None
    def rxdock_docking(self):
        return None
    def read_files(self):
        return None




