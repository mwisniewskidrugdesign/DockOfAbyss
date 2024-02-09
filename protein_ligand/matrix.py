import pandas as pd
import numpy as np
from protein_ligand.docking import Smina, Gnina, RxDock
import os
from typing import List
import settings
from utils.read_files import process_smina_atom_terms_block

class ResultsMatrix:
    '''
    The structure of every matrix is following:
    (x,y,z)
        x is protein index
        y is ligand index
        z is model index
    How (x,y,z) looks like:
        [Real affinity, RMSD value, Predicted affinity, [**atom_terms]]
    Additionaly, thera are generated files with protein and ligand indexing
    '''
    def __init__(self, proteins: List, ligands: List, no_components: int, df: pd.DataFrame):
        settings.init()
        self.proteins = proteins
        self.ligands = ligands
        self.no_components = no_components
        self.df = pd.read_csv(df)
        self.matrix = None
    def create_matrix(self):
        ''' This method can create an empty matrix for all type of programs.'''
        no_proteins = len(self.proteins)
        no_ligands = len(self.ligands)
        self.matrix = np.empty((no_proteins,no_ligands,settings.number_of_models,self.no_components),dtype=object)

        return self.matrix
class SminaResultsMatrix(ResultsMatrix):
    '''
    Create the Smina Result Matrix, where (x,y,z) looks like:
        [
        Real affinity,
        RMSD value,
        ['Gauss1','Gauss2','Repulsion','Hydrophobic','Non_dir_h_bond']
        ]
    '''
    def __init__(self,proteins: List, ligands: List, df: pd.DataFrame):
        super().__init__(proteins,ligands, no_components,df)
        self.protein=None
        self.ligand=None
        self.experimental_affinity = None
        self.predicted_binding_affinity_per_complex = None
        self.scoring_components_per_complex = None
        self.lb_rmsd_per_complex = None
        self.ub_rmsd_per_complex = None

    def read_experimental_affinity(self,affinity_column='value'):
        '''
        This method reads Binding Affinity from the Leak-Proof PDBBind dataframe.
        :param protein: Complex's protein
        :param ligand: Complex's ligand
        :return: Experimental affinity
        '''
        if self.protein == self.ligand:
            self.experimental_affinity = float(self.df[self.df['pdbid']==self.protein][affinity_column].values[0])
        else:
            self.experimental_affinity = np.nan

        return self.experimental_affinity
    def read_RMSD(self):
        '''
        This method reads SMINA RMSD for all models from logfile of one docked complex.
        '''
        log_output_file = settings.datadir+'/docking_scores/smina/logs/'+self.protein+'_'+self.ligand+'.log'
        with open(log_output_file,'r') as log_output:
            lines = [line for line in log_output.readlines()[25:]]
            self.lb_rmsd_per_complex = [line[19:24] for line in lines]
            self.ub_rmsd_per_complex = [line[30:35] for line in lines]

        return self.lb_rmsd_per_complex,self.ub_rmsd_per_complex
    def read_predicted_binding_affinity(self):
        '''
        This method reads SMINA Predicted Binding Affinities for all models from pdbqt file of one docked complex.
        '''
        pdbqt_output_file = settings.datadir+'/docking_scores/smina/pdbqt/'+self.protein+'_'+self.ligand+'.pdbqt'
        with open(pdbqt_output_file,'r') as pdbqt_output:
            self.predicted_binding_affinitiy_per_complex = [float(line.replace('REMARK minimizedAffinity ', '').replace('\n', ''))for line in pdbqt_output.readlines() if 'REMARK minimizedAffinity' in line]
        return self.predicted_binding_affinitiy_per_complex
    def read_atom_terms(self):
        '''
        This method reads Atom Terms for all models from _atom_terms.txt file of one specified complex.
        '''

        atom_terms_output_file = settings.datadir+'/docking_scores/smina/atom_terms/'+self.protein+'_'+self.ligand+'_atom_terms.txt'
        with open(atom_terms_output_file, 'r') as atom_term_output:
            content = atom_term_output.read()
            blocks = content.split('atomid')[1:]
            self.scoring_components_per_complex = []
            for block in blocks:
                atomid, element, coordinates, scoring_components_sum = process_smina_atom_terms_block(block)
                self.scoring_components_per_complex.append(scoring_components_sum)

        return self.scoring_components_per_complex
    def fill_matrix(self,protein_index,ligand_index,model_index):

        self.matrix[protein_index,ligand_index,model_index,0] = self.ub_rmsd_per_complex[model_index]
        self.matrix[protein_index,ligand_index,model_index,1] = self.predicted_binding_affinitiy_per_complex[model_index]
        self.matrix[protein_index,ligand_index,model_index,2:] = np.array(self.scoring_components_per_complex[mode_index])
    def save_matrix(self,output_filename: str):
        np.save(settings.datadir+'/docs/results/'+output_filename+'.npy',self.matrix)
class RxDockResultsMatrix(ResultsMatrix):
    def __init__(self,proteins: List, ligands: List, df: pd.DataFrame):
        super().__init__(proteins,ligands, no_components,df)
        self.protein=None
        self.ligand=None
        self.experimental_affinity = None
        self.predicted_binding_affinity_per_complex = None
        self.scoring_components_per_complex = None
        self.rmsd_per_complex = None
    def read_experimental_affinity(self,affinity_column='value'):
        '''
        This method reads Binding Affinity from the Leak-Proof PDBBind dataframe.
        :param protein: Complex's protein
        :param ligand: Complex's ligand
        :return: Experimental affinity
        '''
        if self.protein == self.ligand:
            self.experimental_affinity = float(self.df[self.df['pdbid']==self.protein][affinity_column].values[0])
        else:
            self.experimental_affinity = np.nan

        return self.experimental_affinity
    def read_data(self):
        sdf_file_path = settings.datadir + '/docking_scores/rxdock/'+self.protein+'_'+self.ligand+'.sdf'
        with open(sdf_file_path, 'r') as output_file:
            list_of_modes = output_file.read().split('$$$$')[:-1]
            modes = []
            for mode_index, mode in enumerate(list_of_modes):
                one_model_values = ['SCORE'] + mode.split('>  <SCORE>')[-1].split('\n')[1:]
                one_model_values = one_model_values[:]
                keys = one_model_values[::3]
                values = one_model_values[1::3]
                keys = [key.replace('>  <', '').replace('>', '') for key in keys]
                together = dict(zip(keys, values))
                modes.append(together)
            unique = set(key for mode in modes for key in mode)
            fill_dicts = {key: 0 for key in unique}
            for mode in modes:
                for key in unique:
                    mode.setdefault(key, 0)
            df = pd.DataFrame(modes)
            self.rmsd_per_complex = df['RMSD'].tolist()
            self.predicted_binding_affinity_per_complex = df['SCORE'].tolist()
            df_components = df.drop(['RMSD', 'SCORE'], axis=1)
            self.scoring_components_per_complex = [row.tolist() for _, row in df_components.iterrows()]
        return self.rmsd_per_complex, self.predicted_binding_affinity_per_complex, self.scoring_components_per_complex
    def fill_matrix(self,protein_index,ligand_index,model_index):

        self.matrix[protein_index,ligand_index,model_index,0] = self.rmsd_per_complex[model_index]
        self.matrix[protein_index,ligand_index,model_index,1] = self.predicted_binding_affinitiy_per_complex[model_index]
        self.matrix[protein_index,ligand_index,model_index,2:] = np.array(self.scoring_components_per_complex[mode_index])
    def save_matrix(self,output_filename: str):
        np.save(settings.datadir+'/docs/results/'+output_filename+'.npy',self.matrix)
class GninaResultsMatrix(ResultsMatrix):
    def __init__(self,proteins: List, ligands: List, df: pd.DataFrame):
        super().__init__(proteins,ligands, no_components,df)
        self.protein=None
        self.ligand=None
        self.experimental_affinity = None
        self.predicted_binding_affinity_per_complex = None
        self.scoring_components_per_complex = None
        self.rmsd_per_complex = None
    def read_rmsd(self):
        return None
    def read_predicted_affinity(self):
        return None
    def read_atom_terms(self):
        return None