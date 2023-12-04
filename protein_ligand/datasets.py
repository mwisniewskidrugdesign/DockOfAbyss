import pandas as pd
from Bio import SeqIO

import settings


def get_list_of_molecules(df,pdb_id_column='pdbid'):

    molecules = df[pdb_id_column].values.tolist()

    return molecules


class DatasetPreparation:
    '''Class to preapre dataset'''
    def __init__(self,df):
        self.df = df
    def create_affinity_columns(self, binding_affinity_column='Affinity Data'):
        '''Create Kd, Ki and IC50 columns from PDBbind.csv Binding Affinity column'''
        print('Creating Affinity Column')
        if 'Kd [nM]' in self.df.columns and 'Ki [nM]' in self.df.columns and 'IC50 [nM]' in self.df.columns:
            print('kurwaaaaa')
        else:
            self.df['Kd [nM]'] = ''
            self.df['Ki [nM]'] = ''
            self.df['IC50 [nM]'] = ''
            self.df['temp_affinity_value'] = ''
            self.df['temp_affinity_unit'] = ''
            self.df[binding_affinity_column] = self.df[binding_affinity_column].astype('str')
            for index, row in self.df.iterrows():
                ''' Change signs to equal sign'''
                if '<=' in row[binding_affinity_column]:
                    self.df.at[index, binding_affinity_column] = row[binding_affinity_column].replace('<=', '=')
                elif '>=' in row[binding_affinity_column]:
                    self.df.at[index, binding_affinity_column] = row[binding_affinity_column].replace('>=', '=')
                elif '<' in row[binding_affinity_column]:
                    self.df.at[index, binding_affinity_column] = row[binding_affinity_column].replace('<', '=')
                elif '>' in row[binding_affinity_column]:
                    self.df.at[index, binding_affinity_column] = row[binding_affinity_column].replace('>', '=')
                elif '~' in row[binding_affinity_column]:
                    self.df.at[index, binding_affinity_column] = row[binding_affinity_column].replace('~', '=')

            self.df[[binding_affinity_column, 'temp_affinity_value']] = self.df[binding_affinity_column].str.split('=', n=1, expand=True)  # split Affinity data into Affinity data Measure and temporary affinity value with unit
            self.df['temp_affinity_unit'] = self.df['temp_affinity_value'].str[-2:]  # column for Affinity Unit: nM, pM, fM, ...
            self.df['temp_affinity_value'] = self.df['temp_affinity_value'].str[:-2].astype('float64')  # column for Affinity value
            self.df['temp_affinity_type'] = self.df[binding_affinity_column]  # column for Affinity type: Kd, Ki or IC50

            for index, row in self.df.iterrows():
                ''' Unify values to [nM]'''
                if 'fM' in row['temp_affinity_unit']:
                    self.df.at[index, 'temp_affinity_value'] = row['temp_affinity_value'] * 1000000
                if 'pM' in row['temp_affinity_unit']:
                    self.df.at[index, 'temp_affinity_value'] = row['temp_affinity_value'] * 1000
                if 'uM' in row['temp_affinity_unit']:
                    self.df.at[index, 'temp_affinity_value'] = row['temp_affinity_value'] * 0.001
                if 'mM' in row['temp_affinity_unit']:
                    self.df.at[index, 'temp_affinity_value'] = row['temp_affinity_value'] / 1000000

            for index, row in self.df.iterrows():
                ''' Transfer to proper column'''
                if 'Kd' in row['temp_affinity_type']:
                    self.df.at[index, 'Kd [nM]'] = row['temp_affinity_value']
                if 'Ki' in row['temp_affinity_type']:
                    self.df.at[index, 'Ki [nM]'] = row['temp_affinity_value']
                if 'IC50' in row['temp_affinity_type']:
                    self.df.at[index, 'IC50 [nM]'] = row['temp_affinity_value']

            self.df = self.df.drop(columns=['temp_affinity_value', 'temp_affinity_unit', 'temp_affinity_type',binding_affinity_column])  # del temporary columns

        return self.df
    def create_seq_columns(self, datadir, pdb_id_column='PDB code',seq_column_name='protein_sequence',sequence_type='protein'):

        '''Create protein or pocket aas sequence columns from PDB fasta files'''
        print('Creating Sequence Column')
        for index, row in self.df.iterrows():

            sequence_lines = []

            sequence_fasta_filepath = datadir + '/'+sequence_type+'/fasta/' + row[pdb_id_column] + '_'+sequence_type+'.fasta'

            for seq_record in SeqIO.parse(sequence_fasta_filepath,'fasta'):
                sequence_lines.append(str(seq_record.seq))

            sequence_lines = ''.join(sequence_lines)
            row[seq_column_name] = sequence_lines
            self.df.at[index, seq_column_name] = row[seq_column_name]

        return self.df
    def dataframe_filters(self,nmr=False,max_resolution=None,affinity_types=['Ki [nM]','IC50 [nM]']):
        if nmr:
            self.df = self.df[self.df['Resolution'].astype(str) != 'NMR']
        if max_resolution is not None:
            self.df= self.df[self.df['Resolution'] <= max_resolution]
        for affinity_type in affinity_types:
            if affinity_type in self.df.columns:
                self.df = self.df[self.df[affinity_type].str.len() == 0]
        return self.df

    def drop_duplicate_ligands(self):
        self.df = self.df.drop_duplicates(subset=['Canonical SMILES'])
        return self.df

    def save_df(self, datadir, output):
        self.df = pd.read_excel(self.df)
        self.df.to_excel(datadir + '/protein_ligand/docs/' + output + '.xlsx')
    def get_molecules(self,pdb_id_column):
        molecules = self.df[pdb_id_column].tolist()
        return molecules

class ReadDocuments:
    '''Read Files'''
    def __init__(self,filepath):
        self.filepath = filepath
    def read_cdhit_file(self,df,pdb_id_column='PDB code'):
        with open(self.filepath,'r') as file:
            lines = file.readlines()
        representants = [line.split('>')[-1][:4].strip() for line in lines if '*' in line]
        df = df[df[pdb_id_column].isin(representants)]

        return df, representants


