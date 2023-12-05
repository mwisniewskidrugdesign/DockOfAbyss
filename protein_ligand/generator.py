import os
import subprocess
import shutil
import subprocess
import logging

import settings

def setup_logger(name,log_file,level=logging.INFO):
    """
    To setup as many loggers as you want
    """

    handler=logging.FileHandler(log_file,'a')
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger
def generate_libraray(datadir):
    '''Create dataset_library'''

    if os.path.exists(datadir):
        dirlist = ['protein', 'native_ligand', 'ligand', 'pocket', 'logs', 'docs', 'docking_scores']
        protein_dirlist = ['mol2','pdb', 'pdbqt', 'fasta']
        pocket_dirlist = ['mol2','pdb', 'pdbqt', 'fasta']
        native_ligand_dirlist = ['sdf', 'mol2', 'pdb', 'pdbqt', 'smi']
        ligand_dirlist = ['mol2','pdb', 'pdbqt']

        for dir in dirlist:
            if not os.path.exists(datadir + '/' + dir):
                result = subprocess.run(['mkdir ' + datadir + '/' + dir], shell=True, capture_output=True,text=True)

        for subdir in protein_dirlist:
            if not os.path.exists(datadir + '/protein/' + subdir):
                result = subprocess.run(['mkdir ' + datadir + '/protein/' + subdir], shell=True,capture_output=True, text=True)

        for subdir in pocket_dirlist:
            if not os.path.exists(datadir + '/pocket/' + subdir):
                result = subprocess.run(['mkdir ' + datadir + '/pocket/' + subdir], shell=True,capture_output=True, text=True)

        for subdir in native_ligand_dirlist:
            if not os.path.exists(datadir + '/native_ligand/' + subdir):
                result = subprocess.run(['mkdir ' + datadir + '/native_ligand/' + subdir], shell=True,capture_output=True, text=True)

        for subdir in ligand_dirlist:
            if not os.path.exists(datadir + '/ligand/' + subdir):
                result = subprocess.run(['mkdir ' + datadir + '/ligand/' + subdir], shell=True,capture_output=True, text=True)
class GetDataset:
    def __init__(self,datadir,df,pdb_id_column='pdbid'):
        self.datadir = datadir              #directory of data
        self.df = df                        #dataframe of dataset
        self.pdb_id_column=pdb_id_column    #pdb_id_columnn = 'pdbid' or diffrent
    def bdb2020plus(self,bdb2020plus_dir,native_ligand=False):
        '''
        Copy files from BindindDB 2020+ dataset
        bdb20plus_dir = directory of BindindDB 2020+ database
        native_ligand = deafult False, if you want native ligand diffrent than docked ligand - True
        '''
        try:
            for index,row in self.df.iterrows():
                print(str(index) + '.' + row[self.pdb_id_column])
                pdb_protein_final_path = self.datadir + '/protein/pdb'
                pdb_ligand_final_path = self.datadir + '/ligand/pdb'
                pdb_native_ligand_final_path = self.datadir + '/native_ligand/pdb'
                sdf_native_ligand_final_path = self.datadir + '/native_ligand/sdf'

                protein_pdb_file = bdb2020plus_dir + '/' + row[self.pdb_id_column] + '/protein.pdb'
                ligand_pdb_file = bdb2020plus_dir + '/' + row[self.pdb_id_column] + '/ligand.pdb'
                native_ligand_pdb_file = bdb2020plus_dir + '/' + row[self.pdb_id_column] + '/ligand.pdb'
                native_ligand_sdf_file = bdb2020plus_dir + '/' + row[self.pdb_id_column] + '/ligand.sdf'

                if not os.path.exists(pdb_protein_final_path + '/' + row[self.pdb_id_column] + '_protein.pdb'):
                    shutil.copy(protein_pdb_file, pdb_protein_final_path + '/' + row[self.pdb_id_column] + '_protein.pdb')

                if not os.path.exists(pdb_native_ligand_final_path + '/' + row[self.pdb_id_column] + '_ligand.pdb'):
                    shutil.copy(native_ligand_pdb_file, pdb_native_ligand_final_path + '/' + row[self.pdb_id_column] + '_ligand.pdb')

                if not os.path.exists(sdf_native_ligand_final_path + '/' + row[self.pdb_id_column] + '_ligand.sdf'):
                    shutil.copy(native_ligand_sdf_file, sdf_native_ligand_final_path + '/' + row[self.pdb_id_column] + '_ligand.sdf')

                if native_ligand == False:
                    if not os.path.exists(pdb_ligand_final_path + '/' + row[self.pdb_id_column] + '_ligand.pdb'):
                        shutil.copy(ligand_pdb_file, pdb_ligand_final_path + '/' + row[self.pdb_id_column] + '_ligand.pdb')

        except IndexError:
            print('Error!')
class Converter:
    def __init__(self,molecule,datadir):
        self.molecule = molecule
        self.datadir = datadir
    def pdb_to_pdbqt(self,protein=False,pocket=False,ligand=False,native_ligand=False):
        '''
        U will need MGL Tools Package
        '''
        if protein==False and pocket==False and ligand == False and native_ligand == False:
            print('There is nothing to do. Precise which molecule type u want to convert (protein/pocket/ligand/native ligand)')
            break
        if protein == True:
            '''Generate Log File'''
            generate_protein_pdbqt_log_file = self.datadir + '/logs/generate_protein_pdbqt_files.log'
            generate_protein_pdbqt_file_logger = setup_logger('generate_protein_pdbqt_file_logger',generate_protein_pdbqt_log_file)

            pdb_protein_file = self.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            pdbqt_protein_file = self.datadir + '/protein/pdbqt/' + self.molecule + '_protein.pdbqt'

            if (not os.path.exists(pdbqt_protein_file) or os.path.getsize(pdbqt_protein_file) == 0) and os.path.exists(pdb_protein_file):
                command = settings.mgltools_dir + '/bin/pythonsh ' + settings.mgltools_dir + '/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ' + pdb_protein_file + ' -o ' + pdbqt_protein_file + ' -A checkhydrogens'
                pdbqt_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(pdbqt_result.stderr)
                generate_protein_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_result.stderr)
            else:
                pdbqt_convert_error = self.molecule+', There is no PDB file to convert.'
                generate_protein_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_convert_error)
        if pocket == True:
            '''Generate Log File'''
            generate_pocket_pdbqt_log_file = self.datadir + '/logs/generate_pocket_pdbqt_files.log'
            generate_pocket_pdbqt_file_logger = setup_logger('generate_pocket_pdbqt_file_logger',generate_pocket_pdbqt_log_file)

            pdb_pocket_file = self.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'
            pdbqt_pocket_file = self.datadir + '/pocket/pdbqt/' + self.molecule + '_pocket.pdbqt'

            if (not os.path.exists(pdbqt_pocket_file) or os.path.getsize(pdbqt_pocket_file) == 0) and os.path.exists(pdb_pocket_file):
                command = settings.mgltools_dir + '/bin/pythonsh ' + settings.mgltools_dir + '/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ' + pdb_pocket_file + ' -o ' + pdbqt_pocket_file + ' -A checkhydrogens'
                pdbqt_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(pdbqt_result.stderr)
                generate_pocket_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_result.stderr)
            else:
                pdbqt_convert_error = self.molecule+', There is no PDB file to convert.'
                generate_protein_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_convert_error)


        if ligand == True:
            '''Generate Log File'''
            generate_ligand_pdbqt_log_file = self.datadir + '/logs/generate_ligand_pdbqt_files.log'
            generate_ligand_pdbqt_file_logger = setup_logger('generate_ligand_pdbqt_file_logger',generate_ligand_pdbqt_log_file)

            pdb_ligand_file = self.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'
            pdbqt_ligand_file = self.datadir + '/ligand/pdbqt/' + self.molecule + '_ligand.pdbqt'

            if (not os.path.exists(pdbqt_ligand_file) or os.path.getsize(pdbqt_ligand_file) == 0) and os.path.exists(pdb_ligand_file):
                command = settings.obabel_dir+' '+pdb_ligand_file+' -O '+pdbqt_ligand_file
                pdbqt_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(pdbqt_result.stderr)
                generate_ligand_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_result.stdout)
            else:
                pdbqt_convert_error = self.molecule+', There is no PDB file to convert.'
                generate_protein_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_convert_error)

        if native_ligand == True:
            '''Generate Log File'''
            generate_native_ligand_pdbqt_log_file = self.datadir + '/logs/generate_native_ligand_pdbqt_files.log'
            generate_native_ligand_pdbqt_file_logger = setup_logger('generate_native_ligand_pdbqt_file_logger',generate_native_ligand_pdbqt_log_file)

            pdb_native_ligand_file = self.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'
            pdbqt_native_ligand_file = self.datadir + '/native_ligand/pdbqt/' + self.molecule + '_ligand.pdbqt'

            if (not os.path.exists(pdbqt_native_ligand_file) or os.path.getsize(pdbqt_native_ligand_file) == 0) and os.path.exists(pdb_native_ligand_file):
                command = settings.obabel_dir +' '+ pdb_native_ligand_file + ' -O ' + pdbqt_native_ligand_file
                pdbqt_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(pdbqt_result.stdout)
                generate_native_ligand_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_result.stdout)
            else:
                pdbqt_convert_error = self.molecule+', There is no PDB file to convert.'
                generate_protein_pdbqt_file_logger.info(self.molecule + ',' + pdbqt_convert_error)
    def pdb_to_seq(self,protein=True,pocket=False):

        if protein == False and pocket == False:
            print('There is no molecule to convert. PLease specify what type of molecules you want to convert (protein/pocket).')
            break
        if protein==True:
            protein_pdb_file = self.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            protein_fasta_file = self.datadir + '/protein/fasta/' + self.molecule + '_protein.fasta'
            if os.path.exists(protein_pdb_file) and not os.path.exists(protein_fasta_file):
                result = subprocess.run([settings.pdb2fasta_path+' '+protein_pdb_file+' > '+protein_fasta_file], shell=True, capture_output=True, text=True)
                print(result.stderr)

        if pocket==True:
            pocket_pdb_file = self.datadir+'/pocket/pdb/'+self.molecule+'_pocket.pdb'
            protein_fasta_file = self.datadir+'/pocket/fasta/'+self.molecule+'_pocket.fasta'
            if os.path.exists(pocket_pdb_file) and not os.path.exists(pocket_fasta_file):
                result = subprocess.run([settings.pdb2fasta_path+' ' + pocket_pdb_file + ' > ' + pocket_fasta_file], shell=True,capture_output=True, text=True)
                print(result.stderr)
    def pdb_to_mol(self,protein=False,pocket=False,ligand=False,native_ligand=False):
        '''Convert pdb files to mol2 format'''
        if protein == True:
            '''Generate Log File'''
            generate_protein_mol2_log_file = self.datadir + '/logs/generate_protein_mol2_files.log'
            generate_protein_mol2_file_logger = setup_logger('generate_protein_mol2_file_logger',generate_protein_mol2_log_file)

            pdb_protein_file = self.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            mol2_protein_file = self.datadir + '/protein/mol2/' + self.molecule + '_protein.mol2'

            if (not os.path.exists(mol2_protein_file) or os.path.getsize(mol2_protein_file) == 0) and os.path.exists(pdb_protein_file):
                command = settings.obabel_dir +' '+ pdb_protein_file + ' -O ' + mol2_protein_file
                mol2_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(mol2_result.stderr)
                generate_protein_mol2_file_logger.info(self.molecule + ',' + mol2_result.stderr)
            else:
                mol2_convert_error = self.molecule+', There is no PDB file to convert.'
                generate_protein_mol2_file_logger.info(self.molecule + ',' + mol2_convert_error)

        if pocket == True:
            '''Generate Log File'''
            generate_pocket_mol2_log_file = self.datadir + '/logs/generate_pocket_mol2_files.log'
            generate_pocket_mol2_file_logger = setup_logger('generate_pocket_mol2_file_logger',generate_pocket_mol2_log_file)

            pdb_pocket_file = self.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'
            mol2_pocket_file = self.datadir + '/pocket/mol2/' + self.molecule + '_pocket.mol2'

            if (not os.path.exists(mol2_pocket_file) or os.path.getsize(mol2_pocket_file) == 0) and os.path.exists(pdb_pocket_file):
                command = settings.obabel_dir +' '+ pdb_pocket_file + ' -O ' + mol2_pocket_file
                mol2_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(mol2_result.stderr)
                generate_pocket_mol2_file_logger.info(self.molecule + ',' + mol2_result.stderr)
            else:
                mol2_convert_error = self.molecule + ', There is no PDB file to convert.'
                generate_protein_mol2_file_logger.info(self.molecule + ',' + mol2_convert_error)


        if ligand == True:
            '''Generate Log File'''
            generate_ligand_mol2_log_file = self.datadir + '/logs/generate_ligand_mol2_files.log'
            generate_ligand_mol2_file_logger = setup_logger('generate_ligand_mol2_file_logger',generate_ligand_mol2_log_file)

            pdb_ligand_file = self.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'
            mol2_ligand_file = self.datadir + '/ligand/mol2/' + self.molecule + '_ligand.mol2'

            if (not os.path.exists(mol2_ligand_file) or os.path.getsize(mol2_ligand_file) == 0) and os.path.exists(pdb_ligand_file):
                command = settings.obabel_dir+' '+pdb_ligand_file+' -O '+mol2_ligand_file
                mol2_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(mol2_result.stderr)
                generate_ligand_mol2_file_logger.info(self.molecule + ',' + mol2_result.stdout)
            else:
                mol2_convert_error = self.molecule + ', There is no PDB file to convert.'
                generate_protein_mol2_file_logger.info(self.molecule + ',' + mol2_convert_error)

        if native_ligand == True:
            '''Generate Log File'''
            generate_native_ligand_mol2_log_file = self.datadir + '/logs/generate_native_ligand_mol2_files.log'
            generate_native_ligand_mol2_file_logger = setup_logger('generate_native_ligand_mol2_file_logger',generate_native_ligand_mol2_log_file)

            pdb_native_ligand_file = self.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'
            mol2_native_ligand_file = self.datadir + '/native_ligand/mol2/' + self.molecule + '_ligand.mol2'

            if (not os.path.exists(mol2_native_ligand_file) or os.path.getsize(mol2_native_ligand_file) == 0) and os.path.exists(pdb_native_ligand_file):
                command = settings.obabel_dir +' '+ pdb_native_ligand_file + ' -O ' + mol2_native_ligand_file
                mol2_result = subprocess.run([command], shell=True, capture_output=True, text=True)
                print(mol2_result.stdout)
                generate_native_ligand_mol2_file_logger.info(self.molecule + ',' + mol2_result.stdout)
            else:
                mol2_convert_error = self.molecule + ', There is no PDB file to convert.'
                generate_protein_mol2_file_logger.info(self.molecule + ',' + mol2_convert_error)
    def sdf_to_mol(self):
        '''To be done if needed'''
        return None
class Documents:
    def __init__(self,datadir):
        self.datadir = datadir
    def merged_fasta(self,df,pdb_id_column,output,protein=True,pocket=False):

        protein_merged_seq_output_file = self.datadir + '/docs/merged_' + output + '_protein.fasta'
        pocket_merged_seq_output_file = self.datadir + '/docs/merged_' + output + '_pocket.fasta'

        protein_records=[]
        pocket_records=[]

        for index, row in df.iterrows():

            pdb_id = row[pdb_id_column]

            if protein == True:
                protein_sequence = row['protein_sequence']
                protein_seq_record = SeqRecord(Seq(protein_sequence),id=pdb_id)
                protein_records.append(protein_seq_record)
            if pocket == True:
                pocket_sequence = row['pocket_sequence']
                pocket_seq_record = SeqRecord(Seq(pocket_sequence),id=pdb_id)
                pocket_records.append(pocket_seq_record)
        if protein == True:
            SeqIO.write(protein_records,protein_merged_seq_output_file,"fasta")
        if pocket == True:
            SeqIO.write(pocket_records,pocket_merged_seq_output_file,"fasta")

        return protein_merged_seq_output_file, pocket_merged_seq_output_file
    def generate_cdhit_cluster_file(self,fasta_file,output,c=0.8,n=5,m=16000,d=0,t=8):
        output_file = self.datadir + '/docs/' + output
        command=settings.cdhit_dir+'/cd-hit -i '+fasta_file+' -o '+output_file+' -c '+str(c)+' -n '+str(n)+' -M '+str(m)+' -d '+str(d)+' -T '+str(t)
        cluster_program=subprocess.run([command],shell=True, capture_output=True, text=True)
        output_file = output_file + '.clstr'
        return output_file