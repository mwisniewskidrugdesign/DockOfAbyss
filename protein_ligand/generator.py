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
        ligand_dirlist = ['sdf','mol2','pdb', 'pdbqt','smi']
        
        for dir in dirlist:
            if not os.path.exists(datadir + '/' + dir):
                result = subprocess.run(['mkdir',datadir+'/'+dir], shell=False, capture_output=True,text=True)

        for subdir in protein_dirlist:
            if not os.path.exists(datadir + '/protein/' + subdir):
                result = subprocess.run(['mkdir',datadir+'/protein/'+subdir], shell=False,capture_output=True, text=True)

        for subdir in pocket_dirlist:
            if not os.path.exists(datadir + '/pocket/' + subdir):
                result = subprocess.run(['mkdir',datadir+'/pocket/'+subdir], shell=False,capture_output=True, text=True)

        for subdir in native_ligand_dirlist:
            if not os.path.exists(datadir + '/native_ligand/' + subdir):
                result = subprocess.run(['mkdir',datadir+'/native_ligand/'+subdir], shell=False,capture_output=True, text=True)

        for subdir in ligand_dirlist:
            if not os.path.exists(datadir + '/ligand/' + subdir):
                result = subprocess.run(['mkdir',datadir + '/ligand/' + subdir], shell=False,capture_output=True, text=True)
        
        if not os.path.exists(datadir+'/doc/temp'):
            result = subprocess.run(['mkdir',datadir+'/docs/temp'],shell=False,capture_output=True,text=True)

class GetDataset:
    def __init__(self,datadir,rawdir,df,pdb_id_column):
        self.datadir = datadir              #directory of data
        self.rawdir = rawdir                #directory of rawdata
        self.df = df                        #dataframe of dataset
        self.pdb_id_column=pdb_id_column    #pdb_id_columnn = 'pdbid' or diffrent
    def bdb2020plus(self,native_ligand=False):
        '''
        Copy files from BindindDB 2020+ dataset
        bdb20plus_dir = directory of BindindDB 2020+ database
        native_ligand = deafult False, if you want native ligand diffrent than docked ligand - True
        '''
        try:
            for index,row in self.df.iterrows():
                print(str(index) + '.' + row[self.pdb_id_column])
                pdb_protein_final_path = self.datadir + '/protein/pdb'
                pdb_native_ligand_final_path = self.datadir + '/native_ligand/pdb'
                sdf_native_ligand_final_path = self.datadir + '/native_ligand/sdf'
                pdb_ligand_final_path = self.datadir + '/ligand/pdb'
                sdf_ligand_final_path = self.datadir + '/ligand/sdf'

                protein_pdb_file = self.rawdir + '/' + row[self.pdb_id_column] + '/protein.pdb'
                native_ligand_pdb_file = self.rawdir + '/' + row[self.pdb_id_column] + '/ligand.pdb'
                native_ligand_sdf_file = self.rawdir + '/' + row[self.pdb_id_column] + '/ligand.sdf'
                ligand_pdb_file = self.rawdir + '/' + row[self.pdb_id_column] + '/ligand.pdb'
                ligand_sdf_file = self.rawdir + '/' + row[self.pdb_id_column] + '/ligand.sdf'

                final_paths = [pdb_protein_final_path, pdb_native_ligand_final_path, sdf_native_ligand_final_path]
                files = [protein_pdb_file,native_ligand_pdb_file,native_ligand_sdf_file]
                strings=['_protein.pdb','_ligand.pdb','_ligand.sdf']

                iterations = range(3)

                if native_ligand == False:
                    final_paths = final_paths+[pdb_ligand_final_path, sdf_ligand_final_path]
                    files = files + [ligand_pdb_file,ligand_sdf_file]
                    strings = strings + ['_ligand.pdb','_ligand.sdf']
                    iterations = range(5)

                for i in iterations:
                    print(final_paths[i])
                    if not os.path.exists(final_paths[i]+'/'+row[self.pdb_id_column]+strings[i]):
                        shutil.copy(files[i],final_paths[i]+'/'+row[self.pdb_id_column]+strings[i])
        except IndexError:
            print('Error!')
    def lp_pdbbind(self,native_ligand=False):
        try:
            for index,row in self.df.iterrows():
                print(str(index) + '.' + row[self.pdb_id_column])

                pdb_protein_final_path = self.datadir + '/protein/pdb'
                pdb_pocket_final_path = self.datadir + '/pocket/pdb'
                mol2_native_ligand_final_path = self.datadir + '/native_ligand/mol2'
                sdf_native_ligand_final_path = self.datadir + '/native_ligand/sdf'
                mol2_ligand_final_path = self.datadir + '/ligand/mol2'
                sdf_ligand_final_path = self.datadir + '/ligand/sdf'

                protein_pdb_file = self.rawdir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_protein.pdb'
                pocket_pdb_file = self.rawdir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_pocket.pdb'
                ligand_mol2_file= self.rawdir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.mol2'
                ligand_sdf_file= self.rawdir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.sdf'
                native_ligand_mol2_file= self.rawdir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.mol2'
                native_ligand_sdf_file= self.rawdir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.sdf'

                final_paths=[pdb_protein_final_path,pdb_pocket_final_path,mol2_native_ligand_final_path,sdf_native_ligand_final_path]
                files = [protein_pdb_file,pocket_pdb_file,native_ligand_mol2_file,native_ligand_sdf_file]
                strings=['_protein.pdb','_pocket.pdb','_ligand.mol2','_ligand.sdf']

                iterations = range(4)

                if native_ligand==False:
                    final_paths = final_paths + [mol2_ligand_final_path,sdf_ligand_final_path]
                    files = files +[ligand_mol2_file,ligand_sdf_file]
                    strings = strings + ['_ligand.mol2','_ligand.sdf']
                    iterations = range(6)

                for i in iterations:
                    if not os.path.exists(final_paths[i]+'/'+row[self.pdb_id_column]+strings[i]):
                        shutil.copy(files[i],final_paths[i]+'/'+row[self.pdb_id_column]+strings[i])
        except IndexError:
            print('Error!')
class Converter:
    def __init__(self,datadir,molecule):
        self.molecule = molecule
        self.datadir = datadir
        settings.init()
    def pdb_to_pdbqt(self,protein=False,pocket=False,ligand=False,native_ligand=False):
        '''
        U will need MGL Tools Package and Obabel
        '''
        if protein==False and pocket==False and ligand == False and native_ligand == False:
            print('There is nothing to do. Precise which molecule type u want to convert (protein/pocket/ligand/native ligand)')
            return None

        if protein == True:

            print('<PROTEIN> MGL PDB TO PDBQT:')

            pdb_protein_file = self.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            pdbqt_protein_file = self.datadir + '/protein/pdbqt/' + self.molecule + '_protein.pdbqt'

            if (not os.path.exists(pdbqt_protein_file) or os.path.getsize(pdbqt_protein_file) == 0) and os.path.exists(pdb_protein_file):
                command = [settings.mgltools_dir+'/bin/pythonsh',settings.mgltools_dir + '/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py','-r',pdb_protein_file,'-o',pdbqt_protein_file,'-A','checkhydrogens']
                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdbqt_result.stderr))+' | '+pdbqt_result.stderr)

                if int(len(pdbqt_result.stderr)) > 0:

                    print('<PROTEIN> OPENBABEL PDB TO PDBQT:')
                    command = [settings.obabel_path,pdb_protein_file,'-O',pdbqt_protein_file]
                    pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                    print('error: ' + str(len(pdbqt_result.stderr)) + ' | ' + pdbqt_result.stderr)

            else:

                pdbqt_convert_error = self.molecule+', There is no PDB Protein file to convert or file already exists.'
                print(pdbqt_convert_error)

        if pocket == True:

            print('<POCKET> MGL PDB TO PDBQT:')

            pdb_pocket_file = self.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'
            pdbqt_pocket_file = self.datadir + '/pocket/pdbqt/' + self.molecule + '_pocket.pdbqt'

            if (not os.path.exists(pdbqt_pocket_file) or os.path.getsize(pdbqt_pocket_file) == 0) and os.path.exists(pdb_pocket_file):
                command = [settings.mgltools_dir + '/bin/pythonsh',settings.mgltools_dir + '/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py','-r', pdb_pocket_file, '-o', pdbqt_pocket_file, '-A', 'checkhydrogens']
                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdbqt_result.stderr))+' | '+pdbqt_result.stderr)

                if int(len(pdbqt_result.stderr)) > 0:

                        print('<POCKET> OPENBABEL PDB TO PDBQT:')

                        command = [settings.obabel_path,pdb_pocket_file,'-O',pdbqt_pocket_file]
                        pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)

                        print('error: ' + str(len(pdbqt_result.stderr)) + ' | ' + pdbqt_result.stderr)
            else:

                pdbqt_convert_error = self.molecule+', There is no PDB Pocket file to convert or file already exists..'
                print(pdbqt_convert_error)

        if ligand == True:

            print('<LIGAND> OPENBABEL PDB TO PDBQT: ')

            pdb_ligand_file = self.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'
            pdbqt_ligand_file = self.datadir + '/ligand/pdbqt/' + self.molecule + '_ligand.pdbqt'

            if (not os.path.exists(pdbqt_ligand_file) or os.path.getsize(pdbqt_ligand_file) == 0) and os.path.exists(pdb_ligand_file):

                command = [settings.obabel_path,pdb_ligand_file,'-O',pdbqt_ligand_file]
                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdbqt_result.stderr))+' | '+pdbqt_result.stderr)

            else:

                pdbqt_convert_error = self.molecule+', There is no PDB Ligand file to convert or file already exists..'
                print(pdbqt_convert_error)

        if native_ligand == True:

            print('<NATIVE LIGAND> PDB TO PDBQT:')

            pdb_native_ligand_file = self.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'
            pdbqt_native_ligand_file = self.datadir + '/native_ligand/pdbqt/' + self.molecule + '_ligand.pdbqt'

            if (not os.path.exists(pdbqt_native_ligand_file) or os.path.getsize(pdbqt_native_ligand_file) == 0) and os.path.exists(pdb_native_ligand_file):

                command = [settings.obabel_path, pdb_native_ligand_file, '-O', pdbqt_native_ligand_file]
                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdbqt_result.stderr))+' | '+pdbqt_result.stderr)

            else:
                pdbqt_convert_error = self.molecule+', There is no PDB Native Ligand file to convert or file already exists..'
                print(pdbqt_convert_error)

    def pdb_to_seq(self,protein=True,pocket=False):

        if protein == False and pocket == False:
            print('There is no molecule to convert. PLease specify what type of molecules you want to convert (protein/pocket).')
            return None

        if protein==True:
            protein_pdb_file = self.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            protein_fasta_file = self.datadir + '/protein/fasta/' + self.molecule + '_protein.fasta'
            if os.path.exists(protein_pdb_file) and not os.path.exists(protein_fasta_file):
                result = subprocess.run([settings.pdb2fasta_path+' '+protein_pdb_file+' > '+protein_fasta_file], shell=False, capture_output=True, text=True, close_fds=True)
                print(result.stderr)

        if pocket==True:
            pocket_pdb_file = self.datadir+'/pocket/pdb/'+self.molecule+'_pocket.pdb'
            protein_fasta_file = self.datadir+'/pocket/fasta/'+self.molecule+'_pocket.fasta'
            if os.path.exists(pocket_pdb_file) and not os.path.exists(pocket_fasta_file):
                result = subprocess.run([settings.pdb2fasta_path+' ' + pocket_pdb_file + ' > ' + pocket_fasta_file], shell=False,capture_output=True, text=True, close_fds=True)
                print(result.stderr)
    def pdb_to_mol(self,protein=False,pocket=False,ligand=False,native_ligand=False):

        '''Convert pdb files to mol2 format'''

        if protein == True:

            pdb_protein_file = self.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            mol2_protein_file = self.datadir + '/protein/mol2/' + self.molecule + '_protein.mol2'

            if (not os.path.exists(mol2_protein_file) or os.path.getsize(mol2_protein_file) == 0) and os.path.exists(pdb_protein_file):

                command = [settings.obabel_path,pdb_protein_file,'-O',mol2_protein_file]
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stderr)

            else:

                mol2_convert_error = self.molecule+', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)

        if pocket == True:

            pdb_pocket_file = self.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'
            mol2_pocket_file = self.datadir + '/pocket/mol2/' + self.molecule + '_pocket.mol2'

            if (not os.path.exists(mol2_pocket_file) or os.path.getsize(mol2_pocket_file) == 0) and os.path.exists(pdb_pocket_file):

                command = [settings.obabel_path,pdb_pocket_file,'-O',mol2_pocket_file]
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stderr)

            else:

                mol2_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)

        if ligand == True:

            pdb_ligand_file = self.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'
            mol2_ligand_file = self.datadir + '/ligand/mol2/' + self.molecule + '_ligand.mol2'

            if (not os.path.exists(mol2_ligand_file) or os.path.getsize(mol2_ligand_file) == 0) and os.path.exists(pdb_ligand_file):
                command = [settings.obabel_path,pdb_ligand_file,'-O',mol2_ligand_file]
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stderr)

            else:

                mol2_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)

        if native_ligand == True:

            pdb_native_ligand_file = self.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'
            mol2_native_ligand_file = self.datadir + '/native_ligand/mol2/' + self.molecule + '_ligand.mol2'

            if (not os.path.exists(mol2_native_ligand_file) or os.path.getsize(mol2_native_ligand_file) == 0) and os.path.exists(pdb_native_ligand_file):
                command = [settings.obabel_path,pdb_native_ligand_file,'-O',mol2_native_ligand_file]
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stdout)

            else:

                mol2_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)
    def mol_to_pdb(self,protein=False,pocket=False,ligand=False,native_ligand=False):
        '''To be done if needed'''

        '''Convert mol2 files to pdb format'''

        if protein == True:

            mol2_protein_file = self.datadir + '/protein/mol2/' + self.molecule + '_protein.mol2'
            pdb_protein_file = self.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'

            if (not os.path.exists(pdb_protein_file) or os.path.getsize(pdb_protein_file) == 0) and os.path.exists(mol2_protein_file):

                command = [settings.obabel_path,mol2_protein_file,'-O',pdb_protein_file]
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule+', There is no MOL2 file to convert or file already exists..'
                print(pdb_convert_error)

        if pocket == True:

            mol2_pocket_file = self.datadir + '/pocket/mol2/' + self.molecule + '_pocket.mol2'
            pdb_pocket_file = self.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'

            if (not os.path.exists(pdb_pocket_file) or os.path.getsize(pdb_pocket_file) == 0) and os.path.exists(mol2_pocket_file):

                command = [settings.obabel_path,mol2_pocket_file,'-O',pdb_pocket_file]
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(pdb_convert_error)

        if ligand == True:

            mol2_ligand_file = self.datadir + '/ligand/mol2/' + self.molecule + '_ligand.mol2'
            pdb_ligand_file = self.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'

            if (not os.path.exists(pdb_ligand_file) or os.path.getsize(pdb_ligand_file) == 0) and os.path.exists(mol2_ligand_file):
                command = [settings.obabel_path,mol2_ligand_file,'-O',pdb_ligand_file]
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(pdb_convert_error)

        if native_ligand == True:

            mol2_native_ligand_file = self.datadir + '/native_ligand/mol2/' + self.molecule + '_ligand.mol2'
            pdb_native_ligand_file = self.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'

            if (not os.path.exists(pdb_native_ligand_file) or os.path.getsize(pdb_native_ligand_file) == 0) and os.path.exists(mol2_native_ligand_file):
                command = [settings.obabel_path,mol2_native_ligand_file,'-O',pdb_native_ligand_file]
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(pdb_convert_error)

    def mol_to_sdf(self,protein=False,pocket=False,ligand=False,native_ligand=False):
        '''To be done if needed'''

        '''Convert sdf files to mol format'''

        if protein == True:

            mol2_protein_file = self.datadir + '/protein/mol2/' + self.molecule + '_protein.mol2'
            sdf_protein_file = self.datadir + '/protein/sdf/' + self.molecule + '_protein.sdf'

            if os.path.exists(mol2_protein_file):

                command = [settings.obabel_path,sdf_protein_file,'-O',mol2_protein_file]
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule+', There is no MOL2 file to convert or file already exists..'
                print(convert_error)

        if pocket == True:

            mol2_pocket_file = self.datadir + '/pocket/mol2/' + self.molecule + '_pocket.mol2'
            sdf_pocket_file = self.datadir + '/pocket/sdf/' + self.molecule + '_pocket.sdf'

            if os.path.exists(mol2_pocket_file):

                command = [settings.obabel_path,mol2_pocket_file,'-O',sdf_pocket_file]
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(convert_error)

        if ligand == True:

            mol2_ligand_file = self.datadir + '/ligand/mol2/' + self.molecule + '_ligand.mol2'
            sdf_ligand_file = self.datadir + '/ligand/sdf/' + self.molecule + '_ligand.sdf'

            if os.path.exists(mol2_ligand_file):
                command = [settings.obabel_path,mol2_ligand_file,'-O',sdf_ligand_file]
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(convert_error)

        if native_ligand == True:

            mol2_native_ligand_file = self.datadir + '/native_ligand/mol2/' + self.molecule + '_ligand.mol2'
            sdf_native_ligand_file = self.datadir + '/native_ligand/sdf/' + self.molecule + '_ligand.sdf'

            if os.path.exists(mol2_native_ligand_file):
                command = [settings.obabel_path,mol2_native_ligand_file,'-O',sdf_native_ligand_file]
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(convert_error)
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
        cluster_program=subprocess.run([command],shell=False, capture_output=True, text=True)
        output_file = output_file + '.clstr'
        return output_file