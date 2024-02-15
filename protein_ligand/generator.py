import os
import subprocess
import shutil
import subprocess
import logging
import settings
from protein_ligand.docking import Smina, Gnina, RxDock

def setup_logger(name,log_file,level=logging.INFO):
    """
    To setup as many loggers as you want
    """

    handler=logging.FileHandler(log_file,'a')
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger
def generate_libraray():
    '''Create dataset_library'''
    settings.init()

    datadir = settings.datadir

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

        if not os.path.exists(datadir + '/doc/results'):
            result = subprocess.run(['mkdir', datadir + '/docs/results'], shell=False, capture_output=True, text=True)
class GetDataset:
    def __init__(self,df,pdb_id_column):
        settings.init()
        self.pdb_id_column=pdb_id_column    #  data frame column with specified pdb_id column
    def bdb2020plus(self,native_ligand=False):
        '''
        Copy files from BindindDB 2020+ dataset
        bdb20plus_dir = directory of BindindDB 2020+ database
        native_ligand = deafult False, if you want native ligand diffrent than docked ligand - True
        '''
        try:
            for index,row in settings.raw_dataframe.iterrows():
                print(str(index) + '.' + row[self.pdb_id_column])
                pdb_protein_final_path = settings.datadir + '/protein/pdb'
                pdb_native_ligand_final_path = settings.datadir + '/native_ligand/pdb'
                sdf_native_ligand_final_path = settings.datadir + '/native_ligand/sdf'
                pdb_ligand_final_path = settings.datadir + '/ligand/pdb'
                sdf_ligand_final_path = settings.datadir + '/ligand/sdf'

                protein_pdb_file = settings.raw_datadir + '/' + row[self.pdb_id_column] + '/protein.pdb'
                native_ligand_pdb_file = settings.raw_datadir + '/' + row[self.pdb_id_column] + '/ligand.pdb'
                native_ligand_sdf_file = settings.raw_datadir + '/' + row[self.pdb_id_column] + '/ligand.sdf'
                ligand_pdb_file = settings.raw_datadir + '/' + row[self.pdb_id_column] + '/ligand.pdb'
                ligand_sdf_file = settings.raw_datadir + '/' + row[self.pdb_id_column] + '/ligand.sdf'

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
        '''
        Copy files from Leak Proof - Protein Data Bank Bind dataset
        lp_pdbbnd = directory of LPPDBBind database
        native_ligand = False (deafult), if you want native ligand be diffrent than docked oen -> True
        '''
        try:

            pdb_protein_final_path = settings.datadir + '/protein/pdb'
            pdb_pocket_final_path = settings.datadir + '/pocket/pdb'
            mol2_native_ligand_final_path = settings.datadir + '/native_ligand/mol2'
            sdf_native_ligand_final_path = settings.datadir + '/native_ligand/sdf'
            mol2_ligand_final_path = settings.datadir + '/ligand/mol2'
            sdf_ligand_final_path = settings.datadir + '/ligand/sdf'

            for index,row in raw_dataframe.iterrows():
                print(str(index) + '.' + row[self.pdb_id_column])

                protein_pdb_file = settings.raw_datadir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_protein.pdb'
                pocket_pdb_file = settings.raw_datadir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_pocket.pdb'
                ligand_mol2_file= settings.raw_datadir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.mol2'
                ligand_sdf_file= settings.raw_datadir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.sdf'
                native_ligand_mol2_file= settings.raw_datadir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.mol2'
                native_ligand_sdf_file= settings.raw_datadir +'/'+row[self.pdb_id_column]+'/'+row[self.pdb_id_column]+'_ligand.sdf'

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
            print('There is no directory!')
class Converter:
    def __init__(self,molecule):
        settings.init()
        self.molecule = molecule    #  specify molecule PDB ID from your dataframe
        self.converting_dictionary = {molecule:{}}
    def pdb_to_pdbqt(self,
                     protein=False,         #  following protein, pocket, ligand and native_ligand arguments are gates to convert specific molecule type from PDB complex
                     pocket=False,
                     ligand=False,
                     native_ligand=False):
        '''
        This method converts the specified molecule from .pdb into .pdbqt format
        You will need:
            * MGL Tools PackageObabel
            * OpenBabel
        '''

        self.converting_dictionary[self.molecule]['pdb_to_pdbqt']={}

        if protein==False and pocket==False and ligand == False and native_ligand == False:

            self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['error'] = 'Argument Error. There is nothing to do. Please specify in settings.py which molecule type you want to convert (protein/pocket/ligand/native_ligand)'

        if protein == True:

            pdb_protein_file = settings.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            pdbqt_protein_file = settings.datadir + '/protein/pdbqt/' + self.molecule + '_protein.pdbqt'

            self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['protein'] = {}

            if (not os.path.exists(pdbqt_protein_file) or os.path.getsize(pdbqt_protein_file) == 0) and os.path.exists(pdb_protein_file):

                command = [settings.mgltools_dir + '/bin/pythonsh', settings.mgltools_dir + '/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py',
                           '-r',pdb_protein_file,
                           '-o',pdbqt_protein_file,
                           '-A','bonds_hydrogens']

                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)

                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['protein']['mgl_tools'] = {
                    'output':pdbqt_result.stdout}

                if int(len(pdbqt_result.stderr)) > 0:

                    self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['protein']['mgl_tools'] = {
                        'error': pdbqt_result.stderr}

                    command = [settings.obabel_path,
                               pdb_protein_file,
                               '-O',pdbqt_protein_file,'-h']

                    pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)

                    self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['protein']['openbabel'] = {
                        'output': pdbqt_result.stdout}

                    if int(len(pdbqt_result.stderr)) > 0:
                        self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['protein']['openbabel'] = {
                            'error': pdbqt_result.err}

            else:
                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['protein']['error'] = 'File error. There is some problem with pdbqt file or pdb file has already exist.'

        if pocket == True:

            pdb_pocket_file = settings.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'
            pdbqt_pocket_file = settings.datadir + '/pocket/pdbqt/' + self.molecule + '_pocket.pdbqt'

            self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['pocket'] = {}

            if (not os.path.exists(pdbqt_pocket_file) or os.path.getsize(pdbqt_pocket_file) == 0) and os.path.exists(pdb_pocket_file):

                command = [settings.mgltools_dir + '/bin/pythonsh',settings.mgltools_dir + '/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py',
                           '-r', pdb_pocket_file,
                           '-o', pdbqt_pocket_file,
                           '-A', 'bonds_hydrogens']

                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)

                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['pocket']['mgl_tools'] = {'output':pdbqt_result.stdout}

                if int(len(pdbqt_result.stderr)) > 0:

                    self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['pocket']['mgl_tools'] = {
                        'error': pdbqt_result.stderr}

                    command = [settings.obabel_path,pdb_pocket_file,
                               '-O',pdbqt_pocket_file,'-h']

                    pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)

                    self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['pocket']['openbabel'] = {
                        'output': pdbqt_result.stdout}

                    if int(len(pdbqt_result.stderr)) > 0:

                        self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['pocket']['openbabel'] = {
                            'error': pdbqt_result.err}
            else:
                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['pocket']['error'] = 'File error. There is some problem with pdbqt file or pdb file has already exist.'

        if ligand == True:

            pdb_ligand_file = settings.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'
            pdbqt_ligand_file = settings.datadir + '/ligand/pdbqt/' + self.molecule + '_ligand.pdbqt'

            self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['ligand'] = {}

            if (not os.path.exists(pdbqt_ligand_file) or os.path.getsize(pdbqt_ligand_file) == 0) and os.path.exists(pdb_ligand_file):

                command = [settings.obabel_path,pdb_ligand_file,'-O',pdbqt_ligand_file,'-h']
                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)

                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['ligand']['openbabel'] = {'output':pdbqt_result.stdout}

            else:

                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['ligand']['error'] = 'File error. There is some problem with pdbqt file or pdb file has already exist.'

        if native_ligand == True:

            pdb_native_ligand_file = settings.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'
            pdbqt_native_ligand_file = settings.datadir + '/native_ligand/pdbqt/' + self.molecule + '_ligand.pdbqt'

            self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['native_ligand'] = {}


            if (not os.path.exists(pdbqt_native_ligand_file) or os.path.getsize(pdbqt_native_ligand_file) == 0) and os.path.exists(pdb_native_ligand_file):

                command = [settings.obabel_path, pdb_native_ligand_file, '-O', pdbqt_native_ligand_file,'-h']
                pdbqt_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)

                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['native_ligand']['openbabel'] = {'output':pdbqt_result.stdout}

            else:
                self.converting_dictionary[self.molecule]['pdb_to_pdbqt']['native_ligand']['error'] = 'File error. There is some problem with pdbqt file or pdb file has already exist.'

    def pdb_to_seq(self,
                   protein=True,         #  following protein, pocket arguments are gates to convert specific molecule type from PDB complex
                   pocket=False):

        self.converting_dictionary[self.molecule]['pdb_to_seq'] = {}

        if protein == False and pocket == False:
            self.converting_dictionary[self.molecule]['pdb_to_seq']['error'] = 'Argument Error. There is nothing to do. Please specify in settings.py which molecule type you want to convert (protein/pocket)'

        if protein==True:

            protein_pdb_file = settings.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            protein_fasta_file = settings.datadir + '/protein/fasta/' + self.molecule + '_protein.fasta'

            self.converting_dictionary[self.molecule]['pdb_to_seq']['protein'] = {}

            if os.path.exists(protein_pdb_file) and not os.path.exists(protein_fasta_file):

                self.converting_dictionary[self.molecule]['pdb_to_seq']['protein']['pdb2fasta'] = {}

                result = subprocess.run([settings.pdb2fasta_path+' '+protein_pdb_file+' > '+protein_fasta_file], shell=False, capture_output=True, text=True, close_fds=True)

                self.converting_dictionary[self.molecule]['pdb_to_seq']['protein']['pdb2fasta']['output'] = result.stdout
                self.converting_dictionary[self.molecule]['pdb_to_seq']['protein']['pdb2fasta']['error'] = result.stderr
            else:
                self.converting_dictionary[self.molecule]['pdb_to_seq']['protein']['error'] = 'File Error.'

        if pocket==True:
            pocket_pdb_file = settings.datadir+'/pocket/pdb/'+self.molecule+'_pocket.pdb'
            protein_fasta_file = settings.datadir+'/pocket/fasta/'+self.molecule+'_pocket.fasta'
            if os.path.exists(pocket_pdb_file) and not os.path.exists(pocket_fasta_file):
                result = subprocess.run([settings.pdb2fasta_path+' ' + pocket_pdb_file + ' > ' + pocket_fasta_file], shell=False,capture_output=True, text=True, close_fds=True)
                print(result.stderr)
    def pdb_to_mol(self,
                   protein=False,         #  following protein, pocket, ligand and native_ligand arguments are gates to convert specific molecule type from PDB complex
                   pocket=False,
                   ligand=False,
                   native_ligand=False):

        '''Convert pdb files to mol2 format'''

        if protein == True:

            pdb_protein_file = settings.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'
            mol2_protein_file = settings.datadir + '/protein/mol2/' + self.molecule + '_protein.mol2'

            if (not os.path.exists(mol2_protein_file) or os.path.getsize(mol2_protein_file) == 0) and os.path.exists(pdb_protein_file):

                command = [settings.obabel_path,pdb_protein_file,'-O',mol2_protein_file,'-h']
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stderr)

            else:

                mol2_convert_error = self.molecule+', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)

        if pocket == True:

            pdb_pocket_file = settings.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'
            mol2_pocket_file = settings.datadir + '/pocket/mol2/' + self.molecule + '_pocket.mol2'

            if (not os.path.exists(mol2_pocket_file) or os.path.getsize(mol2_pocket_file) == 0) and os.path.exists(pdb_pocket_file):

                command = [settings.obabel_path,pdb_pocket_file,'-O',mol2_pocket_file,'-h']
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stderr)

            else:

                mol2_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)

        if ligand == True:

            pdb_ligand_file = settings.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'
            mol2_ligand_file = settings.datadir + '/ligand/mol2/' + self.molecule + '_ligand.mol2'

            if (not os.path.exists(mol2_ligand_file) or os.path.getsize(mol2_ligand_file) == 0) and os.path.exists(pdb_ligand_file):
                command = [settings.obabel_path,pdb_ligand_file,'-O',mol2_ligand_file,'-h']
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stderr)

            else:

                mol2_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)

        if native_ligand == True:

            pdb_native_ligand_file = settings.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'
            mol2_native_ligand_file = settings.datadir + '/native_ligand/mol2/' + self.molecule + '_ligand.mol2'

            if (not os.path.exists(mol2_native_ligand_file) or os.path.getsize(mol2_native_ligand_file) == 0) and os.path.exists(pdb_native_ligand_file):
                command = [settings.obabel_path,pdb_native_ligand_file,'-O',mol2_native_ligand_file,'-h']
                mol2_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print(mol2_result.stdout)

            else:

                mol2_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(mol2_convert_error)
    def mol_to_pdb(self,
                   protein=False,         #  following protein, pocket, ligand and native_ligand arguments are gates to convert specific molecule type from PDB complex
                   pocket=False,
                   ligand=False,
                   native_ligand=False):
        '''To be done if needed'''

        '''Convert mol2 files to pdb format'''

        if protein == True:

            mol2_protein_file = settings.datadir + '/protein/mol2/' + self.molecule + '_protein.mol2'
            pdb_protein_file = settings.datadir + '/protein/pdb/' + self.molecule + '_protein.pdb'

            if (not os.path.exists(pdb_protein_file) or os.path.getsize(pdb_protein_file) == 0) and os.path.exists(mol2_protein_file):

                command = [settings.obabel_path,mol2_protein_file,'-O',pdb_protein_file,'-h']
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule+', There is no MOL2 file to convert or file already exists..'
                print(pdb_convert_error)

        if pocket == True:

            mol2_pocket_file = settings.datadir + '/pocket/mol2/' + self.molecule + '_pocket.mol2'
            pdb_pocket_file = settings.datadir + '/pocket/pdb/' + self.molecule + '_pocket.pdb'

            if (not os.path.exists(pdb_pocket_file) or os.path.getsize(pdb_pocket_file) == 0) and os.path.exists(mol2_pocket_file):

                command = [settings.obabel_path,mol2_pocket_file,'-O',pdb_pocket_file,'-h']
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(pdb_convert_error)

        if ligand == True:

            mol2_ligand_file = settings.datadir + '/ligand/mol2/' + self.molecule + '_ligand.mol2'
            pdb_ligand_file = settings.datadir + '/ligand/pdb/' + self.molecule + '_ligand.pdb'

            if (not os.path.exists(pdb_ligand_file) or os.path.getsize(pdb_ligand_file) == 0) and os.path.exists(mol2_ligand_file):
                command = [settings.obabel_path,mol2_ligand_file,'-O',pdb_ligand_file,'-h']
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(pdb_convert_error)

        if native_ligand == True:

            mol2_native_ligand_file = settings.datadir + '/native_ligand/mol2/' + self.molecule + '_ligand.mol2'
            pdb_native_ligand_file = settings.datadir + '/native_ligand/pdb/' + self.molecule + '_ligand.pdb'

            if (not os.path.exists(pdb_native_ligand_file) or os.path.getsize(pdb_native_ligand_file) == 0) and os.path.exists(mol2_native_ligand_file):
                command = [settings.obabel_path,mol2_native_ligand_file,'-O',pdb_native_ligand_file,'-h']
                pdb_result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(pdb_result.stderr))+' | '+pdb_result.stderr)

            else:

                pdb_convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(pdb_convert_error)

    def mol_to_sdf(self,
                   protein=False,         #  following protein, pocket, ligand and native_ligand arguments are gates to convert specific molecule type from PDB complex
                   pocket=False,
                   ligand=False,
                   native_ligand=False):
        '''To be done if needed'''

        '''Convert sdf files to mol format'''

        if protein == True:

            mol2_protein_file = settings.datadir + '/protein/mol2/' + self.molecule + '_protein.mol2'
            sdf_protein_file = settings.datadir + '/protein/sdf/' + self.molecule + '_protein.sdf'

            if os.path.exists(mol2_protein_file):

                command = [settings.obabel_path,sdf_protein_file,'-O',mol2_protein_file,'-h']
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule+', There is no MOL2 file to convert or file already exists..'
                print(convert_error)

        if pocket == True:

            mol2_pocket_file = settings.datadir + '/pocket/mol2/' + self.molecule + '_pocket.mol2'
            sdf_pocket_file = settings.datadir + '/pocket/sdf/' + self.molecule + '_pocket.sdf'

            if os.path.exists(mol2_pocket_file):

                command = [settings.obabel_path,mol2_pocket_file,'-O',sdf_pocket_file,'-h']
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(convert_error)

        if ligand == True:

            mol2_ligand_file = settings.datadir + '/ligand/mol2/' + self.molecule + '_ligand.mol2'
            sdf_ligand_file = settings.datadir + '/ligand/sdf/' + self.molecule + '_ligand.sdf'

            if os.path.exists(mol2_ligand_file):
                command = [settings.obabel_path,mol2_ligand_file,'-O',sdf_ligand_file,'-h']
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(convert_error)

        if native_ligand == True:

            mol2_native_ligand_file = settings.datadir + '/native_ligand/mol2/' + self.molecule + '_ligand.mol2'
            sdf_native_ligand_file = settings.datadir + '/native_ligand/sdf/' + self.molecule + '_ligand.sdf'

            if os.path.exists(mol2_native_ligand_file):
                command = [settings.obabel_path,mol2_native_ligand_file,'-O',sdf_native_ligand_file,'-h']
                result = subprocess.run(command, shell=False, capture_output=True, text=True, close_fds=True)
                print('error: '+str(len(result.stderr))+' | '+result.stderr)

            else:

                convert_error = self.molecule + ', There is no PDB file to convert or file already exists..'
                print(convert_error)
class Documents:
    def __init__(self):
        settings.init()
    def merged_fasta(self,df,pdb_id_column,output,protein=True,pocket=False):

        protein_merged_seq_output_file = settings.datadir + '/docs/merged_' + output + '_protein.fasta'
        pocket_merged_seq_output_file = settings.datadir + '/docs/merged_' + output + '_pocket.fasta'

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
        output_file = settings.datadir + '/docs/' + output
        command=settings.cdhit_dir+'/cd-hit -i '+fasta_file+' -o '+output_file+' -c '+str(c)+' -n '+str(n)+' -M '+str(m)+' -d '+str(d)+' -T '+str(t)
        cluster_program=subprocess.run([command],shell=False, capture_output=True, text=True)
        output_file = output_file + '.clstr'
        return output_file
    def get_diagonal_unification(self,proteins,ligands):
        '''
        Diffrence Softwares may generate diffrent errors with some complex types.
        To make safe matrix without any wrong data we have to filtered out complexes without all docking outputs.
        Please specify list of programs which u will use to create your matrix in settings.py
        '''

        complexes = [protein + '_' + ligand for protein, ligand in zip(proteins, ligands)]

        def get_program_boolean(program):
            if program == 'smina':
                smina = Smina()
                smina_log_complexes = set([file.replace('.log', '') for file in os.listdir(smina.logs_smina_dir)])
                smina_pdbqt_complexes = set([file.replace('.pdbqt', '') for file in os.listdir(smina.logs_smina_dir)])
                smina_atom_terms_complexes = set(
                    [file.replace('_atom_terms.txt', '') for file in os.listdir(smina.logs_smina_dir)])
                smina_full_outputs = list(
                    smina_pdbqt_complexes.intersection(smina_log_complexes, smina_atom_terms_complexes))
                program_boolean = [complex in smina_full_outputs for complex in complexes]

            if program == 'rxdock':
                rxdock = RxDock()
                rxdock_output_complexes = set(
                    [file.replace('.sd', '') for file in os.listdir(rxdock.rxdock_dir) if file.endswith('.sd')])
                rxdock_rmsd_complexes = set([file.replace('_rmsd.sdf', '') for file in os.listdir(rxdock.rxdock_dir) if
                                             file.endswith('_rmsd.sdf')])
                rxdock_full_outputs = list(rxdock_output_complexes.intersection(rxdock_rmsd_complexes))
                program_boolean = [complex in rxdock_full_outputs for complex in complexes]

            if program == 'gnina':
                gnina = Gnina()
                gnina_sdf_gz_complexes = set(
                    [file.replace('_atom_terms.txt', '') for file in os.listdir(gnina.sdf_gz_output_file) if
                     file.endswith('_atom_terms.txt')])
                gnina_log_complexes = set(
                    [file.replace('.log', '') for file in os.listdir(gnina.logs_gnina_dir) if file.endswith('.log')])
                gnina_atom_terms_complexes = set(
                    [file.replace('.sdf.gz', '') for file in os.listdir(gnina.atom_terms_gnina_dir) if
                     file.endswith('.sdf.gz')])
                gnina_rmsd_complexes = set(
                    [file.replace('_rmsd.txt', '') for file in os.listdir(gnina.rmsd_gnina_dir) if
                     file.endswith('_rmsd.txt')])
                gnina_full_outputs = list(
                    gnina_sdf_gz_complexes.intersection(gnina_log_complexes, gnina_atom_terms_complexes,
                                                        gnina_rmsd_complexes))
                program_boolean = [complex in gnina_full_outputs for complex in complexes]

            return program_boolean

        program_dict = {
            'complex_name': complexes
        }
        for program in settings.matrix_programs:
            program_boolean = get_program_boolean(program)
            program_dict[program] = program_boolean

        df = pd.DataFrame(program_dict)
        df = df[['protein', 'ligand']] = df['complex_name'].str.split('_', expand=True).drop(columns=['complex_name'], inplace=True)

        return df
