def init():
    ##### Software Paths #####

    global station
    global datadir
    global raw_datadir
    global raw_dataframe
    global pdb2fasta_path
    global mgltools_dir
    global cdhit_dir
    global obabel_path
    global smina_tools_dir
    global gnina_container
    global rxdock
    global diffdock_script_directory

    station = 'eden'

    if station == 'inka':
        ###  data directorypaths

        datadir = '/mnt/raid/mwisniewski/PhD/data/lppdbbind/'
        raw_datadir = ''

        ###  dataframe csv filepaths

        raw_dataframe = '/mnt/raid/mwisniewski/PhD/data/LP-PDBBind/dataset/LP_PDBBind.csv'

        ###  software paths

        pdb2fasta_path = 'pdb2fasta'
        mgltools_dir = '/mnt/raid/mwisniewski/software/mgl'
        cdhit_dir = '/mnt/raid/mwisniewski/software/cd-hit-v4.8.1-2019-0228'
        smina_tools_dir = ''
        obabel_path = ''

    if station == 'eden':
        ###  data directorypaths

        datadir = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/lp_pdbbind'
        raw_datadir = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/raw/PDBbind'

        ###  dataframe csv filepaths

        raw_dataframe = '/mnt/evafs/groups/sfglab/mwisniewski/PhD/data/dataframes/LP_PDBBind.csv'

        ###  software paths
        pdb2fasta_path = ''
        mgltools_dir = '/mnt/evafs/groups/sfglab/mwisniewski/software/mgltools'
        cdhit_dir = ''
        smina_tools_dir = '/mnt/evafs/groups/sfglab/mwisniewski/software/smina/smina'
        gnina_container = '/mnt/evafs/groups/sfglab/mwisniewski/software/gnina/gnina.sif'
        obabel_path = '/mnt/evafs/groups/sfglab/mwisniewski/software/openbabel/bin/obabel'
        diffdock_script_directory = '/mnt/evafs/groups/sfglab/mwisniewski/software/diffdock'

    ##########################

    ###### VARIABLES #########

    global number_of_models
    global pipeline
    global steps
    global docking_programs
    global matrix_programs

    number_of_models = 50
    pipeline = 'lp_pdbbind'
    steps = ['docking']  # list of steps to run in your pipeline in list format: generate_library, convert, docking, matrix
    docking_programs = ['diffdock']  # list of docking programs to run in your pipeline
    matrix_programs = []  # list of matrices to generate in your pipeline for specific programs in string format

    if 'convert' in steps:
        global to_pdbqt
        global to_pdb
        global to_mol
        global to_sdf

        to_pdbqt = {
            'protein': True,
            'pocket': False,
            'ligand': True,
            'native_ligand': True
        }
        to_pdb = {
            'protein': False,
            'pocket': False,
            'ligand': False,
            'native_ligand': False
        }
        to_mol = {
            'protein': False,
            'pocket': False,
            'ligand': False,
            'native_ligand': False
        }
        to_sdf = {
            'protein': False,
            'pocket': False,
            'ligand': False,
            'native_ligand': False
        }
    if 'matrix' in steps:
        matrix_programs=[]
