def init():

    ##### Software Paths #####

    global station

    global pdb2fasta_path
    global mgltools_dir
    global cdhit_dir
    global obabel_path
    global smina_tools_dir
    global gnina_container
    global rxdock

    station='eden'

    if station=='inka':
        pdb2fasta_path = 'pdb2fasta'
        mgltools_dir = '/mnt/raid/mwisniewski/software/mgl'
        cdhit_dir = '/mnt/raid/mwisniewski/software/cd-hit-v4.8.1-2019-0228'
        smina_tools_dir = ''
        obabel_path = ''
    if station=='eden':
        pdb2fasta_path = ''
        mgltools_dir = '/home2/sfglab/mwisniewski/Software/mgltools'
        cdhit_dir = ''
        smina_tools_dir = '/mnt/evafs/groups/sfglab/mwisniewski/software/smina/smina'
        gnina_container = '/mnt/evafs/groups/sfglab/mwisniewski/software/gnina/gnina.sif'
        obabel_path = '/mnt/evafs/groups/sfglab/mwisniewski/software/openbabel/bin/obabel'

    ##########################

    ###### VARIABLES #########

    global number_of_models
    global pipeline
    global steps
    global docking_programs

    number_of_models = 50
    pipeline='lp_pdbbind'
    steps=['convert']
    docking_programs=[]

    if 'convert' in steps: #'convert' in steps:
        global to_pdbqt
        global to_pdb
        global to_mol
        global to_sdf

        to_pdbqt={
            'protein':False,
            'pocket': False,
            'ligand':True,
            'native_ligand':False
        }
        to_pdb={
            'protein':False,
            'pocket': False,
            'ligand':False,
            'native_ligand':False
        }
        to_mol={
            'protein':False,
            'pocket': False,
            'ligand':False,
            'native_ligand':False
        }
        to_sdf={
            'protein':False,
            'pocket': False,
            'ligand':False,
            'native_ligand':False
        }

