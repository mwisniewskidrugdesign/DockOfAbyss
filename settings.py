def init():

    ##### Software Paths #####

    global station
    global pdb2fasta_path
    global mgltools_dir
    global cdhit_dir
    global smina_tools_dir
    global obabel_path
    global rxdock

    station='eden'

    if station=='inka':
        pdb2fasta_path = 'pdb2fasta'
        mgltools_dir = '/mnt/raid/mwisniewski/software/mgl'
        cdhit_dir = '/mnt/raid/mwisniewski/software/cd-hit-v4.8.1-2019-0228'
        smina_tools_dir = 'smina'
        obabel_path = 'obabel'


    if station=='eden':
        pdb2fasta_path = '/home2/sfglab/mwisniewski/Software/pdb2fasta'
        mgltools_dir = '/home2/sfglab/mwisniewski/Software/mgltools'
        cdhit_dir = '/home2/sfglab/mwisniewski/Software/cd-hit-v4.8.1-2019-0228'
        smina_tools_dir = 'smina'
        obabel_path = 'obabel'

    ##########################
