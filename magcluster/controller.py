#magcluster

def controller(args, subparser_name=None):
    """user args to run magcluster"""
    from subprocess import run 
    from .maga import maga
    from .magsc_new import magsc
    from .magm import magm
    
    if subparser_name:
        if subparser_name == 'prokka':
            maga(args)
        elif subparser_name == 'mgc_screen':
            magsc(args)
        elif subparser_name == 'clinker':
            magm()
    else:
        usr_arg1 = ['magcluster', '-h']
        run(usr_arg1)

