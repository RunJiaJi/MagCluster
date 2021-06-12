#magcluster

def controller(args, subparser_name=None):
    """user args to run magcluster"""
    from subprocess import run 
    from .maga import maga
    from .magsc import magsc
    from .magm import magm
    
    if subparser_name:
        if subparser_name == 'maga':
            maga(args)
        elif subparser_name == 'magsc':
            magsc(args)
        elif subparser_name == 'magm':
            magm()
    else:
        usr_arg1 = ['magcluster', '-h']
        run(usr_arg1)

