#magcluster
from .change_data import change_data, get_html_path
from subprocess import run 
from .maga import maga
from .magsc_new import magsc
from .magm import magm, get_clinker_cmd

def controller(args, subparser_name=None):
    """user args to run magcluster"""

    
    if subparser_name:
        if subparser_name == 'prokka':
            maga(args)
        elif subparser_name == 'mgc_screen':
            magsc(args)
        elif subparser_name == 'clinker':
            usr_args = get_clinker_cmd(args)
            magm(usr_args)
            html_path = get_html_path(usr_args)
            change_data(html_path)
    else:
        usr_arg1 = ['magcluster', '-h']
        run(usr_arg1)

