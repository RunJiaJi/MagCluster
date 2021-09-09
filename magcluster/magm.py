from .capture_args import capture_args
from subprocess import run

def get_clinker_cmd(args):

    args_dic = args.__dict__

    usr_args = capture_args()
    del usr_args[0:2]
    usr_args.insert(0, 'clinker')

    if args_dic['plot']:
        if args_dic['plot'] == True: #user have -p but no file.html
            if '-p' in usr_args:
                usr_args.insert(usr_args.index('-p')+1, 'genome_align.html')
            if '--plot' in usr_args:
                usr_args.insert(usr_args.index('--plot')+1, 'genome_align.html')
    else:
        usr_args.extend(['-p', 'genome_align.html'])
    return usr_args

def magm(usr_args):
    run(usr_args)

