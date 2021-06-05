

def magm():
    from .capture_args import capture_args
    from subprocess import run

    usr_args = capture_args()
    del usr_args[0:2]
    usr_args.insert(0, 'clinker')
    run(usr_args)