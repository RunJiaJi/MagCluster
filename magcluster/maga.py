

def maga():
    from .capture_args import capture_args
    from subprocess import run

    usr_args = capture_args()
    del usr_args[0:2]
    usr_args.insert(0, 'prokka')
    if '--outdir' not in usr_args:
        usr_args.append('--outdir')
        usr_args.append('maga_annotation')
    if '--prefix' not in usr_args:
        usr_args.append('--prefix')
        usr_args.append('maga_')
    run(usr_args)