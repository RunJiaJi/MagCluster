

def maga(args):
    from subprocess import run
    from .batch_proc import get_prokka_cmd

    prokka_cmd = get_prokka_cmd(args)
    for cmd in prokka_cmd:
        run(cmd)