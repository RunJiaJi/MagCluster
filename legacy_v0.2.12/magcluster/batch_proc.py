from .data import fasta_path
from pathlib import Path, PurePath

def get_files(paths, extensions):

    all_files = []
    for p in paths:
        p_ = Path(p)
        if p_.is_dir():
            for ext in extensions:
                all_files.extend(p_.glob(ext))
        else:
            all_files.append(p_)
    files_ = [str(i) for i in all_files]
    # tmp = []
    # rmo = []
    # for i in files_:
    #     i__ = Path(i)
    #     if i__.is_file():
    #         rmo.append(i)
    #         i_ = './' + i
    #         tmp.append(i_)
    # for i in rmo:
    #     files_.remove(i)
    # files_.extend(tmp)
    return files_

def get_prefix(file):

    prefix_ = PurePath(file).name
    extensions = ['.fa', '.fasta', '.FASTA', '.fna']
    for i in extensions:
        if i in prefix_:
            prefix = prefix_.replace(i, '')
    return prefix

def get_outdir(file):

    outdir = PurePath(file).parent.joinpath(get_prefix(file))
    return str(outdir)+'_annotation'

def get_prokka_cmd(args):

    args_dic = args.__dict__
    del args_dic['subparser_name']
    prokka_cmd_tmp = []
    fafiles_ = []
    extensions = ['*.fa', '*.fasta', '*.FASTA', '*.fna']
    # if len(args.fafile)
    for key, value in args_dic.items():
        if value:
            if value is True:
                key_tmp = '--' + key
                prokka_cmd_tmp.append(key_tmp)
            elif key == 'fafile':
                fafiles_ = get_files(value, extensions)
            else:
                key_tmp = '--' + key
                prokka_cmd_tmp.append(key_tmp)
                prokka_cmd_tmp.append(value)
    prokka_cmd_tmp.insert(0,'prokka')    
    prokka_cmd = []
    tmp = [str(i) for i in prokka_cmd_tmp]
    for i in fafiles_:
        if '--outdir' not in prokka_cmd_tmp:
            tmp.extend(['--outdir', get_outdir(i)])
        if '--prefix' not in prokka_cmd_tmp:
            tmp.extend(['--prefix', get_prefix(i)])
        if '--locustag' not in prokka_cmd_tmp:
            tmp.extend(['--locustag', get_prefix(i)])
        if '--compliant' not in prokka_cmd_tmp:
            tmp.append('--compliant')
        if '--proteins' not in prokka_cmd_tmp:
            tmp.extend(['--proteins', fasta_path()])
        prokka_cmd.append(tmp + [i])

    return prokka_cmd