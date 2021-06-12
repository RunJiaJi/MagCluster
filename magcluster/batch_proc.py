


def get_files(paths, extensions):
    from pathlib import Path
    all_files = []
    for p in paths:
        p_ = Path(p)
        if p_.is_dir():
            for ext in extensions:
                all_files.extend(p_.glob('**/' + ext))
        else:
            all_files.append(p_)
    files_ = [str(i) for i in all_files]
    return files_

def get_prefix(file):
    from pathlib import PurePath
    prefix = PurePath(file).name
    return prefix

def get_outdir(file):
    from pathlib import PurePath
    outdir = PurePath(file).parent.joinpath(get_prefix(file))
    return str(outdir)+'_annotation'

def get_prokka_cmd(args):
    args_dic = args.__dict__
    del args_dic['subparser_name']
    prokka_cmd_tmp = []
    fafile_arg = []
    extensions = ['*.fa', '*.fasta', '*.FASTA']
    # if len(args.fafile)
    for key, value in args_dic.items():
        if value:
            if value == True:
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
    for i in fafiles_:
        tmp = prokka_cmd_tmp.copy()
        if '--outdir' not in prokka_cmd_tmp:
            tmp.extend(['--outdir', get_outdir(i)])
        if '--prefix' not in prokka_cmd_tmp:
            tmp.extend(['--prefix', get_prefix(i)])
        prokka_cmd.append(tmp + [i])
        
        
    return prokka_cmd