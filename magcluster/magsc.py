#magnetosome gene and protein screening

def contig_len(contig):
    len = 0
    for i in contig:
        if i.islower():
            len += 1
    return len

def magene_screen(gbkfile_path, threshold = 2, length = 2000, force=False):
    import re
    import pandas as pd
    import os

    #解析路径，创建文件夹
    mgc_folder = os.path.dirname(gbkfile_path) + '/mgc_screen/'
    if not os.path.exists(mgc_folder):
        os.mkdir(mgc_folder)
    gbkfile_prefix = os.path.basename(gbkfile_path).rstrip('.gbk')
    clean_gbk = mgc_folder + gbkfile_prefix + '_clean.gbk'
    magpro = mgc_folder + gbkfile_prefix + '_magpro.xlsx'

    #读入gbk文件
    with open(gbkfile_path,'r') as f:
        allcontents = f.read()
    allcontigs = allcontents.split('LOCUS') #使用LOCUS分隔不同contig

    mag_contigs = [] #设定包含有mag基因的contig

    #筛选含有mag基因的contig
    for contig in allcontigs:
        if contig_len(contig) >= length:
            if contig.count('agnetosome') >= threshold:
                mag_contig = 'LOCUS' + contig
                mag_contigs.append(mag_contig)
    magtext = ''.join(mag_contigs)

    #写出clean_gbk
    if os.path.exists(clean_gbk):
        if not force:
            raise FileExistsError(gbkfile_prefix + ' _clean.gbk' + 'already exists! Please change --outdir or use --force')
        # elif force:
        #     continue
    with open(clean_gbk,'w') as f:
        f.write(magtext)

    #筛选mag_pro,输出蛋白文件
    locus_tags = []
    protein_names = []
    protein_seqs = []
    lengths = []
    for mag_contig in mag_contigs:
        mag_contig_split = mag_contig.split('gene')
        for i in mag_contig_split:
            if 'agnetosome' in i:
                locus_tag = re.search(r'/locus_tag="(.+)"', i).group(1)
                locus_tags.append(locus_tag)

                protein_name = re.search(r'/product=".*?[Mm]agnetosome protein ([a-zA-Z0-9-]+)', i).group(1)
                protein_names.append(protein_name)

                protein_seq = re.search(r'/translation="([\s\w]+)"', i, re.M).group(1).replace('\n', '').replace(' ','')
                len_ = len(protein_seq)
                protein_seqs.append(protein_seq)
                lengths.append(len_)
    mag_pro_dic = {
        'tag' : locus_tags,
        'protein name' : protein_names,
        'length' : lengths,
        'sequence' : protein_seqs,
    }
    mag_df = pd.DataFrame(
            mag_pro_dic
        )

    if os.path.exists(magpro):
        if not force:
            raise FileExistsError(gbkfile_prefix + ' _magpro.xlsx' + 'already exists! Please change --outdir or use --force')
        # elif force:
        #     continue
    mag_df.to_excel(magpro, sheet_name = 'magpro', index = False)

def magsc(args):
    from .batch_proc import get_files

    gbkfiles = get_files(args.gbkfile, extensions=['*.gbk', '*.gbf'])
    for gbkfile in gbkfiles:
        magene_screen(gbkfile, threshold=args.threshold, length=args.length, force=args.force)
    print('[The protein file is screening...]')
    print("[A xlsx file named as 'magpro.xlsx' is generated.]")
    print('[The genbank file is screening...]')
    
    print("[A .gbk file named as 'XXX_clean.gbk' is produced.]")
    print('[Thank you for using magash.]')