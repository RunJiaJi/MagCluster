#magnetosome gene and protein screening

def magene_screen(gbkfile_path, limit = 1):
    import re
    import pandas as pd

    #读入gbk文件
    with open(gbkfile_path,'r') as f:
        allcontents = f.read()
    allcontigs = allcontents.split('LOCUS') #使用LOCUS分隔不同contig

    mag_contigs = [] #设定包含有mag基因的contig

    #筛选含有mag基因的contig
    for contig in allcontigs:
        if contig.count('Magnetosome') >= limit:
            mag_contig = 'LOCUS' + contig
            mag_contigs.append(mag_contig)

    #使用join将所有mag_contigs整合成一个string/text,写出magene.gbk文件
    magtext = ''.join(mag_contigs)    
    clean_gbk = gbkfile_path.rstrip('.gbk')+'_clean.gbk'
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
            if 'Magnetosome' in i:
                locus_tag = re.search(r'/locus_tag="(.+)"', i).group(1)
                locus_tags.append(locus_tag)

                protein_name = re.search(r'/product="(.+) Magnetosome protein ([a-zA-Z0-9-]+)', i).group(2)
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
    mag_df.to_excel('magpro.xlsx', sheet_name = 'magpro', index = False)

def magsc(args):
    print('[The protein file is screening...]')
    print("[A xlsx file named as 'magpro.xlsx' is generated.]")
    print('[The genbank file is screening...]')
    magene_screen(gbkfile_path = args.gbkfile, limit=args.threshold)
    print("[A .gbk file named as 'XXX_clean.gbk' is produced.]")
    print('[Thank you for using magash.]')