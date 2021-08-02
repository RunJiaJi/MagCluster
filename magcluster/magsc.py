#magnetosome gene and protein screening

import logging

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S"
)
log = logging.getLogger()

def contig_len(contig):
    len = 0
    for i in contig:
        if i.islower():
            len += 1
    return len

def magene_screen(gbkfile_path, threshold = 2, length = 2000, force=False, outdir=None):
    import re
    import pandas as pd
    import os

    #解析路径，创建文件夹
    log.info("Starting mgc_screen...")
    if outdir:
        mgc_folder = os.path.abspath(outdir) + '/'
    else:
        mgc_folder = os.path.dirname(gbkfile_path) + '/mgc_screen/'

    gbkfile_prefix = os.path.basename(gbkfile_path).rstrip('.gbk')
    clean_gbk = mgc_folder + gbkfile_prefix + '_clean.gbk'
    magpro = mgc_folder + gbkfile_prefix + '_magpro.xlsx'

    #读入gbk文件
    log.info("Your file is " + os.path.basename(gbkfile_path))
    log.info("The minimum length of contigs to be considered is " + str(length))
    log.info("The threshold of magnetosome genes in one contig is " + str(threshold))
    log.info("The output directory is " + mgc_folder)
    log.info(" Opening your file...")
    with open(gbkfile_path,'r') as f:
        allcontents = f.read()
    allcontigs = allcontents.split('LOCUS') #使用LOCUS分隔不同contig

    mag_contigs = [] #设定包含有mag基因的contig

    #筛选含有mag基因的contig

    log.info("Starting magnetosome genes screening...")
    for contig in allcontigs:
        if contig_len(contig) >= length:
            if contig.count('agnetosome') >= threshold:
                mag_contig = 'LOCUS' + contig
                mag_contigs.append(mag_contig)
    magtext = ''.join(mag_contigs)
    
    #判断是否为空
    if len(magtext) == 0:
        log.info('No magnetosome genes were found in your file!')
        return False

    log.info("Magnetosome gene cluster containing contigs screening completed!")

    if not os.path.exists(mgc_folder):
        os.mkdir(mgc_folder)
        log.info("Creating output folder: " + mgc_folder)
    elif os.path.exists(mgc_folder):
        log.info(mgc_folder + " folder already exists! Skip to next step")

    log.info("Writing clean.gbk file(s)...")
   
    #写出clean_gbk
    if os.path.exists(clean_gbk):
        if not force:
            log.info(clean_gbk + ' already exists! Please change --outdir or use --force')
            raise SystemExit

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

                protein_name_tmp = re.search(r'/product="([\s\S]+?)"', i).group(1)
                protein_name = re.sub(' +', ' ', protein_name_tmp).replace('\n', '')
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
    log.info("Writing magpro.xlsx file(s)...")
    if os.path.exists(magpro):
        if not force:
            log.info(magpro + ' already exists! Please change --outdir or use --force')
            raise SystemExit
        # elif force:
        #     continue
    mag_df.to_excel(magpro, sheet_name = 'magpro', index = False)
    
    
def magsc(args):
    from .batch_proc import get_files

    gbkfiles = get_files(args.gbkfile, extensions=['*.gbk', '*.gbf'])
    for gbkfile in gbkfiles:
        magene_screen(gbkfile, threshold=args.threshold, length=args.length, force=args.force, outdir=args.outdir)
    log.info("Done! Thank you!")