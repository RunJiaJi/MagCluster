from Bio import SeqIO
import logging
from .batch_proc import get_files
import os
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S"
)
log = logging.getLogger()

def get_loci(contig):
    """Return the mag gene loci in long contig(SeqIO record)"""
    mag_loci = []
    for feature in contig.features:
        if feature.type == 'CDS':
            if 'agnetosome' in feature.qualifiers['product'][0]:
                start_loci = int(feature.location.start)
                end_loci = int(feature.location.end)
                mag_loci.append((start_loci, end_loci))
    return mag_loci

def get_gap(mag_loci):
    gap = []
    for i in range(len(mag_loci)-1):
        gap.append(mag_loci[i+1][0] - mag_loci[i][1])
    return gap

def cut_ctg(contig, windowsize):
    if len(contig) <= windowsize:
        return [contig]
    
    mag_loci = get_loci(contig)
    start, end = mag_loci[0][0], mag_loci[-1][1]
    distance = end - start
    
    if distance <= windowsize:
        return [contig[start:end]]
    
    gap = get_gap(mag_loci)
    max_gap = max(gap)
    max_gap_index = gap.index(max_gap)

    left_cut_start = start
    left_cut_end = mag_loci[max_gap_index][1]
    sub_cotig_left = contig[left_cut_start:left_cut_end]

    right_cut_start = mag_loci[max_gap_index+1][0]
    right_cut_end = mag_loci[-1][1]
    sub_cotig_right = contig[right_cut_start:right_cut_end]

    ctg_list = []
    ctg_list.extend(cut_ctg(sub_cotig_left, windowsize))
    ctg_list.extend(cut_ctg(sub_cotig_right, windowsize))
    return ctg_list

def mag_count(contig):
    mag_count = 0
    for feature in contig.features:#one gene of all genes in one contig
        if  feature.type == 'CDS' and 'agnetosome' in feature.qualifiers['product'][0]:
            mag_count += 1
    return mag_count

def mag_ctg_sc(gbk_file_path, contiglength=2000, windowsize=10000, threshold = 2, outdir=None):
    '''Parse gbk_file
    Screening for magnetosome gene containing contigs with certain minimum length.'''

    #解析路径，创建文件夹
    log.info("Starting mgc_screen...")
    if outdir:
        mgc_folder = os.path.abspath(outdir) + '/'
    else:
        mgc_folder = os.path.dirname(gbk_file_path) + '/mgc_screen/'

    gbkfile_prefix = os.path.basename(gbk_file_path).rstrip('.gbk')
    clean_gbk_outpath = mgc_folder + gbkfile_prefix + '_clean.gbk'
    magpro = mgc_folder + gbkfile_prefix + '_magpro.csv'

    log.info("Your file is " + os.path.basename(gbk_file_path))
    log.info("The minimum length of contigs to be considered is " + str(contiglength))
    log.info("The maxmum length of contigs to be considered is " + str(windowsize))
    log.info("The threshold of magnetosome genes in one contig is " + str(threshold))
    log.info("The output directory is " + mgc_folder)
    log.info("Opening your file...")

    records = SeqIO.parse(gbk_file_path, 'gb')
    MAG_ctg = []
    contigs_list_tmp = []
    log.info("Starting magnetosome genes screening...")
    for record in records:#one contig in all contigs
        if mag_count(record) >= threshold and len(record) >= contiglength:
            MAG_ctg.append(record)
            contig_list = cut_ctg(record, windowsize)
            contigs_list_tmp.extend(contig_list)
    contigs_list = [contig for contig in contigs_list_tmp if mag_count(contig) >= threshold]
    if len(contigs_list)==0:
        log.info('No magnetosome genes were found in your file!')
        return

    log.info("Magnetosome gene cluster containing contigs screening completed!")

    if not os.path.exists(mgc_folder):
        os.mkdir(mgc_folder)
        log.info("Creating output folder: " + mgc_folder)
    elif os.path.exists(mgc_folder):
        log.info(mgc_folder + " folder already exists! Skip to next step")
    log.info("Writing clean.gbk file...")
    SeqIO.write(MAG_ctg, clean_gbk_outpath, "genbank")
    # return len(contigs_list)

    locus_tags = []
    products = []
    translations = []
    magpro_lens = []
    for contig in contigs_list:
        for feature in contig.features:
            if feature.type == 'CDS' and 'agnetosome' in feature.qualifiers['product'][0]:
                locus_tag = feature.qualifiers['locus_tag'][0]
                product = feature.qualifiers['product'][0]
                translation = feature.qualifiers['translation'][0]
                magpro_len = len(translation)

                locus_tags.append(locus_tag)
                products.append(product)
                translations.append(translation)
                magpro_lens.append(magpro_len)
    mag_pro_dic = {
        'locus_tag' : locus_tags,
        'name' : products,
        'length' : magpro_lens,
        'sequence' : translations,
    }
    mag_df = pd.DataFrame(mag_pro_dic)
    log.info("Writing magpro.csv file(s)...")
    mag_df.to_csv(magpro, index = False)

def magsc(args):
    gbkfiles = get_files(args.gbkfile, extensions=['*.gbk'])
    for gbkfile in gbkfiles:
        mag_ctg_sc(gbkfile, contiglength=args.contiglength, windowsize=args.windowsize, threshold=args.threshold, outdir=args.outdir)
    log.info("Done! Thank you!")