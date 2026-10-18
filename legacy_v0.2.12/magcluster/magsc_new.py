from Bio import SeqIO
import logging
from .batch_proc import get_files
from .slice_contig import slice_gbk
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

def cut_ctg(contig, windowsize, threshold):
    if len(contig) <= windowsize:
        return [contig]
    
    mag_loci = get_loci(contig)
    start, end = mag_loci[0][0], mag_loci[-1][1]
    distance = end - start
    
    if distance <= windowsize:
        return [contig[start:end]]
    
    if len(mag_loci)>= threshold:
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
        if len(get_loci(sub_cotig_left))>= threshold:
            ctg_list.extend(cut_ctg(sub_cotig_left, windowsize, threshold))
        if len(get_loci(sub_cotig_right))>= threshold:   
            ctg_list.extend(cut_ctg(sub_cotig_right, windowsize, threshold))
            
        return ctg_list

def mag_count(contig):
    mag_count = 0
    for feature in contig.features:#one gene of all genes in one contig
        if 'product' in feature.qualifiers:
            if  feature.type == 'CDS' and 'agnetosome' in feature.qualifiers['product'][0]:
                mag_count += 1
    return mag_count

def mag_ctg_sc(gbk_file_path, contiglength=2000, windowsize=10000, threshold=3, outdir=None, unmag_num=3):
    '''Parse gbk_file
    Screening for magnetosome gene containing contigs with certain minimum length.'''

    #解析路径，创建文件夹
    log.info("Starting mgc_screen...")
    if outdir:
        mgc_folder = os.path.abspath(outdir) + '/'
    else:
        dir_path = os.path.dirname(gbk_file_path)
        if len(dir_path) == 0:
            mgc_folder = '.' + '/mgc_screen/'
        else:
            mgc_folder = dir_path + '/mgc_screen/'

    gbkfile_prefix = os.path.basename(gbk_file_path).rstrip('.gbk')
    tmp_gbk_outpath = mgc_folder +'tmp/'+ gbkfile_prefix + '_mgc.gbk'
    slice_gbk_outpath = mgc_folder + gbkfile_prefix + '_sliced_mgc.gbk'
    magpro = mgc_folder + gbkfile_prefix + '_magpro.csv'

    log.info("Your file is " + os.path.basename(gbk_file_path))
    log.info("The minimum length of contigs to be considered is " + str(contiglength))
    log.info("The window size for MGCs screening is " + str(windowsize))
    log.info("The threshold of magnetosome genes in a given screening window is " + str(threshold))
    log.info("The output directory is " + mgc_folder)
    log.info("Opening your file...")

    records = SeqIO.parse(gbk_file_path, 'gb')
    MAG_ctg = []
    contigs_list_tmp = []
    log.info("Starting magnetosome genes screening...")
    for record in records:#one contig in all contigs
        if mag_count(record) >= threshold and len(record) >= contiglength:
            # MAG_ctg.append(record)
            contig_list = cut_ctg(record, windowsize, threshold)
            if contig_list:
                MAG_ctg.append(record)
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
    if not os.path.exists(mgc_folder+'tmp/'):
        os.mkdir(mgc_folder+'tmp/')
    log.info("Writing mgc.gbk file...")
    SeqIO.write(MAG_ctg, tmp_gbk_outpath, "genbank")
    
    # return all proteins in a putative MGCs containing contig inferred above
    tmp_gbk_faa_outpath = mgc_folder +'tmp/'+ gbkfile_prefix + '_mgc.faa'
    for ctg in MAG_ctg:
        pros = []
        for f in ctg.features:
            if f.type =='CDS':
                pro = '>' + f.qualifiers['locus_tag'][0]+' '+f.qualifiers['product'][0] + '\n' +f.qualifiers['translation'][0]
                pros.append(pro)
        with open (tmp_gbk_faa_outpath, 'w') as f:
            pros_str = '\n'.join(pros)
            f.write(pros_str)
            
    # return sliced genbank file
    list_sub_records = slice_gbk(MAG_ctg, slice_gbk_outpath, unmag_num)

    # return all proteins in sliced genbank file
    slice_gbk_faa_outpath = mgc_folder + gbkfile_prefix + '_sliced_mgc.faa'
    for ctg in list_sub_records:
        pros = []
        for f in ctg.features:
            if f.type =='CDS':
                pro = '>' + f.qualifiers['locus_tag'][0]+' '+f.qualifiers['product'][0] + '\n' +f.qualifiers['translation'][0]
                pros.append(pro)
        with open (slice_gbk_faa_outpath, 'w') as f:
            pros_str = '\n'.join(pros)
            f.write(pros_str)
    
    # return the predicted putative magnetosome proteins' csv file
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

    # return the predicted putative magnetosome proteins
    magpros = []
    for i in range(len(locus_tags)):
        pro = '>' + locus_tags[i] + ' ' + products[i] + '\n' + translations[i]
        magpros.append(pro)
    magpro_faa = mgc_folder + gbkfile_prefix + '_magpro.faa'
    with open(magpro_faa, 'w') as f:
        magpros_str = '\n'.join(magpros)
        f.write(magpros_str)


def magsc(args):
    gbkfiles = get_files(args.gbkfile, extensions=['*.gbk'])
    for gbkfile in gbkfiles:
        mag_ctg_sc(gbkfile, contiglength=args.contiglength, windowsize=args.windowsize, threshold=args.threshold, outdir=args.outdir, unmag_num=args.unmag_num)
    log.info("Done! Thank you!")

