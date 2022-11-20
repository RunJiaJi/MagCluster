from Bio import SeqIO
from itertools import groupby
from operator import itemgetter
from .modify_mag_name import modify_mag_name

def get_continous_slice(set_mag_idx):
    list_mag_idx = list(set_mag_idx)
    list_mag_idx.sort()
    list_continous_slice = []
    for k, g in groupby(enumerate(list_mag_idx), lambda ix : ix[0] - ix[1]):
        list_continous_slice.append(list(map(itemgetter(1), g)))
    return list_continous_slice

def get_slice_position(list_continous_slice, list_CDS):
    slice_positions = []
    for slice_ in list_continous_slice:
        start_idx = slice_[0]
        start_cds = list_CDS[start_idx]
        start_position = start_cds.location.start.position

        end_idx = slice_[-1]
        end_cds = list_CDS[end_idx]
        end_position = end_cds.location.end.position
        
        posi_pair = (start_position, end_position)
        slice_positions.append(posi_pair)
    return slice_positions

def get_mag_idx(mag_cds_idx, unmag_num):
    return set(range(
        mag_cds_idx - unmag_num, 
        mag_cds_idx + unmag_num + 1
        ))

def find_boundary(record, unmag_num):
    """This function returns the sliced boundary of the input contig"""
    # find magnetosome CDS and extend two non-magnetosome CDS 
    list_features = record.features
    list_CDS = []
    for feature in list_features:
        if feature.type == 'CDS':
            list_CDS.append(feature)
    set_mag_idx = set()
    for cds in list_CDS:
        if 'agnetosome' in cds.qualifiers['product'][0]:# find all magnetosome CDS in the list_CDS
            idx = list_CDS.index(cds)
            set_tmp = get_mag_idx(idx, unmag_num)
            set_mag_idx.update(set_tmp)
    for i in range(unmag_num):
        if len(list_CDS)+i in set_mag_idx:
            set_mag_idx.remove( len(list_CDS)+i )
        if -1 - i in set_mag_idx:
            set_mag_idx.remove( -1-i )
    
    list_continous_slice = get_continous_slice(set_mag_idx)
    slice_positions = get_slice_position(list_continous_slice, list_CDS)
    return slice_positions

def slice_gbk(list_MAG_ctg, slice_gbk_outpath, unmag_num):
    """This function aims to slice the genbank file to extract the MTCs containing part for latter clinker ploting"""

    # slice every contig and concat all sliced contig into a sub_records list
    list_sub_records = []
    for record in list_MAG_ctg:
        list_slice_positions = find_boundary(record, unmag_num)
        for slice_ in list_slice_positions:
            sub_record = record[slice_[0]:slice_[1]]
            renamed_sub_record = modify_mag_name(sub_record)
            list_sub_records.append(renamed_sub_record)
    
    SeqIO.write(list_sub_records, slice_gbk_outpath, 'genbank')

    return list_sub_records
