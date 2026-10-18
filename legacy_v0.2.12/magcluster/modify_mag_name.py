import re

def modify_mag_name(sub_record):
    pattern = re.compile(r'[Mm][am][mdns][A-Za-z0-9\-]+' )
    for feature in sub_record.features:
        ftype = feature.type
        if ftype == 'CDS':
            fpro = feature.qualifiers['product'][0]
            if 'agnetosome' in fpro:
                match = re.search(pattern, fpro)
                if match:
                    match_name = match.group()
                    # modify the protein name
                    feature.qualifiers['product'][0] = match_name
    return sub_record
    