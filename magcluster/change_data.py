import re
import json
import webbrowser

def gene_class(gene_name):
    if 'mam' in gene_name or 'Mam' in gene_name:
        return 'mam'
    if 'mms' in gene_name or 'Mms' in gene_name:
        return "mms"
    if 'mad' in gene_name or 'Mad' in gene_name:
        return 'mad'
    if 'man' in gene_name or 'Man' in gene_name:
        return 'man'
    if 'iron transporter' in gene_name:
        return 'iron metabolism'
    return None


def colour(gene_name):
    if 'mam' == gene_class(gene_name):
        return '#EA4335'
    if 'mms' == gene_class(gene_name):
        return "#FBBC05"
    if 'mad' == gene_class(gene_name):
        return '#34A853'
    if 'man' == gene_class(gene_name):
        return '#4385F4'
    if 'iron metabolism':
        return '#ffff3f'
    return '#FFFF00'

def get_html_path(usr_args):
    if '-p' in usr_args:
        html_path = usr_args[usr_args.index('-p')+1]
    if '--plot' in usr_args:
        html_path = usr_args[usr_args.index('--plot')+1]
    return html_path     

def get_groups(i):
    # groups = []
    # for i in ['mam', 'mms', 'mad', 'man']:
    g = {
        'uid': 'magnetosome_gene_' + i,
        'label': i,
        'genes': [],
        'hidden': False,
        'colour': colour(i)
    }
        # groups.append(g)
    return g

def change_data(html_path):
    with open(html_path, 'r') as f:
        html = f.read()
        data_str = re.search(r'<script>const data=(.*?);function serialise', html).group(1)
        
    data = json.loads(data_str)

    g_mam = get_groups('mam')
    g_mms = get_groups('mms')
    g_mad = get_groups('mad')
    g_man = get_groups('man')
    g_feo = get_groups('iron metabolism')

    #每一个contig
    for genome in data['clusters']:
        for contig in genome['loci']:
            for gene in contig['genes']:
                product = gene['names']['product']
                # if 'agnetosome' in product:
                gc = gene_class(product)
                if gc == 'mam':
                    gene['label'] = product
                    g_mam['genes'].append(gene['uid'])
                elif gc == 'mms':
                    gene['label'] = product
                    g_mms['genes'].append(gene['uid'])
                elif gc == 'mad':
                    gene['label'] = product
                    g_mad['genes'].append(gene['uid'])
                elif gc == 'man':
                    gene['label'] = product
                    g_man['genes'].append(gene['uid'])
                elif gc == 'iron metabolism':
                    gene['label'] = product
                    g_feo['genes'].append(gene['uid'])
                # gene['colour'] = colour(product)
                # gene['label'] = product
                else:
                    gene['label'] = ' '
    groups = [g_mam, g_mms, g_mad, g_man, g_feo]
    data['groups'] = groups
    html = html.replace(data_str, json.dumps(data))

    with open(html_path, 'w') as f:
        f.write(html)
    
    webbrowser.open(html_path)