import os

def fasta_path():
    f = __file__
    p = os.path.dirname(f)
    return os.path.join(p, 'Magnetosome_protein_data.fasta')