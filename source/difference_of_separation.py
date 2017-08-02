from constants import *

import matplotlib.pyplot as plt

from plot import *

import shelve

import separation

ORIGINAL_INTERATOME_RESULTS_FILE = '../data/DataS4_disease_pairs.tsv'

def get_file_name(path):
    return path[path.rfind('/')+1:path.rfind('.')]

def get_non_overlapping_disease_pairs_with_negative_sep_score():
    separations = dict()
    all_genes = separation.load_network(INTERACTOME_DEFAULT_PATH)[1]
    disease_files = get_all_disease_files()
    disease_genes = list(map(lambda path: separation.get_disease_genes(path, all_genes), disease_files))
    key = get_sep_key(len(disease_genes))
    d_A, d_AB = shelve.open(SHELVE_PATH)[key]
    for idx, gene_set_A in enumerate(disease_genes):
        gene_set_A = disease_genes[idx]
        for idx2 in range(idx+1, len(disease_genes)):
            gene_set_B = disease_genes[idx2]
            J = J_score(gene_set_A, gene_set_B)
            C = C_score(gene_set_A, gene_set_B)
            sep = d_AB[idx][idx2-idx-1] - (d_A[idx]+d_A[idx2])/2
            pair = frozenset((get_file_name(disease_files[idx]), get_file_name(disease_files[idx2])))
            separations[pair] = sep
    return separations

def get_separation_from_disease_pairs_original(disease_pairs=None):
    separations = dict()
    all_genes = separation.load_network(INTERACTOME_DEFAULT_PATH)[1]
    disease_files = get_all_disease_files()
    disease_genes = list(map(lambda path: separation.get_disease_genes(path, all_genes), disease_files))
    with open(ORIGINAL_INTERATOME_RESULTS_FILE, 'r') as sep_results_file:
        lines = sep_results_file.readlines()
        for idx, line in enumerate(lines):
            print('{:2.1f}%'.format(100*idx/(len(lines)-1)), end='\r')
            if line.startswith('#'):
                continue
            columns = line.split('\t')
            pair = frozenset(map(lambda s: s.replace(',', ''),columns[:2]))
            if disease_pairs is None or pair in disease_pairs:
                separations[pair] = float(columns[2])
    print('')
    return separations

if __name__ == '__main__':
    differences = list()
    separations = get_non_overlapping_disease_pairs_with_negative_sep_score()
    original_separations = get_separation_from_disease_pairs_original(separations.keys())
    for pair in separations:
        sep = separations[pair]
        original_sep = original_separations[pair]
        differences.append(abs(original_sep-sep))
        pair = tuple(pair)
        if sep < 0 and original_sep >= 0 or original_sep < 0 and sep >= 0:
            print('difference = {: .3g}\tcouple == {}'.format(differences[-1],pair))
    print('for differences: max == {}'.format(max(differences)))
    plt.title(r'Distribution of $s_{AB}$ separation scores between' + '\nprovided results and original interactome')
    plt.xlabel(r'$\left|s_{AB}^{\mathrm{provided}} - s_{AB}^{\mathrm{computed}}\right|$')
    plt.ylabel('Number of disease pairs')
    plt.hist(differences, 50, color='tomato')
    plt.yscale('log')
    plt.ylim(.9, plt.ylim()[1])
    plt.show()
