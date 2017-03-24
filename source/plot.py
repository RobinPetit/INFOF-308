#! /usr/bin/env python3

# std
import matplotlib.pyplot as plt
from glob import glob

# local
import separation
import localization

DISEASES_LOCATION = '../data/diseases/'
SIGNIFICANCE_THRESHOLD = 1.6

def get_all_disease_files():
    """
    returns a list of disease files
    """
    return glob(DISEASES_LOCATION + '*.txt')

def plot_zscore_vs_relative_size(G, all_genes_in_network, disease_file_list,
        nb_simulations_per_disease=1000, show_non_significant=True,
        show_on_screen=False):
    """
    make a plot of relative module size on x axis vs z-score of module size on y axis
    """
    non_significant_counter = 0
    relative_size_list = list()
    zscore_list = list()
    
    nb_diseases = len(disease_file_list)
    for idx, disease_file in enumerate(disease_file_list):
        print('{} out of {}'.format(idx+1, nb_diseases))
        gene_set = separation.get_disease_genes(disease_file, all_genes_in_network)
        mean, std, zscore = localization.get_random_comparison(
                                G, gene_set, nb_simulations_per_disease)
        if zscore < SIGNIFICANCE_THRESHOLD:
            non_significant_counter += 1
            if not show_non_significant:
                continue
        relative_size = localization.get_lcc_size(G,gene_set) / len(gene_set)
        relative_size_list.append(relative_size)
        zscore_list.append(zscore)

    fig = plt.figure()
    plt.plot(relative_size_list, zscore_list, 'ro')
    plt.plot([0, 1], [SIGNIFICANCE_THRESHOLD]*2, 'b--')
    plt.xlabel('relative size = S/N_d')
    plt.ylabel('z-score of module size S')
    plt.title('z-score versus relative module size with {} simulations per disease' \
              .format(nb_simulations_per_disease))
    plt.axis([0, 1, 0, 35])
    if show_on_screen:
        plt.show
    else:
        fig.savefig('../report/images/S4.b.pdf', bbox_inches='tight')
    print('{} non significant z-scores'.format(non_significant_counter))

if __name__ == '__main__':
    # TODO: handle properly different plots and their parameters with program arguments
    interactome_path = separation.INTERACTOME_DEFAULT_PATH
    G, all_genes_in_network = separation.load_network(interactome_path)
    plot_zscore_vs_relative_size(G, all_genes_in_network, get_all_disease_files(),
                                 nb_simulations_per_disease=5000,
                                 show_non_significant=True)
