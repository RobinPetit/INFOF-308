#! /usr/bin/env python3

# std
import matplotlib.pyplot as plt
from glob import glob

import optparse

# local
import separation
import localization
import math_treatment as math

DISEASES_LOCATION = '../data/diseases/'
SIGNIFICANCE_THRESHOLD = 1.6

VERBOSE = False

def set_verbose(option=None, opt=None, value=None, parser=None):
    global VERBOSE
    VERBOSE = True

def print_usage(option=None, opt=None, value=None, parser=None):
    """
    Displays the usage of the script
    TODO
    """
    print("USAGE...")
    exit()

def get_program_arguments():
    """
    parses the arguments given in command line and return them all
    """
    parser = optparse.OptionParser()

    parser.add_option('-u', '--usage',
                      help    ='print more info on how to use this script',
                      action  ="callback",
                      callback=print_usage)
                      
    parser.add_option('-v', '--verbose',
                      help    ='displays details of the steps',
                      action  ='callback',
                      callback=set_verbose)

    parser.add_option('-o', '--output',
                      help    ='output file for the graphs [on screen]',
                      dest    ='output_file',
                      default ='',
                      type    ="string")

    (opts, args) = parser.parse_args()
    return opts

def get_all_disease_files():
    """
    returns a list of disease files
    """
    return glob(DISEASES_LOCATION + '*.txt')

def get_relative_module_size_and_zscore(G, all_genes_in_network,
        disease_file_list, return_non_significant, significance_threshold,
        nb_simulations_per_disease):
    """
    returns a list of relative module sizes and a list of zscores of these module sizes
    
    PARAMETERS:
    -----------
        + G: the interactome network
        + all_genes_in_network: a set of all the genes in G
        + disease_file_list: a list of paths to all desired diseases to analyze
        + return_non_significant: a boolean telling whether non-significant
          modules should be returned or not
        + significance_threshold: a positive value representing the smallest
          z-score considered as significant (in absolute value)
        + nb_simulations_per_disease: the number of simulations to make for each
          disease in order to get the z-score
    
    RETURNS:
    --------
        + relative_size_list: a list of relative module sizes
        + zscore_list: a list of associated z-scores
        + non_significant_counter: number of diseases which don't show
          significant z-scores
    """
    non_significant_counter = 0
    relative_size_list = list()
    zscore_list = list()
    
    nb_diseases = len(disease_file_list)
    for idx, disease_file in enumerate(disease_file_list):
        if VERBOSE:
            print('\r{} out of {}'.format(idx+1, nb_diseases), end='')
        gene_set = separation.get_disease_genes(disease_file, all_genes_in_network)
        mean, std, zscore = localization.get_random_comparison(
                                G, gene_set, nb_simulations_per_disease)
        if zscore < significance_threshold:
            non_significant_counter += 1
            if not return_non_significant:
                continue
        relative_size = localization.get_lcc_size(G,gene_set) / len(gene_set)
        relative_size_list.append(relative_size)
        zscore_list.append(zscore)
    if VERBOSE:
        print('\n{} non significant z-scores'.format(non_significant_counter))
    return relative_size_list, zscore_list, non_significant_counter
    
def plot_linear_regression(x, y, xlabel='x', ylabel='y'):
    """
    makes a linear regression of given data
    """
    # get regression parameters
    m, p = math.linear_regression(x, y)
    # determine plot label
    if xlabel != '' and ylabel != '':
        label = '{} = {:.3f} {} + {:.3f}'.format(ylabel, m, xlabel, p)
    else:
        label = ''
    # plot the regression
    regression_line = plt.plot([0, 1], [p, m+p], label=label)[0]
    # display label if exists
    if label != '':
        plt.legend(handles=[regression_line])

def display_plot(path):
    """
    """
    if path == '':
        plt.show()
    else:
        fig.savefig(path, bbox_inches='tight')

def plot_zscore_vs_relative_size(G, all_genes_in_network, disease_file_list,
        nb_simulations_per_disease=1000, show_non_significant=True,
        show_linear_regression=True, path=''):
    """
    make a plot of relative module size on x axis vs z-score of module size on y axis
    
    PARAMETERS:
    -----------
        + G: the interactome network
        + all_genes_in_network: a set of all the genes in G
        + disease_file_list: a list of paths to all desired diseases to plot
        + nb_simulations_per_disease: the number of simulations used to determine
          the z-score of a module size
        + show_non_significant: a boolean telling whether non-significant module
          sizes should be displayed with the other (with a mark) or simply hidden
        + path: the path of of the pdf to save the plot in (empty string to
          display on screen
          
    RETURNS:
    --------
        + None
    """
    relative_size_list, zscore_list, non_significant_counter = \
        get_relative_module_size_and_zscore(G, all_genes_in_network,
            disease_file_list, show_non_significant, SIGNIFICANCE_THRESHOLD,
            nb_simulations_per_disease)

    fig = plt.figure()
    # dot plot for each disease
    plt.plot(relative_size_list, zscore_list, 'ro')
    # plot significance threshold
    if show_non_significant:
        plt.plot([0, 1], [SIGNIFICANCE_THRESHOLD]*2, 'b--')
    # plot properties
    plt.xlabel('relative size = S/N_d')
    plt.ylabel('z-score of module size S')
    plt.title('z-score versus relative module size with {} simulations per disease' \
              .format(nb_simulations_per_disease))
    plt.axis([0, 1, 0, 35])
    if show_linear_regression:
        plot_linear_regression(relative_size_list, zscore_list)
    display_plot(path)

if __name__ == '__main__':
    opts = get_program_arguments()
    interactome_path = separation.INTERACTOME_DEFAULT_PATH
    G, all_genes_in_network = separation.load_network(interactome_path)
    plot_zscore_vs_relative_size(G, all_genes_in_network, get_all_disease_files(),
                                 nb_simulations_per_disease=100,
                                 show_non_significant=True,
                                 show_linear_regression=True,
                                 path=opts.output_file)

