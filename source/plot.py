#! /usr/bin/env python3

# std
import matplotlib.pyplot as plt
from glob import glob
from time import time, sleep
from multiprocessing import Process, Array, Value
from math import isnan

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

    parser.add_option('-p', '--plot',
                      help   ='the name of the plot to create [zscore_vs_rel_size]',
                      dest   ='plot_name',
                      default='zscore_vs_rel_size',
                      type   ='string')
    parser.add_option('-u', '--usage',
                      help    ='print more info on how to use this script',
                      action  ='callback',
                      callback=print_usage)

    parser.add_option('-v', '--verbose',
                      help    ='displays details of the steps',
                      action  ='callback',
                      callback=set_verbose)

    parser.add_option('-o', '--output',
                      help   ='output file for the graphs [on screen]',
                      dest   ='output_file',
                      default='',
                      type   ='string')

    (opts, args) = parser.parse_args()
    return opts

def get_all_disease_files():
    """
    returns a list of disease files
    """
    return sorted(glob(DISEASES_LOCATION + '*.txt'))

def zscore_process(process_id, counter, relative_size_list, zscore_list, non_significant_counter, all_diseases,
        all_genes_in_network, G, nb_simulations_per_disease, significance_threshold, return_non_significant):
    print('process {} is starting'.format(process_id))
    nb_steps = len(relative_size_list)
    while True:
        with counter.get_lock():
            n = counter.value
            if n > nb_steps:
                break
            counter.value += 1
        print('{} out of {}'.format(n, nb_steps))
        gene_set = separation.get_disease_genes(all_diseases[n-1], all_genes_in_network)
        mean, std, zscore = localization.get_random_comparison(G, gene_set, nb_simulations_per_disease)
        if zscore < significance_threshold:
            non_significant_counter.value += 1
            if not return_non_significant:
                continue
        relative_size = localization.get_lcc_size(G,gene_set) / len(gene_set)
        relative_size_list[n-1] = relative_size
        zscore_list[n-1] = zscore
    print('process {} is finishing'.format(process_id))

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
    nb_diseases = len(disease_file_list)
    non_significant_counter = Value('i', 0)
    relative_size_list = Array('d', [float('nan')] * nb_diseases)
    zscore_list = Array('d', [float('nan')] * nb_diseases)
    counter = Value('i', 1)

    NB_PROCESSES = 4

    processes = list()
    for k in range(NB_PROCESSES):
        processes.append(Process(target=zscore_process, args=(k+1, counter, relative_size_list, zscore_list,
                non_significant_counter, disease_file_list, all_genes_in_network, G,
                nb_simulations_per_disease, significance_threshold, return_non_significant)))
        processes[-1].start()
    for process in processes:
        process.join()
    if VERBOSE:
        print('\n{} non significant z-scores'.format(non_significant_counter.value))
    relative_size_list = [el for el in relative_size_list if not isnan(el)]
    zscore_list = [el for el in zscore_list if not isnan(el)]
    return relative_size_list, zscore_list, non_significant_counter.value

def plot_linear_regression(x, y, xlabel='x', ylabel='y'):
    """
    makes a linear regression of given data
    """
    # get regression parameters
    m, p = math.linear_regression(x, y)
    # determine plot label
    if xlabel != '' and ylabel != '':
        label = '{} = {:.3f} {} {} {:.3f}'.format(ylabel, m, xlabel, '+' if p >= 0 else '-', abs(p))
    else:
        label = ''
    # plot the regression
    regression_line = plt.plot([0, 1], [p, m+p], label=label)[0]
    # display label if exists
    if label != '':
        plt.legend(handles=[regression_line], loc='upper left')

def display_plot(path, fig):
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
        significance_line = plt.plot([0, 1], [SIGNIFICANCE_THRESHOLD]*2, 'b--', label=r'significance threshold: $z$-score = ' + str(SIGNIFICANCE_THRESHOLD))[0]
        legend = plt.legend(handles=[significance_line], loc='lower right')
        plt.gca().add_artist(legend)
    # plot properties
    plt.xlabel(r'relative size: $r = S/N_d$')
    plt.ylabel(r'$z$-score of module size $S$')
    plt.title(r'$z$-score' + ' versus relative module size with {} simulations per disease' \
              .format(nb_simulations_per_disease))
    #plt.axis([0, 1, 0, 35])
    if show_linear_regression:
        plot_linear_regression(relative_size_list, zscore_list, xlabel=r'$r$', ylabel=r'$z$-score')
    plt.annotate('{} diseases in non-significant area'.format(non_significant_counter), xy=(.2, -.5), xytext=(.2, -.5))
    display_plot(path, fig)

##### OVERLAPPING

def J_score(A, B):
    return len(A & B) / len(A | B)

def C_score(A, B):
    return len(A & B) / min(len(A), len(B))

def get_all_d_AB(G, process_idx, counter, disease_genes, matrix):
    print('starting process {}'.format(process_idx))
    while True:
        with counter.get_lock():
            idx = counter.value
            if idx >= len(disease_genes):
                break
            counter.value += 1
        print('{} out of {}'.format(idx+1, len(disease_genes)))
        gene_set_A = disease_genes[idx]
        for idx2 in range(idx+1, len(disease_genes)):
            gene_set_B = disease_genes[idx2]
            matrix[idx][idx2-idx-1] = separation.calc_set_pair_distances(G, gene_set_A, gene_set_B)
    print('ending process {}'.format(process_idx))

def plot_overlapping(G, all_genes_in_network):
    no_overlapping = list()
    subset = list()
    disease_genes = [separation.get_disease_genes(path, all_genes_in_network) for path in get_all_disease_files()]
    a = time()
    d_A = [separation.calc_single_set_distance(G, genes_set) for genes_set in disease_genes]
    processes = list()
    d_AB = [Array('d', [0] * (len(disease_genes)-idx-1)) for idx in range(len(disease_genes)-1)]
    b = time()-a
    print('took {} s to get all distances'.format(int(b)))
    a = time()
    NB_STEPS = 60
    NB_PROCESSES = 4
    counter = Value('i', 0)
    for k in range(NB_PROCESSES):
        processes.append(Process(target=get_all_d_AB, args=(G, k, counter, disease_genes, d_AB)))
        processes[-1].start()
    for process in processes:
        process.join()
    for idx, gene_set_A in enumerate(disease_genes):
        gene_set_A = disease_genes[idx]
        for idx2 in range(idx+1, len(disease_genes)):
            gene_set_B = disease_genes[idx2]
            J = J_score(gene_set_A, gene_set_B)
            C = C_score(gene_set_A, gene_set_B)
            sep = d_AB[idx][idx2-idx-1] - (d_A[idx]+d_A[idx2])/2
            # if no overlapping is observed in genes
            if J == C == 0:
                no_overlapping.append(sep)
            # if one set is fully included in the other
            elif J <= C == 1:
                subset.append(sep)
    b = time()-a
    print('took {} s to compute separation values'.format(int(b)))
    plt.subplot(121)
    plt.title('Complete subset')
    plt.hist(subset, 30, facecolor='g')
    plt.plot([0, 0], [0, 1e3], 'k-.')
    plt.xlabel(r'separation $s_{AB}$')
    plt.ylabel('Number of disease pairs')
    plt.annotate('{} pairs\n({} with'.format(len(subset), len([s_AB for s_AB in subset if s_AB < 0])) + (r' $s_{AB}$ < 0)'), xy=(-2, 100), xytext=(-2, 100))
    plt.yscale('log')
    plt.subplot(122)
    plt.title('No overlapping')
    plt.hist(no_overlapping, 30, facecolor='r')
    plt.plot([0, 0], [0, 1e4], 'k-.')
    plt.xlabel(r'separation $s_{AB}$')
    plt.ylabel('Number of disease pairs')
    plt.annotate('{} pairs\n({} with'.format(len(no_overlapping), len([s_AB for s_AB in no_overlapping if s_AB < 0])) + (r' $s_{AB}$ < 0)'), xy=(.5, 1100), xytext=(.5, 1100))
    plt.yscale('log')
    plt.show()

if __name__ == '__main__':
    a = time()
    opts = get_program_arguments()
    interactome_path = separation.INTERACTOME_DEFAULT_PATH
    G, all_genes_in_network = separation.load_network(interactome_path)
    if opts.plot_name == 'zscore_vs_rel_size':
        plot_zscore_vs_relative_size(G, all_genes_in_network, get_all_disease_files(),
                                     nb_simulations_per_disease=100000,
                                     show_non_significant=True,
                                     show_linear_regression=True,
                                     path=opts.output_file)
    elif opts.plot_name == 'overlapping':
        plot_overlapping(G, all_genes_in_network)
    else:
        print('Error: unknown plot name')
    print('Program took {} seconds to execute'.format(int(time()-a)))

