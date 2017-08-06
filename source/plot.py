#! /usr/bin/env python3

# std
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Rectangle
from glob import glob
from time import time
from multiprocessing import Process, Array, Value, Pool
from concurrent.futures import ProcessPoolExecutor
from math import isnan, log, ceil, exp
import shelve

import optparse

# local
import separation
import localization
import math_treatment as math
from constants import INTERACTOME_PATH, SHELVE_PATH, DISEASES_DIR, NEW_INTERACTOME_PATH, INTERACTOME_DEFAULT_PATH

rc('text', usetex=True)  # Activate LaTeX rendering

SIGNIFICANCE_THRESHOLD = 1.6

VERBOSE = False
FORCE = False

def get_zscore_key(nb_sims, path):
    return path + ';zscore;' + str(nb_sims)

def get_sep_key(subset_length):
    return INTERACTOME_PATH + ';separation;'+ str(subset_length)

def set_verbose(option=None, opt=None, value=None, parser=None):
    global VERBOSE
    VERBOSE = True

def set_force(option=None, opt=None, value=None, parser=None):
    global FORCE
    FORCE = True

def print_usage(option=None, opt=None, value=None, parser=None):
    """
    Displays the usage of the script
    """
    print("USAGE:\n\t./plot.py -h  for help")
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

    parser.add_option('-f', '--force',
                      help='force computation (do not use shelve-stored data',
                      dest='force',
                      action='callback',
                      callback=set_force)

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
    return sorted(glob(DISEASES_DIR + '*.txt'))

counter = Value('i', 0)
def get_zscore(G, all_genes_in_network, disease_file, nb_sims):
    gene_set = separation.get_disease_genes(disease_file, all_genes_in_network)
    mean, std, zscore = localization.get_random_comparison(G, gene_set, nb_sims)
    relative_size = localization.get_lcc_size(G, gene_set) / len(gene_set)
    with counter.get_lock():
        counter.value += 1
        print('{} -- done'.format(counter.value), end='\r')
    return relative_size, zscore

def get_relative_module_size_and_zscore(G, all_genes_in_network, disease_file_list,
        return_non_significant, significance_threshold, interactome_path=INTERACTOME_PATH):
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

    RETURNS:
    --------
        + relative_size_list: a list of relative module sizes
        + zscore_list: a list of associated z-scores
        + non_significant_counter: number of diseases which don't show
          significant z-scores
    """
    nb_diseases = len(disease_file_list)
    with shelve.open(SHELVE_PATH) as db:
        key = get_zscore_key(NB_SIMS, interactome_path)
        try:
            if FORCE:
                raise KeyError()
            relative_size_list, zscore_list = zip(*db[key].values())
        except KeyError:
            with ProcessPoolExecutor() as pool:
                params = [(G, all_genes_in_network, disease_file, NB_SIMS) for disease_file in disease_file_list]
                relative_size_list, zscore_list = zip(*list(pool.map(get_zscore, *zip(*params))))
            d = dict()
            for idx, name in enumerate(disease_file_list):
                d[name] = (relative_size_list[idx], zscore_list[idx])
            db[key] = d
    non_significant_counter = len([zscore for zscore in zscore_list if zscore < significance_threshold])
    if not return_non_significant:
        i = 0
        while i < len(relative_size_list):
            if zscore_list[i] < significance_threshold:
                del zscore_list[i]
                del relative_size_list[i]
            else:
                i += 1
    if VERBOSE:
        print('\n{} non significant z-scores'.format(non_significant_counter))
    return relative_size_list, zscore_list, non_significant_counter

def plot_linear_regression(x, y, xlabel='x', ylabel='y', color='b'):
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
    x_min, x_max = plt.xlim()
    regression_line = plt.plot([x_min, x_max], [m*x_min+p, m*x_max+p], color + '-', label=label)[0]
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
        show_non_significant=True, show_linear_regression=True):
    """
    make a plot of relative module size on x axis vs z-score of module size on y axis

    PARAMETERS:
    -----------
        + G: the interactome network
        + all_genes_in_network: a set of all the genes in G
        + disease_file_list: a list of paths to all desired diseases to plot
        + show_non_significant: a boolean telling whether non-significant module
          sizes should be displayed with the other (with a mark) or simply hidden

    RETURNS:
    --------
        + None
    """
    relative_size_list, zscore_list, non_significant_counter = \
        get_relative_module_size_and_zscore(G, all_genes_in_network,
            disease_file_list, show_non_significant, SIGNIFICANCE_THRESHOLD)

    fig = plt.figure()
    # dot plot for each disease
    plt.plot(relative_size_list, zscore_list, 'ro')
    print('interactome: {}\tmean z-score: {:.3g}\t(min, max) == {}\tmean relative size: {:.3g}\tmean(zscore*rel_size) == {:.3g}' \
          .format(INTERACTOME_PATH, sum(zscore_list)/len(zscore_list), (min(zscore_list), max(zscore_list)), sum(relative_size_list)/len(relative_size_list),
                  sum([zscore_list[i]*relative_size_list[i] for i in range(len(zscore_list))])/len(zscore_list)))
    # plot significance threshold
    if show_non_significant:
        significance_line = plt.plot([0, 1], [SIGNIFICANCE_THRESHOLD]*2, 'b--', label=r'significance threshold: $z$-score = ' + str(SIGNIFICANCE_THRESHOLD))[0]
        legend = plt.legend(handles=[significance_line], loc='lower right')
        plt.gca().add_artist(legend)
    # plot properties
    plt.xlabel(r'relative size: $r = S/N_d$')
    plt.ylabel(r'$z$-score of module size $S$')
    plt.title(r'$z$-score' + ' versus relative module size\nwith {} simulations per disease' \
              .format(NB_SIMS))
    #plt.axis([0, 1, 0, 35])
    if show_linear_regression:
        plot_linear_regression(relative_size_list, zscore_list, xlabel=r'$r$', ylabel=r'$z$-score')
    plt.annotate('{} diseases in non-significant area'.format(non_significant_counter), xy=(.2, -.5), xytext=(.2, -.5), fontsize=14)
    return fig

##### relative size histograms

def plot_relative_size_histogram(G, all_genes_in_network, disease_file_list):
    relative_sizes = get_relative_module_size_and_zscore(G, all_genes_in_network, disease_file_list, True, SIGNIFICANCE_THRESHOLD)[0]
    fig = plt.figure()
    plt.hist(relative_sizes, 30, facecolor='lightblue')
    plt.xlabel('relative size')
    plt.ylabel('number of diseases')
    #plt.yscale('log')
    plt.xlim((0, 1))
    return fig

def plot_relative_sizes_histogram_comparison(disease_file_list):
    return plot_data_histogram_comparison(disease_file_list, r'relative size $S/N_d$',
        'Comparison of the relative size distribution of both interactomes', 0)

def plot_data_histogram_comparison(disease_file_list, xlabel, title, idx, is_log=False):
    data =list()
    for path in (INTERACTOME_DEFAULT_PATH, NEW_INTERACTOME_PATH):
        G, all_genes_in_network = separation.load_network(path)
        data.append(get_relative_module_size_and_zscore(G, all_genes_in_network, disease_file_list, True, SIGNIFICANCE_THRESHOLD, path)[idx])
    fig = plt.figure()
    plt.hist(data, 20, color=('lightblue', 'salmon'), label=('original interactome', 'newer interactome'))
    plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('number of diseases')
    if is_log:
        plt.yscale('log')
    return fig

def plot_zscores_histogram_comparison(disease_file_list):
    return plot_data_histogram_comparison(disease_file_list, r'$z$-score',
        r'Comparison of the $z$-score distribution of both interactomes', 1)

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

def calc_set_dist(args):
    return separation.calc_single_set_distance(*args)

def get_all_distances(G, disease_genes):
    a = time()
    with shelve.open(SHELVE_PATH) as db:
        key = get_sep_key(len(disease_genes))
        try:
            d_A, d_AB = db[key]
        except KeyError:
            with ProcessPoolExecutor() as pool:
                d_A = list(pool.map(calc_set_dist, zip([G]*len(disease_genes), disease_genes)))
            processes = list()
            d_AB = [Array('d', [0] * (len(disease_genes)-idx-1)) for idx in range(len(disease_genes)-1)]
            b = time()-a
            print('took {} s to get all distances'.format(int(b)))
            a = time()
            NB_PROCESSES = 4
            counter = Value('i', 0)
            for k in range(NB_PROCESSES):
                processes.append(Process(target=get_all_d_AB, args=(G, k, counter, disease_genes, d_AB)))
                processes[-1].start()
            for process in processes:
                process.join()
            d_A, d_AB = list(d_A), [list(row) for row in d_AB]
            db[key] = d_A, d_AB
    return d_A, d_AB

def get_all_separations(G, disease_genes):
    return get_scores_and_separations(G, disease_genes)[1]

def get_all_scores(G, disease_genes):
    return get_scores_and_separations(G, disease_genes)[0]

def get_scores_and_separations(G, disease_genes):
    no_overlapping = list()
    subset = list()
    partial_overlap = list()
    all_separations = list()
    Js = list()
    Cs = list()
    d_A, d_AB = get_all_distances(G, disease_genes)
    for idx, gene_set_A in enumerate(disease_genes):
        gene_set_A = disease_genes[idx]
        for idx2 in range(idx+1, len(disease_genes)):
            gene_set_B = disease_genes[idx2]
            J = J_score(gene_set_A, gene_set_B)
            C = C_score(gene_set_A, gene_set_B)
            Js.append(J)
            Cs.append(C)
            sep = d_AB[idx][idx2-idx-1] - (d_A[idx]+d_A[idx2])/2
            # if no overlapping is observed in genes
            if J == C == 0:
                no_overlapping.append(sep)
            # if one set is fully included in the other
            elif J <= C == 1:
                subset.append(sep)
            # if sets partially overlap
            elif 0 < J < 1 and 0 < C < 1:
                partial_overlap.append(sep)
            all_separations.append(sep)
    return ((Js, Cs), (no_overlapping, subset, partial_overlap, all_separations))

def plot_overlapping(G, all_genes_in_network):
    disease_genes = [separation.get_disease_genes(path, all_genes_in_network) for path in get_all_disease_files()]
    a = time()
    no_overlapping, subset, partial_overlap, all_separations = get_all_separations(G, disease_genes)
    b = time()-a
    print('took {} s to compute separation values'.format(int(b)))
    fig = plt.figure()
    NB_BINS = 60
    bins = [-4] + [i/NB_BINS*6 - 4 for i in range(1, NB_BINS+1)]
    plt.subplot(221)
    plt.title(r'({\bf A}) No overlapping')
    plt.hist(no_overlapping, bins, color='tomato', label='mean: {:1.2f}'.format(np.mean(no_overlapping)))
    plt.plot([0, 0], [0, 1e4], 'k-.')
    plt.ylabel('number of disease pairs')
    plt.annotate('{} pairs\n({} with' \
                 .format(len(no_overlapping), len([s_AB for s_AB in no_overlapping if s_AB < 0])) + (r' $s_{AB} < 0$)'),
                 xy=(-3, 400), xytext=(-3, 400), fontsize=18)
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.subplot(222)
    plt.title(r'({\bf B}) Complete subset')
    plt.hist(subset, bins, facecolor='limegreen', label='mean: {:1.2f}'.format(np.mean(subset)))
    plt.plot([0, 0], [0, 1e3], 'k-.')
    plt.annotate('{} pairs\n({} with' \
                 .format(len(subset), len([s_AB for s_AB in subset if s_AB > 0])) + (r' $s_{AB} > 0$)'),
                 xy=(-3, 100), xytext=(-3, 100), fontsize=18)
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.subplot(223)
    plt.title(r'({\bf C}) Partial overlap')
    plt.hist(partial_overlap, bins, facecolor='paleturquoise', label='mean: {:1.2f}'.format(np.mean(partial_overlap)))
    plt.plot([0, 0], [0, 1e4], 'k-.')
    plt.xlabel(r'separation $s_{AB}$')
    plt.ylabel('number of disease pairs')
    plt.yscale('log')
    plt.legend(loc='upper right')
    plt.subplot(224)
    plt.title(r'({\bf D}) All disease pairs')
    plt.xlabel(r'separation $s_{AB}$')
    plt.hist(all_separations, bins, facecolor='orchid', label='mean: {:1.2f}'.format(np.mean(all_separations)))
    plt.plot([0, 0], [0, 1e5], 'k-.')
    plt.yscale('log')
    plt.legend(loc='upper right')
    return fig

def plot_J_C_scores(G, all_genes_in_network):
    disease_genes = [separation.get_disease_genes(path, all_genes_in_network) for path in get_all_disease_files()]
    Js, Cs = get_all_scores(G, disease_genes)
    fig = plt.figure()
    plt.subplot(221)
    plt.title('J-score vs C-score (linear scale)')
    plt.plot(Js, Cs, 'ro')
    plt.plot([0, 1], [0, 1], 'b--')
    plt.xlabel('J-score')
    plt.ylabel('C-score')
    plt.subplot(222)
    plt.title('J-score vs C-score (log-log scale)')
    plt.plot(Js, Cs, 'ro')
    plt.plot([1e-4, 1], [1e-4, 1], 'b--')
    plt.xlabel('J-score')
    plt.ylabel('C-score')
    plt.xscale('log')
    plt.yscale('log')
    plt.subplot(223)
    # legend
    handles = [Rectangle((0,0),1,1,color=c,ec="k") for c in ('tomato', 'steelblue')]
    labels = ['J-score', 'C-score']

    plt.title('J/C-score distribution (linear scale)')
    plt.legend(handles, labels)
    plt.hist((Js, Cs), 30, color=('tomato', 'steelblue'), label=('J-score', 'C-score'))
    plt.xlabel('score')
    plt.ylabel('number of pairs')
    plt.subplot(224)
    plt.title('J/C-score distribution (log scale)')
    plt.legend(handles, labels)
    plt.hist((Js, Cs), 30, color=('tomato', 'steelblue'))
    plt.xlabel('score')
    plt.ylabel('number of pairs')
    plt.yscale('log')
    return fig

def plot_degree_distribution(G):
    fig = plt.figure()
    plt.title('Degree distribution')
    plot_one_degree_distribution(G)
    plt.xlabel(r'degree ($k$)')
    plt.ylabel(r'probability ($P(k)$)')
    plt.xscale('log')
    plt.yscale('log')
    return fig

def plot_one_degree_distribution(G, color='r', label=''):
    degrees_list = list(zip(*G.degree_iter()))[1]
    degrees_count = dict()
    for degree in degrees_list:
        if degree not in degrees_count:
            degrees_count[degree] = 1
        else:
            degrees_count[degree] += 1
    degrees = list()
    frequencies = list()
    accumulated = list()
    accumulated_value = 0
    for degree in sorted(degrees_count.keys()):
        degrees.append(degree)
        frequencies.append(degrees_count[degree] / len(degrees_list))
        accumulated_value += frequencies[-1]
        accumulated.append(accumulated_value)
    mean_degree = round(sum([deg*freq for deg, freq in zip(degrees, frequencies)]))
    if label:
        label += ' --- '
    plt.plot(degrees, frequencies, color + 'o', markeredgecolor=color, markersize=3, label=label + 'mean degree = {}'.format(mean_degree))
    #plt.plot(degrees, accumulated, color + '-', label='Accumulated')
    return degrees, frequencies, accumulated

def get_linear_regression_coefficients_log(x, y):
    if 0 in x:
        del y[x.index(0)]
        del x[x.index(0)]
    x = [log(v, 10) for v in x]
    y = [log(v, 10) for v in y]
    return math.linear_regression(x, y)

def get_power_law_coefficient(x, y):
    return -get_linear_regression_coefficients_log(x, y)[0]

def plot_degree_distribution_comparison():
    fig = plt.figure()
    gathered_data = list()
    plt.subplot(121)
    colors = ('r', 'b')
    labels = ('newer interactome', 'original interactome')
    plt.title(r'({\bf A}) comparison of degree distribution between' + '\noriginal and newer interactomes')
    for idx, path in enumerate((NEW_INTERACTOME_PATH, INTERACTOME_DEFAULT_PATH)):
        G = separation.load_network(path)[0]
        degrees, frequencies, accumulated = plot_one_degree_distribution(G, color=colors[idx], label=labels[idx])
        multiplier = 1
        #multiplier = G.number_of_nodes()
        gathered_data.append((list(degrees), frequencies, [v*multiplier for v in accumulated]))
        power_law_coeff, offset = get_linear_regression_coefficients_log(degrees, frequencies)
        power_law_coeff = -power_law_coeff
        x_min, x_max = plt.xlim()
        if x_min == 0:
            x_min = 1
        NB_DOTS = 10000
        x_list = [x_min + (x_max-x_min)*i/NB_DOTS for i in range(NB_DOTS+1)]
        plt.plot(x_list, [pow(x, -power_law_coeff)*exp(offset) for x in x_list], colors[idx] + '--',
                 label='regression: $P(k) = ' + (str(round(exp(offset), 3))) + ' \\times k^{-' + (str(round(power_law_coeff, 3))) + '}$')
    plt.legend(prop={'size': 16})
    plt.xlabel(r'degree ($k$)')
    plt.ylabel(r'probability ($P(k)$)')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim((1e-5, 1))

    plt.subplot(122)
    y_values = list()
    for idx, (degrees, _, accumulated) in enumerate(gathered_data):
        plt.plot(degrees, accumulated, colors[idx] + '-', label=labels[idx])
        closest_to_100_idx = degrees.index(100 + min([abs(deg-100) for deg in degrees]))
        #print('closest is at index {} and is {}'.format(closest_to_100_idx, degrees[closest_to_100_idx]))
        y_values.append(accumulated[closest_to_100_idx])
    for idx, y in enumerate(y_values):
        plt.plot([0, 1e4], [y, y], colors[idx] + ':')
    plt.plot([100, 100], [0, 1], 'k:')
    plt.xscale('log')
    plt.title(r'({\bf B}) Comparison of cumulative degree frequencies distribution between' + '\noriginal and newer interactome')
    plt.xlabel('degree ($k$)')
    plt.ylabel('Proportion of nodes having degree $\leq k$')
    plt.ylim((0, 1))
    plt.yticks(np.arange(0, 1.1, .1))
    plt.legend(loc='lower right', prop={'size': 22})
    return fig

if __name__ == '__main__':
    a = time()
    opts = get_program_arguments()
    interactome_path = INTERACTOME_PATH
    G, all_genes_in_network = separation.load_network(interactome_path)
    NB_SIMS = 100000
    if opts.plot_name == 'zscore_vs_rel_size':
        fig = plot_zscore_vs_relative_size(G, all_genes_in_network, get_all_disease_files(),
                show_non_significant=True,
                show_linear_regression=True)
    elif opts.plot_name == 'relsize_histogram':
        fig = plot_relative_size_histogram(G, all_genes_in_network, get_all_disease_files())
    elif opts.plot_name == 'relsize_histogram_comparison':
        fig = plot_relative_sizes_histogram_comparison(get_all_disease_files())
    elif opts.plot_name == 'zscore_histogram_comparison':
        fig = plot_zscores_histogram_comparison(get_all_disease_files())
    elif opts.plot_name == 'overlapping':
        fig = plot_overlapping(G, all_genes_in_network)
    elif opts.plot_name == 'scores':
        fig = plot_J_C_scores(G, all_genes_in_network)
    elif opts.plot_name == 'degree':
        fig = plot_degree_distribution(G)
    elif opts.plot_name == 'degrees_comparison':
        fig = plot_degree_distribution_comparison()
    else:
        print('Error: unknown plot name')
        exit()
    print('Program took {} seconds to execute'.format(int(time()-a)))
    display_plot(opts.output_file, fig)

