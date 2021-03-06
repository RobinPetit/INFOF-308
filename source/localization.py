#! /usr/bin/env python3

# -----------------------------------------------------------------------
#
# localization.py
#
# by Joerg Menche
# Last Modified: 2014-12-06
#
# This code determines the network-based distance and sepration for
# two given sets of nodes on given network as described in
#
# Uncovering Disease-Disease Relationships Through The Human
# Interactome
#
# by Joerg Menche, Amitabh Sharma, Maksim Kitsak, Susan Dina
#    Ghiassian, Marc Vidal, Joseph Loscalzo & Albert-Laszlo Barabasi
#
#
# -----------------------------------------------------------------------
#
#
# This program will calculate the size of the largest connected
# component S and mean shortest distance <d_s> for a given gene
# set. It will also compute the expected lcc size for the same number
# of randomly distributed genes.
#
# Check print_usage for further information
# -----------------------------------------------------------------------

import networkx as nx
import random
import numpy as np
import optparse
import sys

import separation as tools


#from lcc_size import LccCounter
#from lcc_size.LccCounter import LccCounter as LccCounterClass

INTERACTOME_DEFAULT_PATH = tools.INTERACTOME_DEFAULT_PATH
DEBUG = tools.DEBUG


# =============================================================================
#
#           S T A R T   D E F I N I T I O N S
#
# =============================================================================

def print_usage(option, opt, value, parser):
    usage_message = """
# ----------------------------------------------------------------------

This program will calculate the network-based localization for a given
gene set

* Required input:

  one files containing a gene set. The file must be in form of a
  table, one gene per line. If the table contains several columns,
  they must be tab-separated, only the first column will be used. See
  the two files MS.txt and PD.txt for valid examples (the contain
  genes for multiple sclerosis and peroxisomal disorders).

* Optional input:

  - file containing an interaction network. If now file is given, the
    default network \"{0}\" will be used instead. The file
    must contain an edgelist provided as a tab-separated table. The
    first two columns of the table will be interpreted as an
    interaction gene1 <==> gene2

  - filename for the output. If none is given,
    \"localiztion_results.txt\" will be used

  - the number or random simulations can be chosen. Default is 1000,
    which should run fast even for large gene sets and typically gives
    good result.


Here's an example that should work, provided the files are in the same
directory as this python script:

./localization.py -n {0} -g PD.txt -o output.txt

# ----------------------------------------------------------------------
    """.format(INTERACTOME_DEFAULT_PATH)

    print(usage_message)
    exit()

# =================================================================================
def get_density(G):
    nb_genes = G.number_of_nodes()
    nb_edges = G.number_of_edges()
    return nb_edges / LccCounter.X(nb_genes)

# =================================================================================
def get_lcc_size(G,seed_nodes):
    """
    return the lcc size
    """

    # getting subgraph that only consists of the black_nodes
    g = nx.subgraph(G,seed_nodes)

    if g.number_of_nodes() != 0:
        # get all components
        components = nx.connected_components(g)
        # return size of biggest connected component
        return max(list(map(len, components)))
    else:
        return 0

def get_random_comparison(G, gene_set, sims):
    return get_random_comparison_simulation(G, gene_set, sims)
    #return get_random_comparison_probability(G, gene_set)

def get_random_comparison_probability(G, gene_set):
    interactome_density = get_density(G)
    print('interactome density: {}'.format(interactome_density))
    # a random graph in the interactome having `gene_set` vertices should have
    # interactome_density * X(|gene_set|) edges
    nb_random_vertices = len(gene_set & set(G.nodes()))
    nb_random_edges = round(interactome_density * LccCounter.X(nb_random_vertices))
    print('random for (n, m) == {}'.format((nb_random_vertices, nb_random_edges)))
    l_mean = LccCounterClass.get_expectation_lcc_m(nb_random_vertices, nb_random_edges)
    l_std = LccCounterClass.get_std_lcc_m(nb_random_vertices, nb_random_edges)
    lcc_observed = get_lcc_size(G, gene_set)
    z_score = (1.*lcc_observed - l_mean)/l_std
    return l_mean, l_std, z_score

# =============================================================================
simulations_cache = dict()
def get_random_comparison_simulation(G,gene_set,sims):
    """
    gets the random expectation for the lcc size for a given gene set
    by drawing the same number of genes at random from the network

    PARAMETERS:
    -----------
        - G       : network
        - gene_set: dito
        - sims    : number of random simulations

    RETURNS:
    --------
        - a tuple containing (mean, std, zscore)
    """
    global simulations_cache
    # getting all genes in the network
    all_genes = G.nodes()

    number_of_seed_genes = len(gene_set & set(all_genes))
    if number_of_seed_genes not in simulations_cache:
        l_list = []

        # simulations with randomly distributed seed nodes
        for i in range(1,sims+1):
            # print out status
            if i % 50 == 0:
                if DEBUG:
                    sys.stdout.write("> random simulation [{} of {}]\r" \
                                     .format(i, sims))
                    sys.stdout.flush()

            # get random seeds
            rand_seeds = set(random.sample(all_genes,number_of_seed_genes))

            # get rand lcc
            lcc = get_lcc_size(G,rand_seeds)
            l_list.append(lcc)

        # get the lcc z-score:
        l_mean = np.mean(l_list)
        l_std  = np.std(l_list)
        simulations_cache[number_of_seed_genes] = (l_mean, l_std)
    else:
        l_mean, l_std = simulations_cache[number_of_seed_genes]

    # get the actual value
    lcc_observed = get_lcc_size(G,gene_set)


    if l_std == 0:
        z_score = 'not available'
    else:
        z_score = (1.*lcc_observed - l_mean)/l_std

    return (l_mean, l_std, z_score)

# =============================================================================
def get_program_arguments():
    """
    parses the arguments given in command line and return them all
    """
    parser = optparse.OptionParser()

    parser.add_option('-u', '--usage',
                      help    ='print more info on how to use this script',
                      action  ="callback",
                      callback=print_usage)

    parser.add_option('-n',
                      help    ='file containing the network edgelist [{}]' \
                               .format(INTERACTOME_DEFAULT_PATH),
                      dest    ='network_file',
                      default =INTERACTOME_DEFAULT_PATH,
                      type    = "string")

    parser.add_option('-g',
                      help    ='file containing gene set',
                      dest    ='gene_file',
                      default ='none',
                      type    = "string")

    parser.add_option('-s',
                      help    ='number of random simulations [1000]',
                      dest    ='sims',
                      default ='1000',
                      type    = "int")

    parser.add_option('-o',
                      help    ='file for results [separation_results.txt]',
                      dest    ='results_file',
                      default ='localization_results.txt',
                      type    = "string")


    (opts, args) = parser.parse_args()
    return opts

def check_arguments():
    """
    checks if the given arguments are correct, and exit if they are not
    """
    if gene_file == 'none':
        error_message = \
        "\tERROR: you must specify an input file with a gene set, for example:\n" + \
        "\t./localization.py -g MS.txt\n\n"  + \
        "\tFor more information, type\n" + \
        "\t./localization.py --usage"

        print(error_message)
        exit(0)


# =============================================================================
#
#           E N D    O F    D E F I N I T I O N S
#
# =============================================================================


if __name__ == '__main__':

    # "Hey Ho, Let's go!" -- The Ramones (1976)

    # --------------------------------------------------------
    #
    # PARSING THE COMMAND LINE
    #
    # --------------------------------------------------------

    opts = get_program_arguments()

    network_file = opts.network_file
    gene_file    = opts.gene_file
    results_file = opts.results_file
    sims         = opts.sims

    check_arguments()

    if network_file == INTERACTOME_DEFAULT_PATH and DEBUG:
        print('> default network from "{}" will be used' \
              .format(INTERACTOME_DEFAULT_PATH))

    # --------------------------------------------------------
    #
    # LOADING NETWORK and DISEASE GENES
    #
    # --------------------------------------------------------

    G, all_genes_in_network = tools.load_network(network_file)
    gene_set = tools.get_disease_genes(gene_file, all_genes_in_network)

    # --------------------------------------------------------
    #
    # CALCULATE NETWORK QUANTITIES
    #
    # --------------------------------------------------------

    # get lcc size S
    lcc = get_lcc_size(G,gene_set)
    if DEBUG:
        print("\n> lcc size = {}".format(lcc))

    # get mean shortest distance
    d_s = tools.calc_single_set_distance(G,gene_set)
    if DEBUG:
        print("> mean shortest distance = {}".format(d_s))

    results_message = ("> gene set from \"{}\": {} genes\n" + \
                       "> lcc size   S = {}\n" + \
                       "> diameter d_s = {}") \
                       .format(gene_file, len(gene_set), lcc, d_s)

    # --------------------------------------------------------
    #
    # CALCULATE RANDOM COMPARISON
    #
    # --------------------------------------------------------

    mean, std, z_score = get_random_comparison(G,gene_set,sims)
    results_message += ("> Random expecation:\n" + \
                      "> lcc [rand] = {}\n" + \
                      "> => z-score of observed lcc = {}") \
                      .format(mean, z_score)

    if DEBUG:
        print(results_message)

    fp = open(results_file,'w')
    fp.write(results_message)
    fp.close()

    if DEBUG:
        print("> results have been saved to {}".format(results_file))

