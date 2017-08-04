#!/usr/bin/python3

from time import time
from sys import argv
from os.path import exists
import networkx as nx

import separation

if __name__ == '__main__':
    if len(argv) != 2 or not exists(argv[1]):
        print('Usage: ./interactome_diameter.py <interactome tsv path>')
    G, genes = separation.load_network(argv[1])
    genes = list(sorted(genes))
    total_nb_of_computations = len(genes)*(len(genes)-1)//2
    sum_of_paths = 0
    counter = 0
    steps = 0
    a = time()
    for idx1 in range(len(genes)):
        for idx2 in range(idx1+1, len(genes)):
            steps += 1
            if steps % 2000 == 0:
                print('estimated remaining time: {:2.2f}s\t{:2.2f}%'.format((time()-a)/steps*(total_nb_of_computations-steps), 100*steps/total_nb_of_computations), end='\r')
            try:
                sum_of_paths += nx.shortest_path_length(G, genes[idx1], genes[idx2])
                counter += 1
            except nx.exception.NetworkXNoPath:
                # sum_of_paths += 0
                pass
    print('\ndiameter: {:.4g}'.format(sum_of_paths / total_nb_of_computations))
