#!/usr/bin/python3

import matplotlib.pyplot as plt

from multiprocessing import Array, Process, Value
from glob import glob
from time import time

import separation
import localization

import math_treatment as math

DISEASES_DIR = '../data/diseases/'
INTERACTOME_PATH = '../data/DataS1_interactome.tsv'
NB_SIMS = 100
STEP = 10

def compute_mean_std(process_id, counter, mean_list, std_list):
    print('process {} starting'.format(process_id))
    amount_tests = len(mean_list)
    while True:
        with counter.get_lock():
            if counter.value > amount_tests:
                break
            n = int(counter.value)
            counter.value += 1
        print('{} out of {}'.format(STEP*n, STEP*amount_tests))
        mean, std = localization.get_random_comparison(G, set(all_genes_as_list[:STEP*n]), NB_SIMS)[:2]
        mean_list[n-1] = mean
        std_list[n-1] = std
    print('process {} finishing'.format(process_id))

if __name__ == '__main__':
    a = time()
    G, all_genes_in_network = separation.load_network(INTERACTOME_PATH)
    all_genes_as_list = list(all_genes_in_network)
    max_value = G.number_of_nodes() // STEP
    NB_PROCESSES = 4
    counter = Value('i', 1)
    mean_list = Array('d', max_value)
    std_list = Array('d', max_value)
    processes = list()
    for k in range(NB_PROCESSES):
        processes.append(Process(target=compute_mean_std, args=(k+1, counter, mean_list, std_list)))
        processes[-1].start()
    for process in processes:
        process.join()

    x_list = [STEP*k for k in range(1, max_value+1)]
    plt.subplot(121)
    plt.title('mean of random module size\nvs size of gene subset ({} simulations)'.format(NB_SIMS))
    mean_plot = plt.plot(x_list, mean_list, 'r-', label=r'$\left\langle S^r\right\rangle$')
    plt.axis([0, STEP*max_value*1.1, 0, max(mean_list)*1.1])
    plt.xlabel('number of genes in subset')
    m, p = math.linear_regression(x_list, mean_list)
    regression_line = plt.plot([0, STEP*max_value*1.1], [p, STEP*max_value*1.1*m + p],
        label=r'$\langle S^r \rangle = $' + ('{:.2f}'.format(m)) + r'$\times$ n' + ('+' if p >= 0 else '-') + ('{:.2f}'.format(abs(p))))
    plt.legend(handles=[mean_plot[0], regression_line[0]], loc='upper left')
    plt.subplot(122)
    plt.title('standard deviation of random module size\nvs size of gene subset ({} simulations)'.format(NB_SIMS))
    std_plot = plt.plot(x_list, std_list, 'b-', label=r'$\sigma\left(S^r\right)$')
    plt.axis([0, STEP*max_value*1.1, 0, max(std_list)*1.1])
    plt.xlabel('number of genes in subset')
    plt.legend(handles=[std_plot[0]], loc='upper left')
    print('\ntook {} seconds to execute'.format(int(time()-a)))
    plt.show()
