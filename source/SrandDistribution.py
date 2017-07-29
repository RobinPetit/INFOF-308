#!/usr/bin/python3

# std
from glob import glob
from time import time
import shelve
from multiprocessing import Value, Pool

from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

# local
import separation
import localization
from constants import *
import math_treatment as math

INTERACTOME_PATH = '../data/DataS1_interactome.tsv'
NB_SIMS = 1000
STEP = 10

def get_key():
    SEPARATOR = ';'
    return SEPARATOR.join(('srand', DISEASES_DIR, INTERACTOME_PATH, str(NB_SIMS), str(STEP)))

nb_processes = Value('i', 0)
def compute_mean_std_process(k):
    with nb_processes.get_lock():
        nb_processes.value += 1
        n = nb_processes.value
    ret = localization.get_random_comparison(G, set(all_genes_as_list[:STEP*k]), NB_SIMS)[:2]
    print('{} out of {}'.format(nb_processes.value, G.number_of_nodes() // STEP + 1), end='\r')
    return ret

def compute_mean_std():
    max_value = G.number_of_nodes() // STEP + 1
    NB_PROCESSES = 4
    results = list()
    with Pool() as pool:
        results = pool.map(compute_mean_std_process, range(max_value))
    mean_list = list()
    std_list = list()
    for idx, result in enumerate(results):
        mean, std = result
        mean_list.append(mean)
        std_list.append(std)
    return mean_list, std_list

def plot_srand_mean_std(means, stds, plot=False, same=True):
    if same:
        plot_mean_std_same_plot(means, stds, plot)
    else:
        plot_mean_std_horizontal(means, stds, plot)

def plot_mean_std_same_plot(means, stds, plot):
    x_list = [STEP*k for k in range(1, len(means)+1)]
    m, p = math.linear_regression(x_list, mean_list)

    host = host_subplot(111, axes_class=AA.Axes)

    right_axis = host.twinx()
    host.set_xlabel('number of genes in subset')
    host.set_ylabel('mean of random module size')
    right_axis.set_ylabel('std of random module size')

    mean_plot = host.plot(x_list, means)
    std_plot = right_axis.plot(x_list, stds)
    regression_line_plot = host.plot([0, STEP*len(means)], [p, STEP*len(means)*m + p],
        label=r'$\langle S^r \rangle = $' + ('{:.2f}'.format(m)) + r'$\times$n ' + ('+' if p >= 0 else '-') + (' {:.2f}'.format(abs(p))))
    host.legend(loc='upper left')
    host.plot([0, STEP*len(means)], [0, STEP*len(means)], label=r'$y=x$')

    host.axis['left'].label.set_color(mean_plot[0].get_color())
    right_axis.axis['right'].label.set_color(std_plot[0].get_color())
    host.set_title('standard deviation of random module size\nvs size of gene subset ({} simulations)'.format(NB_SIMS))
    plt.draw()
    if plot:
        plt.plot()

def plot_mean_std_horizontal(means, stds, plot):
    x_list = [STEP*k for k in range(1, len(mean_list)+1)]
    plt.subplot(121)
    plt.title('mean of random module size\nvs size of gene subset ({} simulations)'.format(NB_SIMS))
    mean_plot = plt.plot(x_list, mean_list, 'r-', label=r'$\left\langle S^r\right\rangle$')
    plt.axis([0, STEP*len(mean_list)*1.1, 0, max(mean_list)*1.1])
    plt.xlabel('number of genes in subset')
    m, p = math.linear_regression(x_list, mean_list)
    regression_line = plt.plot([0, STEP*len(mean_list)*1.1], [p, STEP*len(mean_list)*1.1*m + p],
        label=r'$\langle S^r \rangle = $' + ('{:.2f}'.format(m)) + r'$\times$ n' + ('+' if p >= 0 else '-') + ('{:.2f}'.format(abs(p))))
    plt.legend(handles=[mean_plot[0], regression_line[0]], loc='upper left')
    plt.subplot(122)
    plt.title('mean and standard deviation of random module size vs size of gene subset ({} simulations)'.format(NB_SIMS))
    std_plot = plt.plot(x_list, std_list, 'b-', label=r'$\sigma\left(S^r\right)$')
    plt.axis([0, STEP*len(mean_list)*1.1, 0, max(std_list)*1.1])
    plt.xlabel('number of genes in subset')
    plt.legend(handles=[std_plot[0]], loc='upper left')
    if plot:
        plt.plot()

if __name__ == '__main__':
    a = time()
    with shelve.open(SHELVE_PATH) as db:
        try:
            mean_list, std_list = db[get_key()]
        except KeyError:
            G, all_genes_in_network = separation.load_network(INTERACTOME_PATH)
            all_genes_as_list = list(all_genes_in_network)
            mean_list, std_list = compute_mean_std()
            db[get_key()] = (mean_list, std_list)
    plot_srand_mean_std(mean_list, std_list)
    print('\ntook {} seconds to execute'.format(int(time()-a)))
    plt.show()
