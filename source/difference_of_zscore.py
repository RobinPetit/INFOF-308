from plot import *
from constants import *

import shelve
import matplotlib.pyplot as plt

ORIGINAL_INTERACTOME_RESULTS_FILE = '../data/DataS3_localization.tsv'

ABS = False

def get_file_name(path):
    return path[path.rfind('/')+1:path.rfind('.')]

def get_zscores():
    with shelve.open(SHELVE_PATH) as db:
        zscores = db[get_zscore_key(100000, INTERACTOME_PATH)]
        keys = list(zscores.keys())
        for path in keys:
            zscores[get_file_name(path)] = zscores[path]
            del zscores[path]
    return zscores

if __name__ == '__main__':
    zscores = get_zscores()
    differences = dict()
    with open(ORIGINAL_INTERACTOME_RESULTS_FILE, 'r') as results_file:
        for line in results_file:
            if line.startswith('#'):
                continue
            columns = line.split('\t')
            disease_name = columns[0].replace(',', '')
            diff = zscores[disease_name][1] - float(columns[3])
            if ABS:
                diff = abs(diff)
            differences[disease_name] = diff
    plt.hist(list(differences.values()), 50, color='tomato')
    plt.title('Distribution of $z$-score difference distribution between\nprovided results and computed results')
    plt.xlabel(r'$\left|z^{\mathrm {provided}} - z^{\mathrm {computed}}\right|$')
    plt.ylabel('Number of diseases')
    plt.yscale('log')
    plt.ylim(.9, plt.ylim()[1])
    plt.show()
