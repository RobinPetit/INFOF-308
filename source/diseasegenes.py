#!/usr/bin/python3

from sys import argv
from os.path import isfile
from os.path import isdir

def print_usage():
    """
    print usage of program to help user
    """
    print("The given script is used to generate individual disease genes files " + \
          "named after the diseases in a given directory (which must exist!)\n" +
          "\tRun `{0} <disease genes tsv file> <directory to place the files in>`\n" + \
          "\tExample: {0} ../data/disease_genes.tsv ../data/diseases/"
          . format(argv[0]))

if __name__ == '__main__':
    if len(argv) != 3 or not isfile(argv[1]) or not isdir(argv[2]):
        print_usage()
        print(len(argv) != 3, not isfile(argv[1]), not isdir(argv[2]))
        exit()
    disease_genes_path = argv[1]
    with open(disease_genes_path) as disease_genes_file:
        for line in disease_genes_file.readlines():
            # ignore comments (typically first lines of tsv files)
            if line.strip()[0] == '#':
                continue
            # get each column of the row, supposed to be as follows:
            # disease | number_of_all_genes | number_of_OMIM_genes | number_of_GWAS_genes | OMIM_genes | GWAS_genes
            columns = line.strip().split('\t')
            # remove ',' character which doesn't fit file names
            disease_name = columns[0].replace(',', '')
            # if there exists OMIM genes for the present disease
            if columns[2] != '0':
                # then retrieve them all
                OMIM = columns[4].strip().split(';')
            # otherwise fake them by an empty list
            else:
                OMIM = list()

            # same for GWAS genes
            if columns[3] != '0':
                GWAS = columns[-1].strip().split(';')
            else:
                GWAS = list()
            # assemble these two
            complete_disease_genes_list = list(set(OMIM) | set(GWAS))
            current_disease_path = argv[2] + '/' + disease_name + '.txt'
            # open in write mode to avoid conflicts in files
            current_disease_file = open(current_disease_path, 'w')
            # write all the genes, oen per line, in the right file
            current_disease_file.write('\n'.join(complete_disease_genes_list))
