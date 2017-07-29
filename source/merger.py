from sys import argv

import csv

from separation import load_network

def save_interactions(interactions):
    interactome_file = csv.writer(open('merged_interactome.tsv', 'w'), delimiter='\t')
    interactome_file.writerows(interactions)

def main():
    if len(argv) != 3:
        print('give 2 interactomes to merge as parameters')
        return
    G1 = load_network(argv[1])[0]
    G2 = load_network(argv[2])[0]
    G1_interactions = set([frozenset(edge) for edge in G1.edges()])
    G2_interactions = set([frozenset(edge) for edge in G2.edges()])
    print('G1 has {} nodes and {} edges'.format(G1.number_of_nodes(), len(G1_interactions)))
    print('G2 has {} nodes and {} edges'.format(G2.number_of_nodes(), len(G2_interactions)))
    print('|G1 & G2| == {} and |G1 | G2| == {}'.format(len(G1_interactions & G2_interactions), len(G1_interactions | G2_interactions)))
    all_interactions = G1_interactions | G2_interactions
    save_interactions(all_interactions)

if __name__ == '__main__':
    main()
