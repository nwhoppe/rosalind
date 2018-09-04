#!/usr/bin/env python

# this file is part of the github repository: https://github.com/nwhoppe/rosalind
# author: nwhoppe
# created: 9/3/18

import argparse
import sys


def intro_pattern_matching(input_txt):
    """
    From the input dna sequences, an adjacency dictionary is made corresponding the the underlying trie
    representing the sequences. For now, the trie is not explicitly made. Root node is 1. Unique dna bases at each
    position start a new child node. The edge is labeled with the nucleotide.

    Args:
        input_txt: text file containing one DNA sequence per line.
            each sequence is <= 100 nucleotides
            up to 100 sequences

    Returns:
         Adjacency dictionary corresponding to the trie representing the input sequences. Parent nodes are keys to
         the first level of the dictionary. Edge labels are the keys to the second level, and values at the second
         level are child nodes.

    Raises:
        IOError:
            1. more than 100 sequences
            2. more than 100 nucleotides in a sequence
            3. sequence is not dna (must only have A, C, G, T)
    """

    adjacency_dict = {1: {}}
    with open(input_txt, 'r') as f:
        dna_sequences = f.readlines()
    child_node = 2

    if len(dna_sequences) > 100:
        raise IOError("No more than 100 sequences can be given as input\n")

    for sequence in dna_sequences:
        sequence = sequence.rstrip().upper()
        if len(sequence) > 100:
            raise IOError("DNA sequences cannot be longer than 100 nucleotides\n")
        elif set(sequence) > set("ACGT"):
            raise IOError("DNA sequences must only contain canonical nucleotides A, C, G, T\n")
        else:
            parent_node = 1
            for index, nucleotide in enumerate(sequence):
                if parent_node in adjacency_dict.keys() and nucleotide not in adjacency_dict[parent_node].keys():
                    # add edge and child node to parent that has at least 1 child node
                    adjacency_dict[parent_node][nucleotide] = child_node
                    parent_node = child_node
                    child_node += 1
                elif parent_node not in adjacency_dict.keys():
                    # add edge and child node to leaf node
                    adjacency_dict[parent_node] = {nucleotide: child_node}
                    parent_node = child_node
                    child_node += 1
                else:
                    # parent, edge, and child already existed - set parent to existing child node
                    parent_node = adjacency_dict[parent_node][nucleotide]

    return adjacency_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""""")
    required = parser.add_argument_group('required')
    required.add_argument("-i", "--input_txt", required=True,
                          help="input txt with up to 100 dna sequences up to 100 nucleotides "
                               "in length. one sequence per line")
    args = parser.parse_args()
    output_dict = intro_pattern_matching(args.input_txt)
    for parent in output_dict.iterkeys():
        for edge_label, child in output_dict[parent].iteritems():
            sys.stdout.write("{0} {1} {2}\n".format(parent, child, edge_label))
