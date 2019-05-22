#!/usr/bin/env python

# this file is part of the github repository: https://github.com/nwhoppe/rosalind
# author: nwhoppe
# created: 5/21/19

import argparse


def parse_pdb_id_from_blast_file(blast_txt_file):
    with open(blast_txt_file, 'r') as f:
        significant_hits_bool = False
        line = f.readline()
        while significant_hits_bool is False:
            if line.startswith('Sequences producing significant alignments:'):
                significant_hits_bool = True
            f.readline()
        # pull hits
        


def fetch_pdb_from_blast(blast_txt_file):
    pdb_id_list = parse_pdb_id_from_blast_file(blast_txt_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""takes output from blast alignment text as input and fetches pdb 
    structures using biopython to current directory """)
    required = parser.add_argument_group('required')
    required.add_argument('-i', '--input', required=True)
    args = parser.parse_args()
    fetch_pdb_from_blast(args.input)
