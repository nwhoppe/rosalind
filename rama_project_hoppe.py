#!/usr/bin/env python

# this file is part of the rosalind github repository for UCSF's IPQB program:
# https://github.com/nwhoppe/rosalind
# author: nwhoppe
# created: 8/30/18

import argparse

import collections

import sys

import math


def cartesian_coords_from_pdb(pdb_file, backbone_only=False, chain_to_take="all"):
    """:return nested dictionary with chain as first key, residue position as second key,
    and values as lists of x, y, and z coordinates from input pdb file"""

    cartesian_coords_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                split_line = line.split()
                chain = split_line[4]
                atom_name = split_line[2]
                position = int(split_line[5])
                x, y, z = map(float, split_line[6:9])
                if chain_to_take == 'all' or chain == chain_to_take:
                    if backbone_only and atom_name in {'N', 'CA', 'C'}:
                        cartesian_coords_dict[chain][position].append((x, y, z))
                    elif not backbone_only:
                        cartesian_coords_dict[chain][position].append((x, y, z))
    return cartesian_coords_dict


def undefined_positions_from_list(position_list):
    """:return set of positions with undefined torsion angles
    undefined positions are i - 1, i, and i + 1 for position i missing from position list
    first and last position are also undefined"""
    missing_positions = set(xrange(position_list[0], position_list[-1] + 1)).difference(position_list)
    undefined_positions = {position_list[0]}
    for missing_position in missing_positions:
            undefined_positions.update([missing_position - 1, missing_position, missing_position + 1])
    undefined_positions.add(position_list[-1])
    return undefined_positions


def add_vectors(vector1, vector2):
    return tuple(v1 + v2 for v1, v2 in zip(vector1, vector2))


def subtract_vectors(vector1, vector2):
    return tuple(v1 - v2 for v1, v2 in zip(vector1, vector2))


def normalize_vector(vector):
    magnitude = float(sum([v * v for v in vector])) ** (1.0/2)
    return tuple(v / magnitude for v in vector)


def dot_product(vector1, vector2):
    return sum(v1 * v2 for v1, v2 in zip(vector1, vector2))


def cross_product(vector1, vector2):
    x = vector1[1] * vector2[2] - vector1[2] * vector2[1]
    y = vector1[2] * vector2[0] - vector1[0] * vector2[2]
    z = vector1[0] * vector2[1] - vector1[1] * vector2[0]
    return tuple([x, y, z])


def torsion_from_cartesian(cartesian_tup1, cartesian_tup2, cartesian_tup3, cartesian_tup4):
    """:return torsion angle in degrees calculated from 4 tuples of 3D cartesian coordinates"""
    bond1 = subtract_vectors(cartesian_tup2, cartesian_tup1)
    bond2 = subtract_vectors(cartesian_tup3, cartesian_tup2)
    bond3 = subtract_vectors(cartesian_tup4, cartesian_tup3)

    # normal vectors to the planes defined by (bond1 and bond2) and (bond2 and bond3)
    normal1 = normalize_vector(cross_product(bond1, bond2))
    normal2 = normalize_vector(cross_product(bond2, bond3))

    # normal1, normalized bond2, and normal1 cross normalized bond 2 form an right handed orthonormal basis
    basis_i = normal1
    basis_j = normalize_vector(bond2)
    # dihedral angles use left handed basis
    basis_k = cross_product(basis_j, basis_i)

    # torsion angle is angle between normal2 and basis_i in new orthonormal basis
    component_i = dot_product(normal2, basis_i)
    component_k = dot_product(normal2, basis_k)
    torsion_angle = math.atan2(component_k, component_i)
    return math.degrees(torsion_angle)


def ramachandran_plot(phi_list, psi_list):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.scatter(phi_list, psi_list, )
    plt.tick_params(axis='both', direction='in')
    plt.grid(linestyle='dashed')
    plt.xlim((-180, 180))
    plt.ylim((-180, 180))
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\psi$')
    plt.savefig('rama_plot.pdf')


def main(pdb_files, pdb_chain):
    phi_list = []
    psi_list = []

    for pdb_file in pdb_files:
        cartesian_coord_dict = cartesian_coords_from_pdb(pdb_file, backbone_only=True, chain_to_take=pdb_chain)
        for chain in cartesian_coord_dict.keys():
            position_list = cartesian_coord_dict[chain].keys()
            undefined_positions = undefined_positions_from_list(position_list)
            sys.stderr.write('Undefined torsion angles for pdb {0} chain {1} positions: {2}\n'.format(
                pdb_file, chain, ' '.join(map(str, undefined_positions))
            ))
            defined_positions = set(position_list) - undefined_positions
            for position in defined_positions:
                try:
                    c_minus1 = cartesian_coord_dict[chain][position - 1][2]
                    n, ca, c = cartesian_coord_dict[chain][position]
                    n_plus1 = cartesian_coord_dict[chain][position + 1][0]
                except IndexError:
                    sys.stderr.write('Missing atom for pdb {0} chain {1} position {2}\n'.format(
                        pdb_file, chain, position
                    ))
                    continue
                phi = torsion_from_cartesian(c_minus1, n, ca, c)
                psi = torsion_from_cartesian(n, ca, c, n_plus1)
                phi_list.append(phi)
                psi_list.append(psi)
    ramachandran_plot(phi_list, psi_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""script to make Ramachandran plots for input pdb files""")
    parser.add_argument("-c", "--chain", default='all',
                        help='chain from pdb file to plot. default is all chains.')
    required = parser.add_argument_group('required')
    required.add_argument("-p", "--pdbs", nargs='*', required=True)
    args = parser.parse_args()
    main(args.pdbs, args.chain)
