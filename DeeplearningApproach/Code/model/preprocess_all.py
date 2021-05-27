#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-10-03

import math
import json
import pickle
import numpy as np
from collections import defaultdict
from rdkit import Chem


word_dict = defaultdict(lambda: len(word_dict))
atom_dict = defaultdict(lambda: len(atom_dict))
bond_dict = defaultdict(lambda: len(bond_dict))
fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
edge_dict = defaultdict(lambda: len(edge_dict))

proteins = list()
compounds = list()
adjacencies = list()
regression =list()

def split_sequence(sequence, ngram):
    sequence = '-' + sequence + '='
    # print(sequence)
    words = [word_dict[sequence[i:i+ngram]] for i in range(len(sequence)-ngram+1)]
    return np.array(words)
    # return word_dict

def create_atoms(mol):
    """Create a list of atom (e.g., hydrogen and oxygen) IDs
    considering the aromaticity."""
    # atom_dict = defaultdict(lambda: len(atom_dict))
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    # print(atoms)
    for a in mol.GetAromaticAtoms():
        i = a.GetIdx()
        atoms[i] = (atoms[i], 'aromatic')
    atoms = [atom_dict[a] for a in atoms]
    return np.array(atoms)

def create_ijbonddict(mol):
    """Create a dictionary, which each key is a node ID
    and each value is the tuples of its neighboring node
    and bond (e.g., single and double) IDs."""
    # bond_dict = defaultdict(lambda: len(bond_dict))
    i_jbond_dict = defaultdict(lambda: [])
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond = bond_dict[str(b.GetBondType())]
        i_jbond_dict[i].append((j, bond))
        i_jbond_dict[j].append((i, bond))
    return i_jbond_dict

def extract_fingerprints(atoms, i_jbond_dict, radius):
    """Extract the r-radius subgraphs (i.e., fingerprints)
    from a molecular graph using Weisfeiler-Lehman algorithm."""

    # fingerprint_dict = defaultdict(lambda: len(fingerprint_dict))
    # edge_dict = defaultdict(lambda: len(edge_dict))

    if (len(atoms) == 1) or (radius == 0):
        fingerprints = [fingerprint_dict[a] for a in atoms]

    else:
        nodes = atoms
        i_jedge_dict = i_jbond_dict

        for _ in range(radius):

            """Update each node ID considering its neighboring nodes and edges
            (i.e., r-radius subgraphs or fingerprints)."""
            fingerprints = []
            for i, j_edge in i_jedge_dict.items():
                neighbors = [(nodes[j], edge) for j, edge in j_edge]
                fingerprint = (nodes[i], tuple(sorted(neighbors)))
                fingerprints.append(fingerprint_dict[fingerprint])
            nodes = fingerprints

            """Also update each edge ID considering two nodes
            on its both sides."""
            _i_jedge_dict = defaultdict(lambda: [])
            for i, j_edge in i_jedge_dict.items():
                for j, edge in j_edge:
                    both_side = tuple(sorted((nodes[i], nodes[j])))
                    edge = edge_dict[(both_side, edge)]
                    _i_jedge_dict[i].append((j, edge))
            i_jedge_dict = _i_jedge_dict

    return np.array(fingerprints)

def create_adjacency(mol):
    adjacency = Chem.GetAdjacencyMatrix(mol)
    return np.array(adjacency)

def dump_dictionary(dictionary, filename):
    with open(filename, 'wb') as file:
        pickle.dump(dict(dictionary), file)

def main() :
    with open('../../Data/database/Kcat_combination_0918.json', 'r') as infile :
        Kcat_data = json.load(infile)

    # smiles_all = [data['Smiles'] for data in Kcat_data]

    # print(len(Kcat_data))

    # smiles = "CC1=NC=C(C(=C1O)CO)CO"
    # radius = 3 # The initial setup, I suppose it is 2, but not 2.
    radius = 2
    ngram = 3

    """Exclude data contains '.' in the SMILES format."""
    i = 0
    for data in Kcat_data :
        smiles = data['Smiles']
        sequence = data['Sequence']
        # print(smiles)
        Kcat = data['Value']
        if "." not in smiles and float(Kcat) > 0:
            # i += 1
            # print('This is',i)
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            atoms = create_atoms(mol)
            # print(atoms)
            i_jbond_dict = create_ijbonddict(mol)
            # print(i_jbond_dict)

            fingerprints = extract_fingerprints(atoms, i_jbond_dict, radius)
            # print(fingerprints)
            compounds.append(fingerprints)

            adjacency = create_adjacency(mol)
            adjacencies.append(adjacency)

            words = split_sequence(sequence,ngram)
            # print(words)
            proteins.append(words)

            # print(float(Kcat))

            regression.append(np.array([math.log2(float(Kcat))]))
            print(math.log2(float(Kcat)))

            # regression.append(np.array([math.log10(float(Kcat))]))
            # print(math.log10(float(Kcat)))

    np.save('../../Data/input/'+'compounds', compounds)
    np.save('../../Data/input/'+'adjacencies', adjacencies)
    np.save('../../Data/input/'+'regression', regression)
    np.save('../../Data/input/'+'proteins', proteins)

    dump_dictionary(fingerprint_dict, '../../Data/input/fingerprint_dict.pickle')
    dump_dictionary(atom_dict, '../../Data/input/atom_dict.pickle')
    dump_dictionary(bond_dict, '../../Data/input/bond_dict.pickle')
    dump_dictionary(edge_dict, '../../Data/input/edge_dict.pickle')
    dump_dictionary(word_dict, '../../Data/input/sequence_dict.pickle')


if __name__ == '__main__' :
    main()
