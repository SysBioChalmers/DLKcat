#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import os
import sys
import math
import model
import torch
import requests
import pickle
import numpy as np
from rdkit import Chem
from collections import defaultdict


fingerprint_dict = model.load_pickle('../../Data/input/fingerprint_dict.pickle')
atom_dict = model.load_pickle('../../Data/input/atom_dict.pickle')
bond_dict = model.load_pickle('../../Data/input/bond_dict.pickle')
edge_dict = model.load_pickle('../../Data/input/edge_dict.pickle')
word_dict = model.load_pickle('../../Data/input/sequence_dict.pickle')

def split_sequence(sequence, ngram):
    sequence = '-' + sequence + '='
    # print(sequence)
    # words = [word_dict[sequence[i:i+ngram]] for i in range(len(sequence)-ngram+1)]

    words = list()
    for i in range(len(sequence)-ngram+1) :
        try :
            words.append(word_dict[sequence[i:i+ngram]])
        except :
            word_dict[sequence[i:i+ngram]] = 0
            words.append(word_dict[sequence[i:i+ngram]])

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
    # atoms = list()
    # for a in atoms :
    #     try: 
    #         atoms.append(atom_dict[a])
    #     except :
    #         atom_dict[a] = 0
    #         atoms.append(atom_dict[a])

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
                # fingerprints.append(fingerprint_dict[fingerprint])
                # fingerprints.append(fingerprint_dict.get(fingerprint))
                try :
                    fingerprints.append(fingerprint_dict[fingerprint])
                except :
                    fingerprint_dict[fingerprint] = 0
                    fingerprints.append(fingerprint_dict[fingerprint])

            nodes = fingerprints

            """Also update each edge ID considering two nodes
            on its both sides."""
            _i_jedge_dict = defaultdict(lambda: [])
            for i, j_edge in i_jedge_dict.items():
                for j, edge in j_edge:
                    both_side = tuple(sorted((nodes[i], nodes[j])))
                    # edge = edge_dict[(both_side, edge)]
                    # edge = edge_dict.get((both_side, edge))
                    try :
                        edge = edge_dict[(both_side, edge)]
                    except :
                        edge_dict[(both_side, edge)] = 0
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

def load_tensor(file_name, dtype):
    return [dtype(d).to(device) for d in np.load(file_name + '.npy', allow_pickle=True)]

class Predictor(object):
    def __init__(self, model):
        self.model = model

    def predict(self, data):
        predicted_value = self.model.forward(data)

        return predicted_value

# One method to obtain SMILES by PubChem API using the website
def get_smiles(name):
    # smiles = redis_cli.get(name)
    # if smiles is None:
    try :
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/CanonicalSMILES/TXT' % name
        req = requests.get(url)
        if req.status_code != 200:
            smiles = None
        else:
            smiles = req.content.splitlines()[0].decode()
            # print(smiles)
        # redis_cli.set(name, smiles, ex=None)

        # print smiles
    except :
        smiles = None

    # name_smiles[name] = smiles
    return smiles

def main() :
    name = sys.argv[1:][0]
    print(name)
    # with open('./input.tsv', 'r') as infile :
    with open(name, 'r') as infile :
        lines = infile.readlines()

    fingerprint_dict = model.load_pickle('../../Data/input/fingerprint_dict.pickle')
    atom_dict = model.load_pickle('../../Data/input/atom_dict.pickle')
    bond_dict = model.load_pickle('../../Data/input/bond_dict.pickle')
    word_dict = model.load_pickle('../../Data/input/sequence_dict.pickle')
    n_fingerprint = len(fingerprint_dict)
    n_word = len(word_dict)

    radius=2
    ngram=3

    dim=10
    layer_gnn=3
    side=5
    window=11
    layer_cnn=3
    layer_output=3
    lr=1e-3
    lr_decay=0.5
    decay_interval=10
    weight_decay=1e-6
    iteration=100

    if torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    # torch.manual_seed(1234)
    Kcat_model = model.KcatPrediction(device, n_fingerprint, n_word, 2*dim, layer_gnn, window, layer_cnn, layer_output).to(device)
    Kcat_model.load_state_dict(torch.load('../../Results/output/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration30', map_location=device))
    # print(state_dict.keys())
    # model.eval()
    predictor = Predictor(Kcat_model)

    print('It\'s time to start the prediction!')
    print('-----------------------------------')

    i = 0

    with open('./output.tsv', 'w') as outfile :
        items = ['Substrate Name', 'Substrate SMILES', 'Protein Sequence', 'Kcat value (1/s)']
        outfile.write('\t'.join(items)+'\n')

        for line in lines[1:] :
            line_item = list()
            data = line.strip().split('\t')
            # i += 1
            # print('This is', i, '---------------------------------------')
            # print(data)

            name = data[0]
            smiles = data[1]
            sequence = data[2]
            if smiles and smiles != 'None' :
                smiles = data[1]
            else :
                smiles = get_smiles(name)
            # print(smiles)

            try :
                if "." not in smiles :
                    # i += 1
                    # print('This is',i)

                    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                    atoms = create_atoms(mol)
                    # print(atoms)
                    i_jbond_dict = create_ijbonddict(mol)
                    # print(i_jbond_dict)

                    fingerprints = extract_fingerprints(atoms, i_jbond_dict, radius)
                    # print(fingerprints)
                    # compounds.append(fingerprints)

                    adjacency = create_adjacency(mol)
                    # print(adjacency)
                    # adjacencies.append(adjacency)

                    words = split_sequence(sequence,ngram)
                    # print(words)
                    # proteins.append(words)

                    fingerprints = torch.LongTensor(fingerprints)
                    adjacency = torch.FloatTensor(adjacency)
                    words = torch.LongTensor(words)

                    inputs = [fingerprints, adjacency, words]

                    prediction = predictor.predict(inputs)
                    Kcat_log_value = prediction.item()
                    Kcat_value = '%.4f' %math.pow(2,Kcat_log_value)
                    # print(Kcat_value)
                    line_item = [name,smiles,sequence,Kcat_value]

                    outfile.write('\t'.join(line_item)+'\n')
            except :
                Kcat_value = 'None'
                line_item = [name,smiles,sequence,Kcat_value]
                outfile.write('\t'.join(line_item)+'\n')

    print('Prediction success!')


if __name__ == '__main__' :
    main()

