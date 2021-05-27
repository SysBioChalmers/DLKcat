#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-10-01

import os
import math
import model
import torch
import json
import pickle
import numpy as np
from rdkit import Chem
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.stats import gaussian_kde
from scipy import stats
import seaborn as sns
import pandas as pd
from sklearn.metrics import mean_squared_error,r2_score


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

# def create_ijbonddict(mol):
#     """Create a dictionary, which each key is a node ID
#     and each value is the tuples of its neighboring node
#     and bond (e.g., single and double) IDs."""
#     # bond_dict = defaultdict(lambda: len(bond_dict))
#     i_jbond_dict = defaultdict(lambda: [])
#     for b in mol.GetBonds():
#         i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
#         print(str(b.GetBondType()))
#         bond = bond_dict[str(b.GetBondType())]
#         print(bond)
#         # bond = bond_dict.get(str(b.GetBondType()))
#         # try :
#         #     bond = bond_dict[str(b.GetBondType())]
#         # except :
#         #     bond_dict[str(b.GetBondType())] = 0
#         #     bond = bond_dict[str(b.GetBondType())]

#         i_jbond_dict[i].append((j, bond))
#         i_jbond_dict[j].append((i, bond))
#     return i_jbond_dict

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

def main() :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        Kcat_data = json.load(infile)

    # with open('../species/Saccharomyces_cerevisiaeForKcatPrediction2.txt', 'r') as infile :
    #     lines = infile.readlines()[1:]

    # print(len(lines))  # 6291
    # # print(lines[1])

    # # proteinSeq = get_refSeq()

    fingerprint_dict = model.load_pickle('../../Data/input/fingerprint_dict.pickle')
    atom_dict = model.load_pickle('../../Data/input/atom_dict.pickle')
    bond_dict = model.load_pickle('../../Data/input/bond_dict.pickle')
    word_dict = model.load_pickle('../../Data/input/sequence_dict.pickle')
    n_fingerprint = len(fingerprint_dict)
    n_word = len(word_dict)
    print(n_fingerprint)  # 3958
    print(n_word)  # 8542

    radius=2
    ngram=3
    # n_fingerprint = 3958
    # n_word = 8542

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

    # prediction = predictor.predict(inputs)

    i = 0
    # x = list()
    # y = list()
    experimental_values = list()
    predicted_values = list()
    for data in Kcat_data :
        # print(data)
        # print(data['Substrate'])
        if data['Type'] == 'mutant' :
            # print(data)
            i += 1
            print('This is', i, '---------------------------------------')
            smiles = data['Smiles']
            sequence = data['Sequence']
            print(smiles)
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

                value = float(data['Value'])
                print(value)
                print(type(value))
                # y1.append(value)
                experimental_values.append(math.log10(value))

                prediction = predictor.predict(inputs)
                Kcat_log_value = prediction.item()
                Kcat_value = math.pow(2,Kcat_log_value)
                print(Kcat_value)
                print(type(Kcat_value))
                # x1.append(Kcat_value)
                predicted_values.append(math.log10(Kcat_value))

    # correlation, p_value = stats.pearsonr(x, y)
    correlation1, p_value1 = stats.pearsonr(experimental_values, predicted_values)

    # https://blog.csdn.net/u012735708/article/details/84337262?utm_medium=distribute.pc_relevant.none-
    # task-blog-BlogCommendFromMachineLearnPai2-1.pc_relevant_is_cache&depth_1-utm_source=
    # distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-1.pc_relevant_is_cache
    r2 = r2_score(experimental_values,predicted_values)
    rmse = np.sqrt(mean_squared_error(experimental_values,predicted_values))
    print("---------------------")
    print('\n\n')
    # print(correlation)
    print(correlation1)
    print(p_value1)
    print('R2 is', r2)
    print('RMSE is', rmse)

    # Results:
    # 0.8970561077126646
    # 0.0
    # R2 is 0.8031064639769758
    # RMSE is 0.6683890205006177


    allData = pd.DataFrame(list(zip(experimental_values,predicted_values)))
    allData.columns = ['Experimental value', 'Predicted value']

    plt.figure(figsize=(1.5,1.5))

    # To solve the 'Helvetica' font cannot be used in PDF file
    # https://stackoverflow.com/questions/59845568/the-pdf-backend-does-not-currently-support-the-selected-font
    # rc('text', usetex=True) 
    rc('font',**{'family':'serif','serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42
    # plt.rc('text', usetex=True)

    plt.axes([0.12,0.12,0.83,0.83])

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    # http://showteeth.tech/posts/24328.html
    # https://stackoverflow.com/questions/49662964/density-scatter-plot-for-huge-dataset-in-matplotlib
    kcat_values_vstack = np.vstack([experimental_values,predicted_values])
    experimental_predicted = gaussian_kde(kcat_values_vstack)(kcat_values_vstack)

    # plt.scatter(data = allData, x = 'Predicted value', y = 'Experimental value')
    # sns.regplot(data = allData, x = 'Experimental value', y = 'Predicted value', color='#2166ac', scatter_kws={"s": 1})
    ax = plt.scatter(x = experimental_values, y = predicted_values, c=experimental_predicted, s=3, edgecolor='')

    plt.text(-6.7, 6.0, 'r = 0.90', fontweight ="normal", fontsize=6)
    plt.text(-6.7, 5.0, 'p value = 0', fontweight ="normal", fontsize=6)
    plt.text(-6.7, 3.9, 'N = 7,481', fontweight ="normal", fontsize=6)

    plt.text(2, -6, 'Mutant', fontweight ="normal", fontsize=6)

    plt.rcParams['font.family'] = 'Helvetica'

    plt.xlabel("Experimental $k$$_\mathregular{cat}$ value", fontdict={'weight': 'normal', 'fontname': 'Helvetica', 'size': 7}, fontsize=7)
    plt.ylabel('Predicted $k$$_\mathregular{cat}$ value',fontdict={'weight': 'normal', 'fontname': 'Helvetica', 'size': 7},fontsize=7)

    plt.xticks([-8, -6, -4, -2, 0, 2, 4, 6, 8])
    plt.yticks([-8, -6, -4, -2, 0, 2, 4, 6, 8])

    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    # plt.rcParams['text.usetex'] = True

    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    plt.savefig("../../Results/figures/Fig3b.pdf", dpi=400, bbox_inches='tight')


if __name__ == '__main__' :
    main()
