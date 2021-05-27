#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import os
import csv
import math
import subsequence_model
import torch
import json
import pickle
import numpy as np
from rdkit import Chem
from Bio import SeqIO
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.legend_handler import HandlerPathCollection
from scipy import stats
import seaborn as sns
import pandas as pd
from sklearn.metrics import mean_squared_error,r2_score


fingerprint_dict = subsequence_model.load_pickle('../../Data/input/fingerprint_dict.pickle')
atom_dict = subsequence_model.load_pickle('../../Data/input/atom_dict.pickle')
bond_dict = subsequence_model.load_pickle('../../Data/input/bond_dict.pickle')
edge_dict = subsequence_model.load_pickle('../../Data/input/edge_dict.pickle')
word_dict = subsequence_model.load_pickle('../../Data/input/sequence_dict.pickle')

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

# To set a constant marker size in the legend
# https://stackoverflow.com/questions/47115869/how-do-i-change-the-size-of-the-scatter-markers-in-the-legend
marker_size = 25
def update_prop(handle, orig):
    handle.update_from(orig)
    handle.set_sizes([marker_size])

def plot_attention_weights(attention_profiles, wildtype_like_positions, wildtype_decreased_positions, wildtype_like, wildtype_decreased) :
    positions = list()
    weights = list()
    i = 0
    for attention in attention_profiles :
        i += 1
        positions.append(i)
        weights.append(float(attention))

    plt.figure(figsize=(2.0,1.5))

    # To solve the 'Helvetica' font cannot be used in PDF file
    # https://stackoverflow.com/questions/59845568/the-pdf-backend-does-not-currently-support-the-selected-font
    rc('font',**{'family':'serif','serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42

    plt.axes([0.12,0.12,0.83,0.83])

    # plt.rcParams['xtick.direction'] = 'in'
    # plt.rcParams['ytick.direction'] = 'in'

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    plt.plot(positions, weights, color='#A65628', linestyle='--', linewidth=0.75)  # color='k'  color='#A65628'

    # s = [10*4**n for n in range(len(x))]  # change the size according to the detailed number

    # plt.scatter(wildtype_like_positions, wildtype_like, s=3, color='#2166ac', marker='^')  # markersize=1
    # plt.scatter(wildtype_decreased_positions, wildtype_decreased, s=3, color='#b2182b', marker='^')

    print(Counter(wildtype_like_positions))
    print(Counter(wildtype_decreased_positions))
    # Counter({88: 1, 159: 1, 89: 1, 244: 1, 219: 1})
    # Counter({88: 14, 86: 2, 242: 2, 243: 1, 33: 1, 84: 1, 201: 1, 200: 1, 257: 1})
    # for position in wildtype_like_positions :
    #     print(Counter(wildtype_like_positions)[position])

    # Increase scatter marker size
    # https://www.delftstack.com/howto/matplotlib/how-to-set-marker-size-of-scatter-plot-in-matplotlib/
    plt.scatter(wildtype_like_positions, wildtype_like, s=[Counter(wildtype_like_positions)[position]*5 for position in wildtype_like_positions], color='#2166ac', marker='o', label='Wildtype_like')  # marker='o' or marker='s'
    sc = plt.scatter(wildtype_decreased_positions, wildtype_decreased, s=[Counter(wildtype_decreased_positions)[position]*5 for position in wildtype_decreased_positions], color='#b2182b', marker='o', label='Wildtype_decreased')
 
    plt.rcParams['font.family'] = 'Helvetica'
    plt.xlabel('Sequence position', fontsize=7)
    plt.ylabel('Attention weight', fontsize=7)
    # plt.ylabel('Importance contribution', fontsize=7)

    plt.xticks([0,50,100,150,200,250,300])
    plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])

    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    # To set a constant marker size in the legend
    # https://stackoverflow.com/questions/47115869/how-do-i-change-the-size-of-the-scatter-markers-in-the-legend
    plt.legend(handler_map={type(sc): HandlerPathCollection(update_func=update_prop)}, frameon=False, markerscale=2.0, numpoints=1, prop={"size":6})

    # labels = ax.get_legend_handles_labels()[1]
    # plt.legend(handles[0:2], labels[0:2], loc=1, markerscale=2.0, frameon=False, prop={'size': 6})
    # plt.legend(frameon=False, markerscale=2.0, numpoints=1, prop={"size":6})

    plt.savefig("../../Results/figures/Fig3e.pdf", dpi=400, bbox_inches = 'tight')

class Predictor(object):
    def __init__(self, model):
        self.model = model

    def predict(self, data):
        predicted_value,attention_profiles = self.model.forward(data)

        return predicted_value, attention_profiles

def extract_wildtype_mutant() :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        Kcat_data = json.load(infile)

    entry_keys = list()
    for data in Kcat_data :
        # print(data['ECNumber'])
        # print(data['Substrate'])
        # print(data['Organism'])

        substrate = data['Substrate']
        organism = data['Organism']
        EC = data['ECNumber']
        entry_key = substrate + '&' + organism + '&' + EC
        # print(entry_key.lower())
        entry_keys.append(entry_key)

    entry_dict = dict(Counter(entry_keys))
    # print(entry_dict)

    duplicated_keys = [key for key, value in entry_dict.items() if value > 1]
    # print(duplicated_keys)

    duplicated_dict = {key:value for key, value in entry_dict.items() if value > 1}
    # print(duplicated_dict)
    # https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value
    # print(sorted(duplicated_dict.items(), key=lambda x: x[1], reverse=True)[:30])
    duplicated_list = sorted(duplicated_dict.items(), key=lambda x: x[1], reverse=True)[:30]

    for duplicated in duplicated_list[:1] :
        # print('The subtrate name:', duplicated[0])
        for data in Kcat_data :
            # duplicated_one_entry = duplicated_list[0].split('&')
            substrate = data['Substrate']
            organism = data['Organism']
            EC = data['ECNumber']
            one_entry = substrate + '&' + organism + '&' + EC
            if one_entry == duplicated[0] :
                enzyme_type = data['Type']
                Kcat_value = data['Value']
                # print('Substrate:', substrate)
                # print('%s enzyme: %s' %(enzyme_type, Kcat_value))
        # print('----'*15+'\n')

    return duplicated_list

def compare_list(mutant, wildtype) :
    different_attentions = list()
    for i in range(0, len(wildtype)) :
        if mutant[i]  !=  wildtype[i] :
            different_attentions.append(mutant[i])
        else :
            continue

    return different_attentions

def compare_mutant_wildtype_sequence(mutant, wildtype) :
    different_positions = list()
    for i in range(0, len(wildtype)) :
        if mutant[i]  !=  wildtype[i] :
            # different_positions.append(mutant[i])
            different_positions.append(i)
        else :
            continue

    return different_positions

def extract_wildtype_kcat(entry) :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        Kcat_data = json.load(infile)

    for data in Kcat_data :
        substrate = data['Substrate']
        organism = data['Organism']
        EC = data['ECNumber']
        one_entry = substrate + '&' + organism + '&' + EC
        if one_entry == entry :
            enzyme_type = data['Type']
            if enzyme_type == 'wildtype' :
                wildtype_kcat = float(data['Value'])

    if wildtype_kcat :
        return wildtype_kcat
    else :
        return None

def extract_wildtype_sequence(entry) :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        Kcat_data = json.load(infile)

    for data in Kcat_data :
        substrate = data['Substrate']
        organism = data['Organism']
        EC = data['ECNumber']
        one_entry = substrate + '&' + organism + '&' + EC
        if one_entry == entry :
            enzyme_type = data['Type']
            if enzyme_type == 'wildtype' :
                wildtype_sequence = data['Sequence']

    if wildtype_sequence :
        return wildtype_sequence
    else :
        return None

def extract_wildtype_attention(wildtype_entry) :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        Kcat_data = json.load(infile)

    fingerprint_dict = subsequence_model.load_pickle('../../Data/input/fingerprint_dict.pickle')
    atom_dict = subsequence_model.load_pickle('../../Data/input/atom_dict.pickle')
    bond_dict = subsequence_model.load_pickle('../../Data/input/bond_dict.pickle')
    word_dict = subsequence_model.load_pickle('../../Data/input/sequence_dict.pickle')
    n_fingerprint = len(fingerprint_dict)
    n_word = len(word_dict)
    # print(n_fingerprint)  # 3958
    # print(n_word)  # 8542

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
    Kcat_model = subsequence_model.KcatPrediction(device, n_fingerprint, n_word, 2*dim, layer_gnn, window, layer_cnn, layer_output).to(device)
    Kcat_model.load_state_dict(torch.load('../../Results/output/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration30', map_location=device))
    # print(state_dict.keys())
    # subsequence_model.eval()
    predictor = Predictor(Kcat_model)

    for data in Kcat_data :
        substrate = data['Substrate']
        organism = data['Organism']
        EC = data['ECNumber']
        enzyme_type = data['Type']
        entry = substrate + '&' + organism + '&' + EC

        if entry == wildtype_entry and enzyme_type == 'wildtype' :
            smiles = data['Smiles']
            sequence = data['Sequence']
            Kcat = data['Value']
            if "." not in smiles and float(Kcat) > 0:

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
                # print('Current kcat value:', value)

                prediction, wildtype_attention_profiles = predictor.predict(inputs)
                # Kcat_log_value = prediction.item()
                # Kcat_value = math.pow(2,Kcat_log_value)
                # plot_attention_weights(attention_profiles)

            return wildtype_attention_profiles, sequence

def output_wildtype_enzyme(wildtype_attention_profiles, sequence) :
    sequence_length = len(sequence)
    attention_weights_length = len(wildtype_attention_profiles)

    print('The length of wildtype enzyme:', sequence_length)
    print('The length of attention weights:', attention_weights_length)
    print(sequence)
    print(wildtype_attention_profiles)

    with open('../../Results/output/supple_wildtype_PNP_attention_weights.tsv', 'w') as outfile :
        i = 0
        items = ['Sequence position', 'Amino acid', 'Attention weight']
        outfile.write('\t'.join(items) + '\n')
        for attention in wildtype_attention_profiles :
            i += 1
            line = [str(i), sequence[i-1], attention]
            outfile.write('\t'.join(line) + '\n')

def wildtype_like_decreased_info() :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        Kcat_data = json.load(infile)

    wildtype_mutant_entries = extract_wildtype_mutant()

    fingerprint_dict = subsequence_model.load_pickle('../../Data/input/fingerprint_dict.pickle')
    atom_dict = subsequence_model.load_pickle('../../Data/input/atom_dict.pickle')
    bond_dict = subsequence_model.load_pickle('../../Data/input/bond_dict.pickle')
    word_dict = subsequence_model.load_pickle('../../Data/input/sequence_dict.pickle')
    n_fingerprint = len(fingerprint_dict)
    n_word = len(word_dict)
    # print(n_fingerprint)  # 3958
    # print(n_word)  # 8542

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
    Kcat_model = subsequence_model.KcatPrediction(device, n_fingerprint, n_word, 2*dim, layer_gnn, window, layer_cnn, layer_output).to(device)
    Kcat_model.load_state_dict(torch.load('../../Results/output/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration30', map_location=device))
    # print(state_dict.keys())
    # subsequence_model.eval()
    predictor = Predictor(Kcat_model)

    print('It\'s time to start the prediction!')
    print('-----------------------------------')

    i = 0
    alldata = dict()
    alldata['type'] = list()
    alldata['entry'] = list()
    alldata['weights'] = list()


    for wildtype_mutant_entry in wildtype_mutant_entries :
        entry_names = wildtype_mutant_entry[0].split('&')
        # print('This entry is:', entry_names)
        # print('The total amount of wildtype and variant enzymes in the entry is:', wildtype_mutant_entry[1])

        experimental_values = list()
        predicted_values = list()
        wildtype_like = list()
        wildtype_decreased = list()
        wildtype_like_positions = list()
        wildtype_decreased_positions = list()

        if entry_names[0] == 'Inosine' :
            print('This entry is:', entry_names)
            for data in Kcat_data :
                # print(data)
                # print(data['Substrate'])
                substrate = data['Substrate']
                organism = data['Organism']
                EC = data['ECNumber']
                entry = substrate + '&' + organism + '&' + EC

                if entry == wildtype_mutant_entry[0] :
                    wildtype_kcat = extract_wildtype_kcat(entry)
                    wildtype_sequence = extract_wildtype_sequence(entry)
                    wildtype_attention_profiles =extract_wildtype_attention(entry)[0]
                    # print(len(wildtype_attention_profiles))
                    # print(wildtype_attention_profiles)
                    # print('wildtype kcat:', wildtype_kcat)
                    # print(data)
                    # if wildtype_kcat :
                    i += 1
                    # print('This is', i, '---------------------------------------')
                    smiles = data['Smiles']
                    sequence = data['Sequence']
                    enzyme_type = data['Type']
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
                        # print('Current kcat value:', value)
                        normalized_value = value/wildtype_kcat
                        # print('%.2f' % normalized_value)
                        # print(type(value))
                        # print(type(normalized_value))
                        experimental_values.append(math.log10(value))

                        prediction, attention_profiles = predictor.predict(inputs)

                        # different_attentions = compare_list(attention_profiles, wildtype_attention_profiles)
                        # different_weights = compare_list(attention_profiles, wildtype_attention_profiles)
                        different_positions = compare_mutant_wildtype_sequence(sequence, wildtype_sequence)

                        entry_name = entry_names[0]
                        if normalized_value >= 0.5 and normalized_value < 2.0 :
                            # for weight in different_weights :
                            for position in different_positions :
                                wildtype_like_positions.append(position+1)
                                wildtype_like.append(float(wildtype_attention_profiles[position]))
                                alldata['type'].append('Wildtype_like')
                                alldata['entry'].append(entry_name)
                                alldata['weights'].append(float(wildtype_attention_profiles[position]))

                        if normalized_value < 0.5 :
                            # for weight in different_weights :
                            for position in different_positions :
                                wildtype_decreased_positions.append(position+1)
                                wildtype_decreased.append(float(wildtype_attention_profiles[position]))
                                alldata['type'].append('Wildtype_decreased')
                                alldata['entry'].append(entry_name)
                                alldata['weights'].append(float(wildtype_attention_profiles[position]))

            print('wildtype_like_positions:', wildtype_like_positions)
            print('wildtype_decreased_positions:', wildtype_decreased_positions)
            # print(set(wildtype_decreased_positions))
            # print(len(wildtype_like_positions))
            # print(len(wildtype_decreased_positions))
            print('Attention weights in wildtype_like:', wildtype_like)
            print('Attention weights in wildtype_decreased:', wildtype_decreased)
            # print(set(wildtype_decreased))

            return wildtype_like_positions, wildtype_decreased_positions, wildtype_like, wildtype_decreased

def main() :
    substrate, organism, EC = ('Inosine', 'Homo sapiens', '2.4.2.1') 
    entry = substrate + '&' + organism + '&' + EC
    wildtype_attentions, sequence = extract_wildtype_attention(entry)
    # output_wildtype_enzyme(wildtype_attentions, sequence)
    wildtype_like_positions, wildtype_decreased_positions, wildtype_like, wildtype_decreased = wildtype_like_decreased_info()
    plot_attention_weights(wildtype_attentions, wildtype_like_positions, wildtype_decreased_positions, wildtype_like, wildtype_decreased)


if __name__ == '__main__' :
    main()
