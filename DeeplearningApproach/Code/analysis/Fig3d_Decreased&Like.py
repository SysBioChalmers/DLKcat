#!/usr/bin/python
# coding: utf-8

import os
import math
import model
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
from scipy import stats
import seaborn as sns
import pandas as pd
from scipy.stats import ranksums
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

def compare_prediction_wildtype_mutant() :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        Kcat_data = json.load(infile)

    wildtype_mutant_entries = extract_wildtype_mutant()

    fingerprint_dict = model.load_pickle('../../Data/input/fingerprint_dict.pickle')
    atom_dict = model.load_pickle('../../Data/input/atom_dict.pickle')
    bond_dict = model.load_pickle('../../Data/input/bond_dict.pickle')
    word_dict = model.load_pickle('../../Data/input/sequence_dict.pickle')
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
    Kcat_model = model.KcatPrediction(device, n_fingerprint, n_word, 2*dim, layer_gnn, window, layer_cnn, layer_output).to(device)
    Kcat_model.load_state_dict(torch.load('../../Results/output/all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50', map_location=device))
    # print(state_dict.keys())
    # model.eval()
    predictor = Predictor(Kcat_model)

    print('It\'s time to start the prediction!')
    print('-----------------------------------')

    # prediction = predictor.predict(inputs)

    i = 0
    alldata = dict()
    alldata['type'] = list()
    alldata['entry'] = list()
    alldata['kcat_value'] = list()


    for wildtype_mutant_entry in wildtype_mutant_entries :
        entry_names = wildtype_mutant_entry[0].split('&')
        # print('This entry is:', entry_names)
        # print('The total amount of wildtype and variant enzymes in the entry is:', wildtype_mutant_entry[1])

        experimental_values = list()
        predicted_values = list()
        wildtype_like = list()
        wildtype_decreased = list()

        if entry_names[0] in ['7,8-Dihydrofolate', 'Glycerate 3-phosphate', 'L-Aspartate', 'Penicillin G', 'Inosine', 'Isopentenyl diphosphate'] :
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

                        prediction = predictor.predict(inputs)
                        Kcat_log_value = prediction.item()
                        Kcat_value = math.pow(2,Kcat_log_value)
                        # print(Kcat_value)
                        # print('%.2f' % normalized_value)
                        # print(type(Kcat_value))
                        predicted_values.append(math.log10(Kcat_value))

                        # entry_names = wildtype_mutant_entry[0].split('&')
                        # entry_name = entry_names[0] + '&' + entry_names[2]
                        entry_name = entry_names[0]
                        if normalized_value >= 0.5 and normalized_value < 2.0 :
                            wildtype_like.append(math.log10(Kcat_value))
                            alldata['type'].append('Wildtype_like')
                            alldata['entry'].append(entry_name)
                            alldata['kcat_value'].append(math.log10(Kcat_value))
                        if normalized_value < 0.5 :
                            wildtype_decreased.append(math.log10(Kcat_value))
                            alldata['type'].append('Wildtype_decreased')
                            alldata['entry'].append(entry_name)
                            alldata['kcat_value'].append(math.log10(Kcat_value))

            if wildtype_like and wildtype_decreased :
                p_value = ranksums(wildtype_like, wildtype_decreased)[1]
                print('The amount of wildtype_like:', len(wildtype_like))
                print('The amount of wildtype_decreased:', len(wildtype_decreased))
                print('P value is:', p_value)
                print('\n')

            correlation1, p_value1 = stats.pearsonr(experimental_values, predicted_values)

            # https://blog.csdn.net/u012735708/article/details/84337262?utm_medium=distribute.pc_relevant.none-
            # task-blog-BlogCommendFromMachineLearnPai2-1.pc_relevant_is_cache&depth_1-utm_source=
            # distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-1.pc_relevant_is_cache
            r2 = r2_score(experimental_values,predicted_values)
            rmse = np.sqrt(mean_squared_error(experimental_values,predicted_values))
            # print("---------------------")
            # print('\n\n')
            # print(correlation)
            print('r is %.4f' % correlation1)
            print('P value is', p_value1)
            # print('R2 is %.4f' % r2)
            # print('RMSE is %.4f' % rmse)
            # print('-----'*10 + '\n')


    # Plot the boxplot figures between the wildtype_like and wildtype_decreased
    allData = pd.DataFrame(alldata)
    # print(type(allData))

    plt.figure(figsize=(1.5, 1.5))
    # To solve the 'Helvetica' font cannot be used in PDF file
    # https://stackoverflow.com/questions/59845568/the-pdf-backend-does-not-currently-support-the-selected-font
    rc('font',**{'family':'serif','serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42

    plt.axes([0.12,0.12,0.83,0.83])
    
    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)
    plt.tick_params(which='major',width=0.4)

    # rectangular box plot
    palette = {"Wildtype_like": '#2166ac', "Wildtype_decreased": '#b2182b'}

    # for ind in allData.index:
    #     allData.loc[ind,'entry'] = '${0}$'.format(allData.loc[ind,'entry'])

    ax = sns.boxplot(data=alldata, x="entry", y="kcat_value", hue="type",
            palette=palette, showfliers=False, linewidth=0.5)  # boxprops=dict(alpha=1.0)

    ax = sns.stripplot(data=alldata, x="entry", y="kcat_value", hue="type", jitter=0.3,  
            palette=palette, size=1.3, dodge=True)  

    # https://stackoverflow.com/questions/58476654/how-to-remove-or-hide-x-axis-label-from-seaborn-boxplot
    # plt.xlabel(None) will remove the Label, but not the ticks. 
    ax.set(xlabel=None)

    for patch in ax.artists:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.3))

    # print(ax.artists)
    # print(ax.lines)
    # print(len(ax.lines))
    # https://cduvallet.github.io/posts/2018/03/boxplots-in-python
    for i, artist in enumerate(ax.artists):
        # print(i)

        if i % 2 == 0:
            col = '#2166ac'
        else:
            col = '#b2182b'

        # This sets the color for the main box
        artist.set_edgecolor(col)

        # Each box has 5 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        for j in range(i*5,i*5+5):
            # print(j)
            line = ax.lines[j]
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)
    handles = [ax.artists[0], ax.artists[1]]

    # for tick in ax.get_xticklabels() :
    #     tick.set_rotation(30)

    plt.rcParams['font.family'] = 'Helvetica'

    plt.text(-0.2, 1.5, '***', fontweight ="normal", fontsize=6)
    plt.text(1, 1.0, '*', fontweight ="normal", fontsize=6)
    plt.text(1.9, 1.0, '**', fontweight ="normal", fontsize=6)
    plt.text(2.9, -0.7, '**', fontweight ="normal", fontsize=6)
    plt.text(3.9, 1.4, '**', fontweight ="normal", fontsize=6)
    plt.text(5, -0.3, '*', fontweight ="normal", fontsize=6)

    plt.ylabel("$k$$_\mathregular{cat}$ value", fontname='Helvetica', fontsize=7)

    plt.xticks(rotation=30,ha='right')
    plt.ylim(-4,4)
    plt.yticks([-4,-2,0,2,4])

    plt.xticks(fontsize=7)
    plt.yticks(fontsize=6)

    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    ax = plt.gca()
    # handles,labels = ax.get_legend_handles_labels()
    labels = ax.get_legend_handles_labels()[1]
    # print(handles)
    # print(labels)
    # specify just one legend
    lgd = plt.legend(handles[0:2], labels[0:2], loc=1, frameon=False, prop={'size': 6})

    # plt.rcParams['font.family'] = 'Helvetica'

    plt.savefig("../../Results/figures/Fig3d.pdf", dpi=400, bbox_inches = 'tight')


if __name__ == '__main__' :
    # extract_wildtype_mutant()
    compare_prediction_wildtype_mutant()
