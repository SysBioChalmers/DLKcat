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
from scipy import stats

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

def get_refSeq() :
    # get the protein sequence accoding to protein sequence id
    # Note that the file 343taxa_proteins.fasta could be downloaded from url: https://figshare.com/articles/dataset/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692?file=13083299
    # Then change the following directory to the directory of 343 protein fasta file
    with open("/directory/to/343taxa_proteins.fasta", "r") as handleGene :
        proteinSeq = dict()
        for record in SeqIO.parse(handleGene, "fasta") :
    # ['__add__', '__bool__', '__class__', '__contains__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__',
    # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__le___', '__len__', '__lt__', 
    # '__module__', '__ne__', '__new__', '__nonzero__', '__radd__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', 
    # '__subclasshook__', '__weakref__', '_per_letter_annotations', '_seq', '_set_per_letter_annotations', '_set_seq', 'annotations', 'dbxrefs', 'description', 
    # 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'translate', 'upper']
            # if record.id.startswith("Candida_albicans") :
            # if record.id == gene :
            proteinSeq[record.id] = str(record.seq)
        # print("The protein number of %s is: %d" % (gene,len(proteinSeq)))

    return proteinSeq

def get_organisms() :
    filenames = os.listdir('../species/MLKCATRESULT/')
    filenames = [filename.split('ForKcat')[0] for filename in filenames if filename.endswith('.txt')]
    print(len(filenames)) # 343
    # print(filenames[:3])  # ['yHMPu5000035645_Yarrowia_divulgata', 'Saccharomyces_uvarum', 'Cyberlindnera_jadinii']
    return filenames

class Predictor(object):
    def __init__(self, model):
        self.model = model

    def predict(self, data):
        predicted_value = self.model.forward(data)

        return predicted_value

def main() :

    proteinSeq = get_refSeq()

    fingerprint_dict = model.load_pickle('../../Data/input/fingerprint_dict.pickle')
    atom_dict = model.load_pickle('../../Data/input/atom_dict.pickle')
    bond_dict = model.load_pickle('../../Data/input/bond_dict.pickle')
    edge_dict = model.load_pickle('../../Data/input/edge_dict.pickle')
    word_dict = model.load_pickle('../../Data/input/sequence_dict.pickle')
    n_fingerprint = len(fingerprint_dict)
    n_word = len(word_dict)
    n_edge = len(edge_dict)

    radius=2
    ngram=3

    dim=10
    # dim=5
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

    organisms = get_organisms()

    i = 0
    for organism in organisms : 
    # if organism == 'Saccharomyces_paradoxus' :   # Saccharomyces_cerevisiae
        i += 1
        print('This is', i, '---------------------------------------')
        print(organism)

        # with open('../species/MLKCATRESULT/%sForKcatPrediction.txt' % organism, 'r') as infile :
        with open('../../Data/input/kcatpredictionfile/%sForKcatPrediction.txt' % organism, 'r') as infile :
            lines = infile.readlines()

        print(len(lines))  # 6291
        print(lines[:2])
        print('--'*20+'\n')

        # Create the directory '343 species' under the 'Results/output' directory 
        # The generated prediction results for 343 species are stored in our zenodo: https://doi.org/10.5281/zenodo.5797013
        file =open('../../Results/output/343species/%s_PredictionResults.txt' % organism, 'w')

        file.write(lines[0].strip() + '\t%s\n' % 'Kcat value (substrate first)')


        for line in lines[1:] :  # [1:]

            data = line.strip('\n').split('\t')
            # print(data)
            # print(len(data))
            # print(data[4])  # Smiles
            # print(data[5])  # protein ID
            smiles_info = data[4].split(';')
            sequence_info = data[5].split(';')
            # print(smiles_info)
            # print(sequence_info)

            if len(smiles_info) :
                file.write(line.strip() + '\t')

            n = 0
            for smiles in smiles_info :
                if n :
                    file.write(';')
                Kcat_values = list()
                n += 1
                if smiles and "." not in smiles :
                    for sequence_id in sequence_info :
                        if sequence_id :
                            print(sequence_id)
                            # print(proteinSeq[sequence_id])
                            if organism != 'Saccharomyces_cerevisiae' :
                                sequence = proteinSeq[sequence_id]
                            else :
                                sequence = proteinSeq['Saccharomyces_cerevisiae@'+sequence_id]
                            # print(smiles)
                            if "." not in smiles :
                                # i += 1
                                # print('This is',i)
                                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                                try :
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
                                    # try :
                                    prediction = predictor.predict(inputs)
                                    Kcat_log_value = prediction.item()
                                    Kcat_value = math.pow(2,Kcat_log_value)
                                    # Kcat_value = math.pow(10,Kcat_log_value)
                                    
                                    # print(Kcat_log_value)
                                    # print(Kcat_value)
                                    # print(type(Kcat_value))
                                    print('%.4f' %(Kcat_value))
                                except :
                                    Kcat_value = '#'

                                # Kcat_values.append('%.4f' %(Kcat_value))

                                if type(Kcat_value) == float :
                                    Kcat_values.append('%.4f' %(Kcat_value))
                                else :
                                    Kcat_values.append(Kcat_value)

                    # if len(Kcat_values) > 1 :
                    added_content = ','.join(Kcat_values)
                    file.write(added_content)

            file.write('\n')
        file.close()


if __name__ == '__main__' :
    main()
