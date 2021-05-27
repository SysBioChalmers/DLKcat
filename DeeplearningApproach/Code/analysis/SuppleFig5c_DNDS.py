#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import os
import csv
import math
import xlrd
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from Bio import SeqIO
from collections import defaultdict
from scipy import stats
from scipy.stats import ranksums
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import rc


def getIndex() :
    # Downloaded the data orthomcl_SeqIDs_index.txt from the Figshare data repository (10.6084/m9.figshare.5854692; 
    # https://figshare.com/articles/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692)
    # get the ortholog accoding to protein sequence id, that means Alloascoidea_hylecoeti@Seq_1 as the key, 0_0 as the value
    with open("../../Data/directory/to/orthomcl_SeqIDs_index.txt", "r") as indexFile :
        indexs = indexFile.readlines()

    indexSeqId = dict()
    for index in indexs :
        index_Seq = index.strip().split(": ")
        indexSeqId[index_Seq[0]] = index_Seq[1]

    return indexSeqId

def getIndex2() :
    # Downloaded the data orthomcl_SeqIDs_index.txt from the Figshare data repository (10.6084/m9.figshare.5854692; 
    # https://figshare.com/articles/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692)
    # get the ortholog accoding to protein sequence id, that means Alloascoidea_hylecoeti@Seq_1 as the key, 0_0 as the value
    with open("../../Data/directory/to/orthomcl_SeqIDs_index.txt", "r") as indexFile :
        indexs = indexFile.readlines()

    indexSeqId = dict()
    for index in indexs :
        index_Seq = index.strip().split(": ")
        indexSeqId[index_Seq[1]] = index_Seq[0]

    return indexSeqId 

def getOrthologIndex() :
    # Downloaded the data orthomcl_clusters.txt from the Figshare data repository (10.6084/m9.figshare.5854692; 
    # https://figshare.com/articles/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692)
    with open("../../Data/directory/to/orthomcl_clusters.txt", "r") as orthologFile :
        orthologs = orthologFile.readlines()

    orthologIndex = dict()
    for ortholog in orthologs :
        ortholog_Index = ortholog.strip().split(" ")
        # orthologIndex = {'OG1001': {'328_2397', '189_1696', '279_256',.....}}
        ortholog = ortholog_Index[0][:-1]
        
        orthologIndex[ortholog] = ortholog_Index[1:]

    return orthologIndex

def getOrthologIndex2() :
    # Downloaded the data orthomcl_clusters.txt from the Figshare data repository (10.6084/m9.figshare.5854692; 
    # https://figshare.com/articles/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692)
    with open("../../Data/directory/to/orthomcl_clusters.txt", "r") as orthologFile :
        orthologs = orthologFile.readlines()

    orthologIndex = dict()
    for ortholog in orthologs :
        ortholog_Index = ortholog.strip().split(" ")
        # orthologIndex = {'OG1001': {'328_2397', '189_1696', '279_256',.....}}
        ortholog = ortholog_Index[0][:-1]
        
        for index in ortholog_Index[1:] :
            orthologIndex[index] = ortholog
    # print(orthologIndex)  # {'302_3224': 'OG1000', '317_1502': 'OG1000', '318_1938': 'OG1001', '320_301': 'OG1001', '325_5347': 'OG1001'}

    return orthologIndex

def get_organisms() :
    filenames = os.listdir('../../Data/MLKCATRESULT/')
    filenames = [filename.split('ForKcat')[0] for filename in filenames if filename.endswith('.txt')]
    print(len(filenames)) # 343
    # print(filenames[:3])  # ['yHMPu5000035645_Yarrowia_divulgata', 'Saccharomyces_uvarum', 'Cyberlindnera_jadinii']
    return filenames

def getDNDS_all() :
    with open('../../Data/gene_dn_ds_03_02.csv', 'r') as infile :
        lines = infile.readlines()[1:]
    # print(len(lines1))

    dnds_dict = dict()
    for line in lines :
        data = line.strip().split(',')
        # print(data)
        if data[2] :
            OG_line = line.strip().split(',')[1].split('.')[0]
            dnds_score = line.strip().split(',')[2]
            # print(dnds_score)

            dnds_dict[OG_line] = float(dnds_score)

    return dnds_dict

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
   
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def species_clade() :
    with open("../../../BayesianApporach/Data/343_phenotype_clade.tsv", 'r') as infile :
        lines = infile.readlines()[1:]

    species = list()
    clade = list()

    for line in lines :
        data = line.strip().split('\t')
        species.append(data[0])
        clade.append(data[1])

    print(species[-3:])
    print(clade[-3:])

    species_clade = dict(zip(species,clade))
    # print(len(species_clade))
    return species_clade

def main() :

    SeqIdIndex = getIndex2()
    IndexOrtholog = getOrthologIndex2()
    ortholog_DNDS = getDNDS_all()
    organisms = get_organisms()
    species_clades = species_clade()

    # organisms = ['Saccharomyces_cerevisiae','Yarrowia_lipolytica','Kluyveromyces_marxianus','Kluyveromyces_lactis','Komagataella_pastoris','Lachancea_kluyveri','Candida_albicans']
    # organisms = ['Saccharomyces_cerevisiae','Yarrowia_lipolytica','Kluyveromyces_marxianus','Kluyveromyces_lactis','Lachancea_kluyveri', 
    #              'Saccharomyces_uvarum']

    all_clades = ['Outgroup', 'Lipomycetaceae', 'Trigonopsidaceae', 'Dipodascaceae/Trichomonascaceae', 'Alloascoideaceae', 'Sporopachydermia clade',
                    'Pichiaceae', 'CUG-Ala', 'CUG-Ser1', 'CUG-Ser2', 'Phaffomycetaceae', 'Saccharomycodaceae', 'Saccharomycetaceae']

    all_clades_order = {'Outgroup':1, 'Lipomycetaceae':2, 'Trigonopsidaceae':3, 'Dipodascaceae/Trichomonascaceae':4, 'Alloascoideaceae':5, 'Sporopachydermia clade':6,
                    'Pichiaceae':7, 'CUG-Ala':8, 'CUG-Ser1':9, 'CUG-Ser2':10, 'Phaffomycetaceae':11, 'Saccharomycodaceae':12, 'Saccharomycetaceae':13}

    alldata = dict()
    alldata['type'] = list()
    alldata['clade'] = list()
    alldata['Kcat_value'] = list()
    counts_cluster_1 = list()
    counts_cluster_2 = list()

    for clade in all_clades :
        for organism in organisms :
            if species_clades[organism.lower()] == clade :
                # print('This is', organism)

                with open('../prediction/343species_0115/%s_PredictionResults.txt' % organism, 'r') as infile :
                    lines = infile.readlines()

                # seqIds_values = dict()
                seq_kcat = list()

                for line in lines[1:] : 
                    seqIds_values = dict()
                    # seq_kcat = list()
                    data = line.strip('\n').split('\t')
                    smiles = data[4].split(';')
                    seqIds = data[5].split(';')
                    values = data[-1].split(';')

                    if values :
                        for i, seqId in enumerate(seqIds) :
                            for value in values :
                                if value :
                                    try :
                                        Kcats = value.split(',')
                                        if Kcats[i] != '#' :
                                            seqIds_values[seqId].append(float(Kcats[i]))
                                        else :
                                            pass
                                    except :
                                        Kcat = list()
                                        if Kcats[i] != '#' :
                                            Kcat.append(float(Kcats[i]))
                                            seqIds_values[seqId] = Kcat
                                        else :
                                            pass
                    # print(seqIds_values)

                    for seqId, value in seqIds_values.items() :
                        max_value = max(value)
                        seq_kcat.append((seqId, max_value))

                # print(len(seq_kcat))  # 5876
                # print(seq_kcat[:3])
                seq_kcat_no_copy = list(set(seq_kcat))
                # print(len(seq_kcat_no_copy))  # 2992
                # print(seq_kcat_no_copy[:3])
                # print(len(seqIds_values))

                for item in seq_kcat_no_copy :
                    seqId = item[0]
                    max_value = item[1]
                    index = SeqIdIndex[seqId]
                    ortholog = IndexOrtholog[index]
                    try :
                        dnds = float(ortholog_DNDS[ortholog])
                        # print(dnds)
                        kcatValue = math.log10(max_value)

                        if dnds>0 and dnds<=0.15 :
                            # alldata['type'].append('Conserved')
                            alldata['type'].append('dN/dS <= 0.15')
                            # alldata['clade'].append(clade)
                            alldata['clade'].append(all_clades_order[clade])
                            alldata['Kcat_value'].append(kcatValue)
                        else :
                            # alldata['type'].append('Non-conserved')
                            alldata['type'].append('dN/dS > 0.15')
                            # alldata['clade'].append(clade)
                            alldata['clade'].append(all_clades_order[clade])
                            alldata['Kcat_value'].append(kcatValue)

                    except :
                        continue

    # All clades:
    # ['Saccharomycodaceae', 'CUG-Ser1', 'CUG-Ser2', 'Dipodascaceae/Trichomonascaceae', 'Pichiaceae', 'Lipomycetaceae', 'Alloascoideaceae', 
    # 'Sporopachydermia clade', 'Saccharomycetaceae', 'Trigonopsidaceae', 'Phaffomycetaceae', 'CUG-Ala', 'Outgroup']

    # print(alldata['type'][:3])
    # print(alldata['clade'][:3])
    # print(alldata['Kcat_value'][:3])

    # print(len(alldata['type']))
    # print(len(alldata['clade']))
    # print(len(alldata['Kcat_value']))

    allData = pd.DataFrame(alldata)
    # print(type(allData))

    # for clade in all_clades :
    #     print('This is the clade:', clade)
    #     cluster_1 = list()
    #     cluster_2 = list()
    #     # types = allData.iloc[:,1]
    #     # print(len(types))
    #     # print(types[:3])
    #     # for clade_type in types :
    #     #     if clade_type == clade :
    #     for row_index, row in allData.iterrows() :
    #         if row['clade'] == clade and row['type'] == 'dN/dS <= 0.15' :
    #             # print(row['Kcat_value'])
    #             cluster_1.append(row['Kcat_value'])
    #         if row['clade'] == clade and row['type'] == 'dN/dS > 0.15' :
    #             # print(row['Kcat_value'])
    #             cluster_2.append(row['Kcat_value'])

    #     stat, p_value = ranksums(cluster_1,cluster_2)
    #     print('The P_value between the two dN/dS clusters is:', p_value)

        # Results :
        # This is the clade: Outgroup
        # The P_value between the two dN/dS clusters is: 1.6243302130328922e-61
        # This is the clade: Lipomycetaceae
        # The P_value between the two dN/dS clusters is: 7.879651158117646e-67
        # This is the clade: CUG-Ser1
        # The P_value between the two dN/dS clusters is: 0.0
        # This is the clade: Phaffomycetaceae
        # The P_value between the two dN/dS clusters is: 3.6142539325434596e-75
        # This is the clade: Dipodascaceae/Trichomonascaceae
        # The P_value between the two dN/dS clusters is: 8.512690063502762e-117
        # This is the clade: Trigonopsidaceae
        # The P_value between the two dN/dS clusters is: 4.157606980523744e-25
        # This is the clade: Saccharomycodaceae
        # The P_value between the two dN/dS clusters is: 7.633228849794443e-58
        # This is the clade: Sporopachydermia clade
        # The P_value between the two dN/dS clusters is: 1.2098972408782565e-07
        # This is the clade: Pichiaceae
        # The P_value between the two dN/dS clusters is: 1.8480664486765028e-291
        # This is the clade: CUG-Ser2
        # The P_value between the two dN/dS clusters is: 2.801972561349211e-19
        # This is the clade: CUG-Ala
        # The P_value between the two dN/dS clusters is: 3.1431166089138013e-38
        # This is the clade: Saccharomycetaceae
        # The P_value between the two dN/dS clusters is: 3.553336840154913e-298
        # This is the clade: Alloascoideaceae
        # The P_value between the two dN/dS clusters is: 0.0002450681317830253

    plt.figure(figsize=(2.5, 2.0))
    # To solve the 'Helvetica' font cannot be used in PDF file
    # https://stackoverflow.com/questions/59845568/the-pdf-backend-does-not-currently-support-the-selected-font
    rc('font',**{'family':'serif','serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42

    plt.axes([0.12,0.12,0.83,0.83])
    
    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)
    plt.tick_params(which='major',width=0.4)

    palette = {"dN/dS <= 0.15": '#b2182b', "dN/dS > 0.15": '#2166ac'}

    ax = sns.boxplot(data=alldata, x="clade", y="Kcat_value", hue="type",
            palette=palette, showfliers=False, linewidth=0.5)

    # https://stackoverflow.com/questions/58476654/how-to-remove-or-hide-x-axis-label-from-seaborn-boxplot
    # plt.xlabel(None) will remove the Label, but not the ticks. 
    ax.set(xlabel=None)
    # ax.set(xticks=None)

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

        # if i % 2 == 0:
        #     col = '#b2182b'
        # else:
        #     col = '#2166ac'

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

    for i in range(len(all_clades)) :
        plt.text(i-0.3, 2.6, '***', fontweight ="normal", fontsize=6)

    plt.ylabel("$k$$_\mathregular{cat}$ value", fontname='Helvetica', fontsize=7)

    plt.xticks(rotation=30,ha='right')
    plt.ylim(-2,5)
    plt.yticks([-2,-1,0,1,2,3,4,5])
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

    # https://blog.csdn.net/weixin_38314865/article/details/88633880
    plt.savefig("../../Results/figures/SuppleFig5c.pdf", dpi=400, bbox_inches = 'tight')


if __name__ == '__main__' :
    species_clade()
    # main()
