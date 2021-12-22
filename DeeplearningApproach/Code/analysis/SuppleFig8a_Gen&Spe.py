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

def plot_spe_gen() :
    with open('../../../BayesianApporach/Results/kcat_gen_spe.txt') as infile :
        alllines = infile.readlines()[1:]

    # all_clades = ['Outgroup', 'Lipomycetaceae', 'Trigonopsidaceae', 'Dipodascaceae/Trichomonascaceae', 'Alloascoideaceae', 'Sporopachydermia clade',
    #                 'Pichiaceae', 'CUG-Ala', 'CUG-Ser1', 'CUG-Ser2', 'Phaffomycetaceae', 'Saccharomycodaceae', 'Saccharomycetaceae']

    alldata = dict()
    alldata['type'] = list()
    alldata['clade'] = list()
    alldata['Kcat_value'] = list()

    for line in alllines :
        kcatValue = float(line.strip().split('\t')[0])
        data_type = line.strip().split('\t')[1]
        clade_order = int(line.strip().split('\t')[2])

        if data_type == 'gen' :
            alldata['type'].append('Generalist')
            alldata['clade'].append(clade_order)
            alldata['Kcat_value'].append(kcatValue)

    for line in alllines :
        kcatValue = float(line.strip().split('\t')[0])
        data_type = line.strip().split('\t')[1]
        clade_order = int(line.strip().split('\t')[2])

        if data_type == 'spe' :
            alldata['type'].append('Specialist')
            alldata['clade'].append(clade_order)
            alldata['Kcat_value'].append(kcatValue)

    allData = pd.DataFrame(alldata)
    # print(type(allData))

    for clade in range(1,14) :
        print('This is the clade:', clade)
        cluster_1 = list()
        cluster_2 = list()
        # types = allData.iloc[:,1]
        # print(len(types))
        # print(types[:3])
        # for clade_type in types :
        #     if clade_type == clade :
        for row_index, row in allData.iterrows() :
            if row['clade'] == clade and row['type'] == 'Specialist' :
                # print(row['Kcat_value'])
                cluster_1.append(row['Kcat_value'])
            if row['clade'] == clade and row['type'] == 'Generalist' :
                # print(row['Kcat_value'])
                cluster_2.append(row['Kcat_value'])

        stat, p_value = ranksums(cluster_1,cluster_2)
        print('The P_value between the two clusters is:', p_value)

        # Results :
        # This is the clade: 1
        # The P_value between the two clusters is: 1.7089850523355335e-33
        # This is the clade: 2
        # The P_value between the two clusters is: 5.4226066209879584e-24
        # This is the clade: 3
        # The P_value between the two clusters is: 9.280288602757499e-20
        # This is the clade: 4
        # The P_value between the two clusters is: 7.655748855812308e-177
        # This is the clade: 5
        # The P_value between the two clusters is: 1.2522285545479606e-08
        # This is the clade: 6
        # The P_value between the two clusters is: 2.2092138972003236e-10
        # This is the clade: 7
        # The P_value between the two clusters is: 6.9148612803985e-238
        # This is the clade: 8
        # The P_value between the two clusters is: 4.720528655096346e-15
        # This is the clade: 9
        # The P_value between the two clusters is: 0.0
        # This is the clade: 10
        # The P_value between the two clusters is: 3.7033670928277406e-22
        # This is the clade: 11
        # The P_value between the two clusters is: 1.2089161920329716e-127
        # This is the clade: 12
        # The P_value between the two clusters is: 4.4365830029574404e-39
        # This is the clade: 13
        # The P_value between the two clusters is: 0.0

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

    palette = {"Specialist": '#b2182b', "Generalist": '#2166ac'}

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

    for i in range(13) :
        plt.text(i-0.3, 2.95, '***', fontweight ="normal", fontsize=6)

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
    plt.savefig("../../Results/figures/SuppleFig8a.pdf", dpi=400, bbox_inches = 'tight')


if __name__ == '__main__' :
    plot_spe_gen()
