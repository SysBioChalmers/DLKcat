#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2021-05-14

# This python script is to classify subsystem base on EC number, Kcat and subsystem mapping for 343 species

import os
import json
import math
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import numpy as np
import statsmodels.api as sm


def EC_Kcat() :
    organisms = get_organisms()
    EC_Kcat = dict()

    for organism in organisms:
        # print(organism)
        with open('../../Results/343species_0314/%s_PredictionResults.txt' % organism, 'r') as infile1 :  # 0314  0115
            lines1 = infile1.readlines()

        with open('../../Data/ecres/%s_ec.txt' % organism, 'r') as infile2 :
            lines2 = infile2.readlines()

        rxnseq_ec = dict()
        for line in lines2 :
            data = line.strip('\n').split(',')
            rxnseq = data[0]
            # ec = [EC[2:] for EC in data[2].split('; ') if EC]
            # print(data[2])
            ec = list()
            for EC in data[2].split('; ') :
                if EC.endswith(';') :
                    ec.append(EC[2:-1])
                else :
                    ec.append(EC[2:])
            # print(list(set(ec)))
            # rxnseq_ec[rxnseq] = list(set(ec))
            rxnseq_ec[rxnseq] = ec
        # print(rxnseq_ec)

        # seqIds_values = dict()
        seq_kcat = list()

        for line in lines1[1:] : 
            seqIds_values = dict()
            # seq_kcat = list()
            data = line.strip('\n').split('\t')
            rxnid = data[0]
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
                                    seqIds_values[rxnid].append(float(Kcats[i]))
                                else :
                                    pass
                            except :
                                Kcat = list()
                                if Kcats[i] != '#' :
                                    Kcat.append(float(Kcats[i]))
                                    seqIds_values[rxnid] = Kcat
                                else :
                                    pass
            # print(seqIds_values)

            for seqId, value in seqIds_values.items() :
                max_value = max(value)
                seq_kcat.append((seqId, max_value))

            # for seqId, values in seqIds_values.items() :
            #     for value in values :
            #         seq_kcat.append((seqId, value))

        # print(len(seq_kcat))  # 4957
        # print(seq_kcat[:5])
        seq_kcat_no_copy = list(set(seq_kcat))
        # print(len(seq_kcat_no_copy))  # 4957
        # print(seq_kcat_no_copy[:3])
        # print(len(seqIds_values))

        # for item in seq_kcat_no_copy :
        for item in seq_kcat :
            seqId = item[0]
            max_value = item[1]
            # values = item[1]
            try : 
                EC_Numbers = rxnseq_ec[seqId]
            except :
                continue

            for EC_Number in EC_Numbers :
                try :
                    if EC_Kcat[EC_Number] :
                        value = math.log10(max_value)
                        EC_Kcat[EC_Number].append(value)
                except :
                    Kcat = list()
                    value = math.log10(max_value)
                    Kcat.append(value)
                    EC_Kcat[EC_Number] = Kcat


    # print(len(EC_Kcat))
    # print(list(EC_Kcat.values())[:3])
    return EC_Kcat

def get_organisms() :
    filenames = os.listdir('../../Data/MLKCATRESULT/')
    filenames = [filename.split('ForKcat')[0] for filename in filenames if filename.endswith('.txt')]
    print(len(filenames)) # 343
    # print(filenames[:3])  # ['yHMPu5000035645_Yarrowia_divulgata', 'Saccharomyces_uvarum', 'Cyberlindnera_jadinii']
    return filenames

def EC_subsystem() :
    with open('../../Data/subsystem/module_ec.txt', 'r') as infile :
        datasets = infile.readlines()

    print(len(datasets)) # 2200

    metabolism_types = list()
    types_abbre = {
    'Primary - Carbohydrate & Energy Metabolism': 'Primary-CE',
    'Secondary_other': 'Secondary_other',
    'Intermediate': 'Intermediate', 
    'Secondary': 'Secondary',
    'Primary - amino acids, fatty acids and nucleotides': 'Primary-AFN',
    'x': 'x'
    }

    types_EC = dict()
    for data in datasets :
        # print(data)
        line = data.strip().split('\t')
        # print(line)
        metabolism_types.append(line[2])
        abbre = types_abbre[line[2]]
        # EC_types[line[1][2:]] = line[2]
        try :
            if types_EC[abbre] :
                types_EC[abbre].append(line[1][2:])
        except :
            EC_Number = list()
            EC_Number.append(line[1][2:])
            types_EC[abbre] = EC_Number

    # metabolism_types = list(set(metabolism_types))
    # print(len(metabolism_types)) # 6
    # print(metabolism_types)  
    # # ['Primary - Carbohydrate & Energy Metabolism', 'Secondary_other', 'Intermediate', 'Secondary', 'Primary - amino acids, fatty acids and nucleotides', 'x']

    # print(types_EC)

    print(len(types_EC))

    i = 0
    new_types_EC = dict()
    for types, EC_Number in types_EC.items() :
        # print(len(set(EC_Number)))
        i += len(set(EC_Number))
        new_types_EC[types] = list(set(EC_Number))
        # print('The type of %s has %s unique EC Number.' % (types, len(set(EC_Number))))

    # print('Total EC number is:', i) # 2200, the same with all entries in module.txt

    return new_types_EC

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
   
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def Kcat_subsystem() :
    EC_Kcat_relation = EC_Kcat()
    types_EC = EC_subsystem()
    # print(EC_Kcat_relation)

    types_Kcat = dict()
    for types, EC_Number in types_EC.items() :
        print(types)
        Kcat_values = list()
        for EC in EC_Number :
            try :
                Kcat_values += EC_Kcat_relation[EC]
            except :
                continue
        types_Kcat[types] = Kcat_values

    # # print(len(types_Kcat))
    # for types, Kcat in types_Kcat.items() :
    #     # print('The type of %s has %s Kcat values.' % (types, len(Kcat)))
    #     # print('---'*15)
    #     # print('The median value of %s is %s' %(types, math.log10(median(Kcat))))
    #     print('The median value of %s is %s' %(types, median(Kcat)))

    return types_Kcat

def plot_subsystem_Kcat_counts() :
    types_Kcat = Kcat_subsystem()

    for types, Kcat in types_Kcat.items() :
        print('The type of %s has %s Kcat values.' % (types, len(Kcat)))
        # print('---'*15)
        # print('The median value of %s is %s' %(types, math.log10(median(Kcat))))
        # print('The median value of %s is %s' %(types, median(Kcat)))

    # types = ['Primary-CE', 'Primary-AFN', 'Intermediate', 'Secondary']
    types = ['Primary-CE', 'Primary-AFN', 'Intermediate', 'Secondary', 'Secondary_other']
    counts = [len(types_Kcat[subsystem]) for subsystem in types]

    print(types)
    print(counts)

    plt.figure(figsize=(3.4,2.5))

    # https://juejin.im/post/6858230839986421767
    plt.bar(range(len(types)), counts, tick_label=types, width=0.5, alpha=0.8, color='pink', edgecolor='r')

    # ax = plt.axes()
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    # plt.ylim(0,600)
    # plt.yticks([0,100,200,300,400,500,600])

    # plt.xlabel('Subsystem type',fontsize=12)
    plt.ylabel("Counts", fontsize=12)
    plt.xticks(rotation=30, ha='right')
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    plt.savefig("../../Results/figures/subsystem_Kcat_counts.pdf", dpi=400, bbox_inches='tight')

# https://my.oschina.net/u/4349448/blog/3448306 python code
def plot_subsystem_distribution() :
    types_Kcat = Kcat_subsystem()

    plt.figure(figsize=(2.5,2.0))
    # To solve the 'Helvetica' font cannot be used in PDF file
    # https://stackoverflow.com/questions/59845568/the-pdf-backend-does-not-currently-support-the-selected-font
    rc('font',**{'family':'serif','serif':['Helvetica']})
    plt.rcParams['pdf.fonttype'] = 42

    plt.axes([0.12,0.12,0.83,0.83])

    plt.tick_params(direction='in')
    plt.tick_params(which='major',length=1.5)
    plt.tick_params(which='major',width=0.4)

    plt.rcParams['font.family'] = 'Helvetica'

    types_color = {'Primary-CE': '#F781BF', 'Intermediate': '#4DAF4A', 'Primary-AFN': '#A65628', 'Secondary': '#3182BD'}

    for types, Kcat in types_Kcat.items() :
        if types in ['Primary-CE', 'Intermediate', 'Primary-AFN', 'Secondary'] :
            # print('The type of %s has %s Kcat values.' % (types, len(Kcat)))
            # print('---'*15)  '%.4f' %(Kcat_value)
            # print('The median value on log10 scale of %s is %.4f' %(types, math.log10(median(Kcat))))

            print('The median value of %s is %.2f' %(types, math.pow(10, median(Kcat))))
            # print('The median value in log10 of %s is %.2f' %(types, median(Kcat)))

        # if types == 'Primary-CE' :
        if types in ['Primary-CE', 'Intermediate', 'Primary-AFN', 'Secondary'] :
            ecdf = sm.distributions.ECDF(Kcat)

            x = np.linspace(min(Kcat),max(Kcat),50000)  # 10000
            y = ecdf(x)

            # plt.xlim(1e-4, 1e4)

            # plt.xscale('log')
            # plt.xticks([-2, -1, 0, 1, 2, 3])
            # plt.plot(x,y,linewidth='2',label='Primary-CE')

            plt.plot(x,y,linewidth='1.00',label=types,color=types_color[types])

            # https://blog.csdn.net/weixin_38314865/article/details/115173371?utm_medium=distribute.pc_relevant.none-task-blog-baidujs_baidulandingword-5&spm=1001.2101.3001.4242
            # https://blog.csdn.net/dta0502/article/details/83827345
            plt.axvline(x=median(Kcat),ymin=0,ymax=0.5,linewidth='1.00',linestyle='--',color=types_color[types])

    plt.text(-2.5, 0.89, 'Primary-CE', fontweight ="normal", fontsize=6, color='#F781BF')
    plt.text(-2.5, 0.82, 'Primary-AFN', fontweight ="normal", fontsize=6, color='#A65628')
    plt.text(-2.5, 0.75, 'Secondary', fontweight ="normal", fontsize=6, color='#3182BD')
    plt.text(-2.5, 0.68, 'Intermediate', fontweight ="normal", fontsize=6, color='#4DAF4A')

    plt.rcParams['font.family'] = 'Helvetica'

    plt.xlabel('Predicted $k$$_\mathregular{cat}$ value', fontsize=7)
    plt.ylabel('Cumulative distribution', fontsize=7)

    plt.xlim([-3,3])
    plt.xticks([-3,-2,-1,0,1,2,3])

    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    # plt.legend(loc="lower right", frameon=False, prop={"size":6})

    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    plt.savefig("../../Results/figures/SuppleFig5a.pdf", dpi=400, bbox_inches='tight')


if __name__ == "__main__" :
    # EC_subsystem()
    # main()
    # EC_Kcat()
    # Kcat_subsystem()
    # plot_subsystem_Kcat_counts()
    plot_subsystem_distribution()

    # Results :
    # The median value in log10 of Primary-CE is 1.01
    # The median value in log10 of Intermediate is 0.73
    # The median value in log10 of Primary-AFN is 0.86
    # The median value in log10 of Secondary is 0.86

    # The median value of Primary-CE is 10.26
    # The median value of Intermediate is 5.52
    # The median value of Primary-AFN is 7.23
    # The median value of Secondary is 7.18

