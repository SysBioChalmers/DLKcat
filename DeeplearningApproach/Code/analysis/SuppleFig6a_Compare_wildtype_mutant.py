#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2021-06-14

import os
import json
import math
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import numpy as np
import statsmodels.api as sm


def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2
   
    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

def Kcat_enzymeType() :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        entries = json.load(infile)

    types_Kcat = {
    'Wildtype': list(),
    'Mutant': list()
    }
    for entry in entries :
        enzymeType = entry['Type']
        value = entry['Value']

        if enzymeType == 'wildtype' and float(value) > 0 :
            value_log10 = math.log10(float(value))
            types_Kcat['Wildtype'].append(value_log10)
        if enzymeType == 'mutant' and float(value) > 0 :
            value_log10 = math.log10(float(value))
            types_Kcat['Mutant'].append(value_log10)

    return types_Kcat

# https://my.oschina.net/u/4349448/blog/3448306 python code
def plot_enzymeType_distribution() :
    types_Kcat = Kcat_enzymeType()

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

    types_color = {'Wildtype': '#b2182b', 'Mutant': '#2166ac'}
    for types, Kcat in types_Kcat.items() :
        if types in ['Wildtype', 'Mutant'] :

            # print('The median value of %s is %.2f' %(types, math.pow(10, median(Kcat))))
            print('The median value in log10 of %s is %.2f' %(types, median(Kcat)))

        if types in ['Wildtype', 'Mutant']:
            ecdf = sm.distributions.ECDF(Kcat)

            x = np.linspace(min(Kcat),max(Kcat),50000)  # 10000
            y = ecdf(x)

            # plt.xscale('log')
            # plt.xticks([-2, -1, 0, 1, 2, 3])
            # plt.plot(x,y,linewidth='2',label='Primary-CE')

            plt.plot(x,y,linewidth='1.00',label=types,color=types_color[types])

            # https://blog.csdn.net/weixin_38314865/article/details/115173371?utm_medium=distribute.pc_relevant.none-task-blog-baidujs_baidulandingword-5&spm=1001.2101.3001.4242
            # https://blog.csdn.net/dta0502/article/details/83827345
            plt.axvline(x=median(Kcat),ymin=0,ymax=0.5,linewidth='1.00',linestyle='--',color=types_color[types])

    plt.text(3.0, 0.25, 'Wildtype', fontweight ="normal", fontsize=6, color='#b2182b')
    plt.text(3.0, 0.18, 'Mutant', fontweight ="normal", fontsize=6, color='#2166ac')

    plt.rcParams['font.family'] = 'Helvetica'

    plt.xlabel('Experimental $k$$_\mathregular{cat}$ value', fontsize=7)
    plt.ylabel('Cumulative distribution', fontsize=7)

    plt.xlim([-6,6])
    plt.xticks([-6,-4,-2,0,2,4,6])
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    plt.savefig("../../Results/figures/SuppleFig6a.pdf", dpi=400, bbox_inches='tight')


if __name__ == "__main__" :
    plot_enzymeType_distribution()

    # Results :
    # The median value of Wildtype is 8.60
    # The median value of Mutant is 1.90

    # The median value in log10 of Wildtype is 0.93
    # The median value in log10 of Mutant is 0.28

