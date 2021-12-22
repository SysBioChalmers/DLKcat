#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import os
import math
import json
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import rc
import seaborn as sns
import pandas as pd
from scipy.stats import gaussian_kde
from scipy import stats
from sklearn.metrics import mean_squared_error,r2_score


def main() :
    experimental_values = list()
    predicted_values = list()
    with open('../../Data/test_dataset/test_out.txt', 'r') as testfile :
        testData = testfile.readlines()[1:]

    number = 0
    for data in testData :
        line = data.strip().split('\t')
        # print(line)

        enzymeType = line[2]
        if enzymeType == 'wildtype' :
            number += 1
            experimental, predicted = float(line[0]), float(line[1])
            experimental_values.append(experimental)
            predicted_values.append(predicted)

    # correlation, p_value = stats.pearsonr(x, y)
    correlation, p_value = stats.pearsonr(experimental_values, predicted_values)

    # https://blog.csdn.net/u012735708/article/details/84337262?utm_medium=distribute.pc_relevant.none-
    # task-blog-BlogCommendFromMachineLearnPai2-1.pc_relevant_is_cache&depth_1-utm_source=
    # distribute.pc_relevant.none-task-blog-BlogCommendFromMachineLearnPai2-1.pc_relevant_is_cache
    r2 = r2_score(experimental_values,predicted_values)
    rmse = np.sqrt(mean_squared_error(experimental_values,predicted_values))

    print('The data point number: %s' % number)
    print('r is: %.2f' % correlation)
    # print('p value is: %s' % p_value)
    print('p value is: %.4f' % p_value)
    print('R2 is: %.2f' % r2)
    print('RMSE is: %.2f' % rmse)
    print('\n')
        
    # Results:
    # The data point number: 940
    # r is: 0.65
    # p value is: 0.0000
    # R2 is: 0.40
    # RMSE is: 1.11

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

    plt.text(-4.7, 6.9, 'r = %.2f' % correlation, fontweight ="normal", fontsize=6)
    plt.text(-4.7, 5.9, 'p value = 0', fontweight ="normal", fontsize=6)
    plt.text(3, -5.0, 'Wildtype', fontweight ="normal", fontsize=6)
    plt.text(-4.7, 4.8, 'N = 940', fontweight ="normal", fontsize=6)

    plt.rcParams['font.family'] = 'Helvetica'

    plt.xlabel("Experimental $k$$_\mathregular{cat}$ value", fontdict={'weight': 'normal', 'fontname': 'Helvetica', 'size': 7}, fontsize=7)
    plt.ylabel('Predicted $k$$_\mathregular{cat}$ value',fontdict={'weight': 'normal', 'fontname': 'Helvetica', 'size': 7},fontsize=7)

    plt.xticks([-6, -4, -2, 0, 2, 4, 6, 8])
    plt.yticks([-6, -4, -2, 0, 2, 4, 6, 8])

    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)

    # plt.rcParams['text.usetex'] = True

    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)

    plt.savefig("../../Results/figures/SuppleFig6b.pdf", dpi=400, bbox_inches='tight')

if __name__ == '__main__' :
    main()

