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

def write_Kcat_enzymeType() :
    with open('../../Data/database/Kcat_combination_0918_wildtype_mutant.json', 'r') as infile :
        entries = json.load(infile)

    types_Kcat = {
    'Wildtype': list(),
    'Mutant': list()
    }

    write_items = list()
    for entry in entries :
        enzymeType = entry['Type']
        value = entry['Value']

        if enzymeType == 'wildtype' and float(value) > 0 :
            value_log10 = math.log10(float(value))
            types_Kcat['Wildtype'].append(value_log10)
            write_items.append(['Wildtype', str(value_log10)])
        if enzymeType == 'mutant' and float(value) > 0 :
            value_log10 = math.log10(float(value))
            types_Kcat['Mutant'].append(value_log10)
            write_items.append(['Mutant', str(value_log10)])

    with open('./source_data/SuppleFig6a.txt', 'w') as outfile :
        items = ['Type', 'Kcat value (log10)']
        outfile.write('\t'.join(items)+'\n')

        i = 0
        for item in write_items :
            outfile.write('\t'.join(item)+'\n')

        print("Writing to file done!")


if __name__ == "__main__" :
    write_Kcat_enzymeType()


