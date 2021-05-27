#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN

import csv
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc


with open('../../Data/database/Kcat_combination_0918.json', 'r') as infile :
    entries = json.load(infile)

print(len(entries))

Kcat = [float(entry['Value']) for entry in entries]

plt.figure(figsize=(3,3))

# To solve the 'Helvetica' font cannot be used in PDF file
# https://stackoverflow.com/questions/59845568/the-pdf-backend-does-not-currently-support-the-selected-font
rc('font',**{'family':'serif','serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42

# plt.axes([0.12,0.12,0.83,0.83])

plt.tick_params(direction='in')
plt.tick_params(which='major',length=1.5)
plt.tick_params(which='major',width=0.4)

plt.hist(Kcat,5000,color='#2166ac')
plt.xlabel('$k$$_\mathregular{cat}$ value', fontsize=7)
plt.ylabel('Counts', fontsize=7)

plt.rcParams['font.family'] = 'Helvetica'

# plt.xlim(0,500000)
# plt.xticks([0,10,100,1000,10000,100000])

ax = plt.gca()
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)

plt.yscale('log')
plt.xscale('log')

plt.xticks(fontsize=6)
plt.yticks(fontsize=6)

plt.tight_layout()

plt.savefig("../../Results/figures/SuppleFig1a.pdf", dpi=400)



