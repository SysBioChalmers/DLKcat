#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# https://blog.csdn.net/roguesir/article/details/77839721

import matplotlib.pyplot as plt
from matplotlib import rc

with open('../../Results/output/MAEs--all--radius2--ngram3--dim20--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration30.txt', 'r') as infile :
    lines = infile.readlines()[1:]

# print(len(lines))

epoch_dev = list()
RMSE_dev = list()

for line in lines[:17] :
	data = line.strip().split('\t')
	epoch_line = int(data[0])
	RMSE_line = float(data[-4])
	if epoch_line%2 == 0 or epoch_line in [1,99] :
		epoch_dev.append(epoch_line)
		RMSE_dev.append(RMSE_line)

epoch_test = list()
RMSE_test = list()

for line in lines[:17] :
	data = line.strip().split('\t')
	epoch_line = int(data[0])
	RMSE_line = float(data[-3])
	if epoch_line%2 == 0 or epoch_line in [1,99] :
		epoch_test.append(epoch_line)
		RMSE_test.append(RMSE_line)

# fig=plt.figure(figsize=(1.5,1.5))
# # fig.add_axes([0.2,0.2,0.6,0.6])
# # fig.add_axes([6.8/39.6,6.8/39.6,31.7/39.6,31.7/39.6])
# fig.add_axes([0.12,0.12,0.83,0.83])

plt.figure(figsize=(1.5,1.5))

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

plt.plot(epoch_dev,RMSE_dev,color='#b2182b',linestyle='dashed',linewidth=0.75,marker='o',markerfacecolor='#b2182b', markersize=3,label='Validation')
plt.plot(epoch_test,RMSE_test,color='#2166ac',linestyle='dashed',linewidth=0.75,marker='^',markerfacecolor='#2166ac', markersize=3,label='Test') 

plt.rcParams['font.family'] = 'Helvetica'
# plt.rc('font', family='Helvetica')
plt.xticks([0,3,6,9,12,15,18])
plt.yticks([1.00,1.05,1.10,1.15,1.20])

plt.xlabel('Epoch', fontsize=7)
plt.ylabel('RMSE', fontsize=7)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(frameon=False, prop={"size":6})

ax = plt.gca()
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)

plt.savefig("../../Results/figures/Fig2a.pdf", dpi=400, bbox_inches='tight')

