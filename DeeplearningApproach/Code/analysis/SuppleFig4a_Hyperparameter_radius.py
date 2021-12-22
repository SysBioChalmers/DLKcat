#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-11-06
# https://blog.csdn.net/roguesir/article/details/77839721

import matplotlib.pyplot as plt
from matplotlib import rc


with open('../../Data/output_hyper/MAEs--all--radius0--ngram1--dim10--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50.txt', 'r') as infile1 :
    lines1 = infile1.readlines()[1:]

with open('../../Data/output_hyper/MAEs--all--radius1--ngram2--dim10--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50.txt', 'r') as infile2 :
    lines2 = infile2.readlines()[1:]

with open('../../Data/output_hyper/MAEs--all--radius2--ngram3--dim10--layer_gnn3--window11--layer_cnn3--layer_output3--lr1e-3--lr_decay0.5--decay_interval10--weight_decay1e-6--iteration50.txt', 'r') as infile3 :
    lines3 = infile3.readlines()[1:]

epoch_1 = list()
R2_1 = list()
for line in lines1[:30] :
	data = line.strip().split('\t')
	# print(data)
	epoch_line = int(data[0])
	R2_line = float(data[-2])
	if epoch_line%2 == 0 or epoch_line in [1,30] :
		epoch_1.append(epoch_line)
		R2_1.append(R2_line)

epoch_2 = list()
R2_2 = list()
for line in lines2[:30] :
	data = line.strip().split('\t')
	# print(data)
	epoch_line = int(data[0])
	R2_line = float(data[-2])
	if epoch_line%2 == 0 or epoch_line in [1,30] :
		epoch_2.append(epoch_line)
		R2_2.append(R2_line)

epoch_3 = list()
R2_3 = list()
for line in lines3[:30] :
	data = line.strip().split('\t')
	# print(data)
	epoch_line = int(data[0])
	R2_line = float(data[-2])
	if epoch_line%2 == 0 or epoch_line in [1,30] :
		epoch_3.append(epoch_line)
		R2_3.append(R2_line)

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

plt.plot(epoch_1,R2_1,color='#FC9E05',linestyle='dashed',linewidth=0.75,marker='s',markerfacecolor='#FC9E05', markersize=1,label='0-radius & 1-gram')
plt.plot(epoch_2,R2_2,color='#2166ac',linestyle='dashed',linewidth=0.75,marker='^',markerfacecolor='#2166ac', markersize=1,label='1-radius & 2-gram')
plt.plot(epoch_3,R2_3,color='#b2182b',linestyle='dashed',linewidth=0.75,marker='o',markerfacecolor='#b2182b', markersize=1,label='2-radius & 3-gram')

plt.rcParams['font.family'] = 'Helvetica'

plt.xticks([0,5,10,15,20,25,30])
plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])
# plt.yticks([0,0.2,0.4,0.6,0.8])

plt.xlabel('Epoch', fontsize=7)
# plt.ylabel('R2', fontsize=7)
plt.ylabel('R$^2$', fontsize=7)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.legend(frameon=False, prop={"size":6})

ax = plt.gca()
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.spines['top'].set_linewidth(0.5)
ax.spines['right'].set_linewidth(0.5)

plt.savefig("../../Results/figures/SuppleFig4a.pdf", dpi=400, bbox_inches='tight')


