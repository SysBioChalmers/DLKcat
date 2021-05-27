#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-06-25
# This script is to extract Kcat data from all EC files into one file.

import os
import csv
import re

outfile = open("../../Data/database/Kcat_brenda.tsv", "wt")
# with open("./Kcat_sabio.tsv", "wt") as outfile :
tsv_writer = csv.writer(outfile, delimiter="\t")
tsv_writer.writerow(["EntryID", "Type", "ECNumber", "Substrate", 'EnzymeType', "Organism","Value", "Unit"])

filenames = os.listdir('../../Data/database/Kcat_brenda')
# print(len(filenames)) # 3339 EC files

i = 0
j = 0
for filename in filenames :
    print(filename[2:-4])
    if filename != '.DS_Store' :
        with open("../../Data/database/Kcat_brenda/%s" %(filename), 'r', encoding="utf-8") as file :
            lines = file.readlines()

    for line in lines[1:] :
        data = line.strip().split('\t')
        value = float(data[2])
        desc = data[4]
        # print(desc)
        if value > 0 :  # Kcat value should not be less than 0, but there exist some weird values, less than 0.
            # print(value)
            # print(type(value))
            i += 1
            if 'mutant' in desc or 'mutated' in desc:
                # print(desc)
                mutant = re.findall('[A-Z]\d+[A-Z]', desc)  # re is of great use
                # print(mutant)
                if len(mutant) >=1 :
                    enzymeType = '/'.join(mutant)
                else :
                    continue
                    # print(desc)
                    # j += 1
                    # print(j)
                    # enzymeType = 'wildtype'
                    # print('------------------------------------------------')
            else :
                enzymeType = 'wildtype'
            tsv_writer.writerow([i, 'kcat', filename[2:-4], data[3], enzymeType, data[1], str(value), 's^(-1)'])

outfile.close()
