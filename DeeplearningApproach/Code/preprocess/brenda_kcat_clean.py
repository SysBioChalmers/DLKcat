#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-07-09  Run in python 3.7
# This script is to clean Kcat data extracted from BRENDA database

import csv

with open("../../Data/database/Kcat_brenda.tsv", "r", encoding='utf-8') as file :
    lines = file.readlines()[1:]

Kcat_data = list()
Kcat_data_include_value = list()
for line in lines :
    # print(line)
    data = line.strip().split('\t')
    Type = data[1]
    ECNumber = data[2]
    Substrate = data[3]
    EnzymeType = data[4]
    Organism =data[5]
    Value = data[6]
    Unit = data[7]
    Kcat_data_include_value.append([Type, ECNumber, Substrate, EnzymeType, Organism, Value, Unit])
    Kcat_data.append([Type, ECNumber, Substrate, EnzymeType, Organism])

print(len(Kcat_data))  # 69140, in which 22723 mutant 46417 wildtype

new_lines = list()
for line in Kcat_data :
    if line not in new_lines :
        new_lines.append(line)

print(len(new_lines))  # 67566 included all elements, 52390 included all except for Kcat value and unit, 32305 if further not include enzymeType

i = 0
clean_Kcat = list()
for new_line in new_lines :
    # print(new_line)
    i += 1
    print(i)
    value_unit = dict()
    Kcat_values = list()
    for line in Kcat_data_include_value :
        if line[:-2] == new_line :
            value = line[-2]
            value_unit[str(float(value))] = line[-1]
            # print(type(value))  # <class 'str'>
            Kcat_values.append(float(value))
    # print(value_unit)
    # print(Kcat_values)
    max_value = max(Kcat_values) # choose the maximum one for duplication Kcat value under the same entry as the data what we use
    unit = value_unit[str(max_value)]
    # print(max_value)
    # print(unit)

    new_line.append(str(max_value))
    new_line.append(unit)
    if new_line[-1] == 's^(-1)' :
        clean_Kcat.append(new_line)

# print(clean_Kcat)
print(len(clean_Kcat))  # 52390


with open("../../Data/database/Kcat_brenda_clean.tsv", "w") as outfile :
    records = ['Type', 'ECNumber', 'Substrate', 'EnzymeType', 'Organism', 'Value', 'Unit']
    outfile.write('\t'.join(records) + '\n')
    for line in clean_Kcat :
        outfile.write('\t'.join(line) + '\n')

