#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-07-13


import os
import csv


# with open("./Kcat_sabio_4_new/%s" %('1.1.1.184.txt'), 'r', encoding="utf-8") as file :
#     lines = file.readlines()

# for line in lines[1:] :
#     data = line.strip().split('\t')
#     print(data)


outfile = open("../../Data/database/Kcat_sabio_4_unisubstrate.tsv", "wt")
# with open("./Kcat_sabio.tsv", "wt") as outfile :
tsv_writer = csv.writer(outfile, delimiter="\t")
tsv_writer.writerow(["EntryID", "Type", "ECNumber", "Substrate", "EnzymeType", "PubMedID", 
    "Organism", "UniprotID", "Value", "Unit"])

filenames = os.listdir('../../Data/database/Kcat_sabio_4')
# print(len(filenames)) # 1741 EC files
i = 0
j=0
for filename in filenames :
    print(filename[1:-4])
#     # if filename == '1.1.1.184.txt' :

    if filename != '.DS_Store' :
        with open("../../Data/database/Kcat_sabio_4/%s" % filename, 'r', encoding="utf-8") as file :
            lines = file.readlines()

        for line in lines[1:] :
            data = line.strip().split('\t')
            # print(data)
            try :
                if data[7] == 'kcat' and data[9] :
                    i += 1
                    print(i)
                    print(data)
                    entryID = data[0]
                    for line in lines[1:] :
                        data2 = line.strip().split('\t')
                        if data2[0] == entryID and data2[7] == 'Km' :
                            j += 1
                            print(j)
                            tsv_writer.writerow([j, data[7], data[6], data2[8], data[2], data[3], data[4], data[5], data[9], data[-1]])
            except :
                continue

outfile.close()

