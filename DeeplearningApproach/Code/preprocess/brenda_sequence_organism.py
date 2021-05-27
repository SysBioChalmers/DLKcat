#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-06-16


# E-mail in BRENDA:
email = 'youremail'
# Password in BRENDA:
password = 'yourpassword'


# #Construct BRENDA client:
import string
import hashlib
import os
import json
from SOAPpy import SOAPProxy ## for usage without WSDL file
endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
client      = SOAPProxy(endpointURL)
password    = hashlib.sha256(password).hexdigest()
credentials = email + ',' + password


filenames = os.listdir('../../Data/database/Kcat_brenda')
# print(len(filenames)) # 1741 EC files
i = 0

EC_organisms = dict()
for filename in filenames :
    EC = filename[2:-4]
    print(EC)
    if filename != '.DS_Store' :
        with open("../../Data/database/Kcat_brenda/%s" %(filename), 'r') as file :
            lines = file.readlines()
    organisms = list()
    for line in lines[1:] :
        data = line.strip().split('\t')
        organism = data[1]
        organisms.append(organism)
    organisms = list(set(organisms))

    organism_seqcounts = dict()
    for organism in organisms :
        print(organism)

        # parameters = "j.doe@example.edu,"+password+","+"ecNumber*1.1.1.1#organism*Mus musculus"
        # resultString = client.getSequence(parameters)

        # parameters = credentials+","+"ecNumber*1.1.1.1#organism*Homo sapiens"
        # parameters = credentials+","+"ecNumber*3.1.3.17#organism*Oryctolagus cuniculus"
        parameters = credentials+","+"ecNumber*%s#organism*%s" %(EC, organism)
        sequence = client.getSequence(parameters)
        split_sequences = sequence.strip().split('#!') #noOfAminoAcids #!
        # sequence = client.getSequence("ecNumber*1.1.1.1#organism*Mus musculus")

        # for seq in split_sequences :
        #     print(seq)
        #     print('--------------------------------')
        organism_seqcounts[organism] = len(split_sequences)
        # print(len(split_sequences))

    EC_organisms[EC] = organism_seqcounts

print(EC_organisms)
with open('../../Data/database/brenda_EC_organims_try.json', 'w') as outfile:
    json.dump(EC_organisms, outfile, indent=4)
    