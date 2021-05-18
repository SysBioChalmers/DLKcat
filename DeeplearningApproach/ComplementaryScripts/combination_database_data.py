#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-08-08

# This python script is to obtain protein sequence for each Kcat entries

import os
import re
import json
import requests
import time
from urllib import request
from zeep import Client
import hashlib
# import string
# import hashlib
# from SOAPpy import WSDL
# from SOAPpy import SOAPProxy ## for usage without WSDL file


# This function is to obtain the protein sequence according to the protein id from Uniprot API
# https://www.uniprot.org/uniprot/A0A1D8PIP5.fasta
# https://www.uniprot.org/help/api_idmapping
def uniprot_sequence(id) :
    url = "https://www.uniprot.org/uniprot/%s.fasta" % id
    IdSeq = dict()

    try :
        data = request.urlopen(url)
        respdata = data.read().decode("utf-8").strip()
        IdSeq[id] =  "".join(respdata.split("\n")[1:])
    except :
        print(id, "can not find from uniprot!")
        IdSeq[id] = None
    print(IdSeq[id])
    return IdSeq[id]
    
def uniprotID_entry() :
    # uniprot_sequence('P18314')
    with open("../complementaryData/Kcat_combination_0731.tsv", "r", encoding='utf-8') as file :
        combination_lines = file.readlines()[1:]

    uniprotID_list = list()
    uniprotID_seq = dict()
    uniprotID_noseq = list()

    i=0
    for line in combination_lines :
        data = line.strip().split('\t')
        uniprotID = data[5]

        if uniprotID :
        #     seq = uniprot_sequence('P49384')
            if ' ' in uniprotID :
                # i += 1  # 561
                # print(i)
                # print(uniprotID.split(' '))
                uniprotID_list += uniprotID.split(' ')
            else :
                # print(uniprotID)
                uniprotID_list.append(uniprotID)

    # print(len(uniprotID_list))    # 14045
    uniprotID_unique = list(set(uniprotID_list))
    # print(len(uniprotID_unique)) # 1776
    # print(uniprotID_unique[-6:])

    for uniprotID in uniprotID_unique :
        i += 1
        print(i)
        sequence = uniprot_sequence(uniprotID)
        if sequence :
            uniprotID_seq[uniprotID] = sequence
        else :
            uniprotID_noseq.append(uniprotID)


    print(len(uniprotID_seq))  # 1755
    print(len(uniprotID_noseq))  # 21
    print(uniprotID_noseq)
    # ['P0A5R0', 'P0C5C1', 'P51698', 'P96807', 'Q01745', 'P00892', 'D0B556', 'V5MWQ6', 'Q02469', 'P96223', 'P0A4Z2', 
    # 'P0A4X4', 'P96420', 'Q47741', 'O05783', 'A3S939', 'P0A4X6', 'P56967', 'O60344', 'P04804', 'O52310']
    # check one by one

    with open('../complementaryData/uniprotID_entry.json', 'w') as outfile :
        json.dump(uniprotID_seq, outfile, indent=4)

def uniprotID_noseq() :
    with open('../complementaryData/uniprotID_entry.json', 'r') as infile :
        uniprotID_seq = json.load(infile)

    print(len(uniprotID_seq))
    # uniprotID_noseq = ['P0A5R0', 'P0C5C1', 'P51698', 'P96807', 'Q01745', 'P00892', 'D0B556', 'V5MWQ6', 'Q02469', 'P96223', 'P0A4Z2', 
    # 'P0A4X4', 'P96420', 'Q47741', 'O05783', 'A3S939', 'P0A4X6', 'P56967', 'O60344', 'P04804', 'O52310']

    uniprotID_noseq = {'P0A5R0':'P9WIL4', 'P0C5C1':'P9WKD2', 'P51698':'A0A1L5BTC1', 'P96807':'P9WNP2', 'Q01745':'I1S2N3', 'P00892':'P0DP89', 
    'Q02469':'P0C278', 'P96223':'P9WNF8', 'P0A4Z2':'P9WPY2', 'P0A4X4':'P9WQ86', 'P96420':'P9WQB2', 'Q47741':'F2MMN9', 'O05783':'P9WIQ2', 
    'P0A4X6':'P9WQ80', 'P56967':'F2MMP0', 'O60344':'P0DPD6', 'P04804':'P60906', 'O52310':'P0CL72'}
    # 'D0B556', 'A3S939', 'V5MWQ6'  On April 1, 2015 this entry was made redundant.

    for uniprotID, mappedID in uniprotID_noseq.items() :
        sequence = uniprot_sequence(mappedID)
        print(uniprotID)
        print(sequence)
        if sequence :
            uniprotID_seq[uniprotID] = sequence
        else :
            print('No sequence found!---------------------------')

    print(len(uniprotID_seq))  # 1773   'D0B556', 'A3S939', 'V5MWQ6' no sequence found!

    with open('../complementaryData/uniprotID_entry_all.json', 'w') as outfile :
        json.dump(uniprotID_seq, outfile, indent=4)

# You can try to retrieve sequences from uniprot using rest interface.
# Example: (ec: 1.1.1.1 , organisms: Homo sapiens)
# http://www.uniprot.org/uniprot/?query=ec:1.1.1.1+AND+organism:"Homo sapiens"&format=fasta
# full information abut syntax you can find here: http://www.uniprot.org/help/programmatic_access
def seq_by_ec_organism(ec, organism) :
    IdSeq = dict()
    # https://www.biostars.org/p/356687/
    params = {"query": "ec:%s AND organism:%s AND reviewed:yes" % (ec, organism), "format": "fasta"}
    response = requests.get("http://www.uniprot.org/uniprot/", params=params)
    # print(type(response.text)) # <class 'str'>

    try :
        # respdata = response.text.strip()
        # # print(respdata) 
        # IdSeq[ec+'&'+organism] =  "".join(respdata.split("\n")[1:])

        respdata = response.text
        # print(respdata)
        sequence = list()
        seq = dict()
        i = 0
        for line in respdata.split('\n') :
            if line.startswith('>') :
                name=line
                seq[name] = ''
            else :
                seq[name] += line.replace('\n', '').strip()
        IdSeq[ec+'&'+organism] =  list(seq.values())

    except :
        print(ec+'&'+organism, "can not find from uniprot!")
        IdSeq[ec+'&'+organism] = None

    print(IdSeq[ec+'&'+organism])
    return IdSeq[ec+'&'+organism]

# Run in python 2.7
def seq_by_brenda(ec, organism) :
    # # E-mail in BRENDA:
    # email = 'leyu@chalmers.se'
    # # Password in BRENDA:
    # password = 'yuanle13579'

    # endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
    # client      = SOAPProxy(endpointURL)
    # password    = hashlib.sha256(password).hexdigest()
    # credentials = email + ',' + password

    # parameters = credentials+","+"ecNumber*%s#organism*%s" %(ec, organism)
    # content = client.getSequence(parameters)

    # # E-mail in BRENDA:
    # email = 'leyu@chalmers.se'
    # # Password in BRENDA:
    # password = 'yuanle13579'

    # wsdl = "https://www.brenda-enzymes.org/soap/brenda.wsdl"
    # client      = WSDL.Proxy(wsdl)
    # password    = hashlib.sha256(password).hexdigest()
    # credentials = email + ',' + password

    # parameters = credentials+","+"ecNumber*%s#organism*%s" %(ec, organism)
    # content = client.getSequence(parameters)

    # split_sequences = content.strip().split('!') #noOfAminoAcids #!
    # # UniProtKB/TrEMBL is a computer-annotated protein sequence database complementing the UniProtKB/Swiss-Prot Protein Knowledgebase.
    # sequences = list()
    # # print(split_sequences)
    # for sequence in split_sequences :
    #     dict_entry = dict()
    #     # print(sequence)
    #     list_one = sequence.split('#')
    #     # print(list_one)
    #     for one in list_one[:-1] :
    #         # print(one)
    #         dict_entry[one.split('*')[0]] = one.split('*')[1]
    #     # try :
    #     #     if dict_entry['source'] == 'Swiss-Prot' :
    #     #         sequences.append(dict_entry['sequence'])
    #     #     else :
    #     #         continue
    #     # except :
    #     #     sequences = None
    #     try :
    #         sequences.append(dict_entry['sequence'])
    #     except :
    #         sequences = None

    # print(sequences)

    #New method using Python 3 because using Python 2 method provided by BRENDA could just run less than 10 hits as above
    # E-mail in BRENDA:
    email = 'leyu@chalmers.se'
    # Password in BRENDA:
    password = 'yuanle13579'

    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    password    = hashlib.sha256(password.encode("utf-8")).hexdigest()
    client = Client(wsdl)
    # credentials = email + ',' + password

    # parameters = credentials+","+"ecNumber*%s#organism*%s" %(ec, organism)
    parameters = ( email,password,"ecNumber*%s" % ec,"organism*%s" % organism, "sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*Swiss-Prot", "id*" ) # *Swiss-Prot
    entries = client.service.getSequence(*parameters)

    # print(entries)

    sequences = list()
    # print(split_sequences)
    if entries :
        for entry in entries :
            sequences.append(entry['sequence'])

    print(sequences)
    print(len(sequences))
    return sequences

def nouniprotID_entry_uniprot() :
    # ec = '1.1.1.206'
    # organism = 'Datura stramonium'
    # seq_by_ec_organism(ec, organism)

    with open("../complementaryData/Kcat_combination_0731.tsv", "r", encoding='utf-8') as file :
        combination_lines = file.readlines()[1:]

    IdSeq = dict()
    entries = list()
    i=0
    for line in combination_lines :
        data = line.strip().split('\t')
        ec = data[0]
        organism = data[2]
        uniprotID = data[5]

        if not uniprotID :
            entries.append((ec,organism))

    # print(len(entries))  # 28104
    entries_unique = set(entries)
    # print(len(entries_unique)) # 7258

    for entry in list(entries_unique) :
        # print(entry)
        ec, organism = entry[0], entry[1]
        i += 1
        print('This is', str(i)+'------------')
        IdSeq[ec+'&'+organism] = seq_by_ec_organism(ec, organism)
    # print(len(IdSeq)
        if i%10 == 0 :
            time.sleep(3)

    with open('../complementaryData/nouniprotID_entry_all.json', 'w') as outfile :
        json.dump(IdSeq, outfile, indent=4)

# Run in python 2.7
def nouniprotID_entry_brenda() :
    with open("../complementaryData/Kcat_combination_0731.tsv", "r") as file :
        combination_lines = file.readlines()[1:]

    IdSeq = dict()
    entries = list()
    i=0
    for line in combination_lines :
        data = line.strip().split('\t')
        ec = data[0]
        organism = data[2]
        uniprotID = data[5]

        if not uniprotID :
            entries.append((ec,organism))

    # print(len(entries))  # 28104
    entries_unique = set(entries)
    # print(len(entries_unique)) # 7258

    for entry in list(entries_unique) :
        # print(entry)
        ec, organism = entry[0], entry[1]
        i += 1
        print('This is', str(i)+'------------')
        # print(ec)
        # print(organism)
        IdSeq[ec+'&'+organism] = seq_by_brenda(ec,organism)

    with open('../complementaryData/nouniprotID_entry_brenda.json', 'w') as outfile :
        json.dump(IdSeq, outfile, indent=4)

def combine_sequence() :
    with open('../complementaryData/uniprotID_entry_all.json', 'r') as file1:
        uniprot_file1 = json.load(file1)

    with open('../complementaryData/nouniprotID_entry_all.json', 'r') as file2:  # By Uniprot API
        nouniprot_file2 = json.load(file2)

    with open('../complementaryData/nouniprotID_entry_brenda.json', 'r') as file3:  # By BRENDA API
        nouniprot_file3 = json.load(file3)

    with open("../complementaryData/Kcat_combination_0731.tsv", "r", encoding='utf-8') as file4 :
        Kcat_lines = file4.readlines()[1:]

    # i = 0
    # for proteinKey, sequence in nouniprot_file2.items() :
    #     if sequence :
    #         if len(sequence) == 1 :  # 1178 BRENDA  1919 Uniprot
    #         # if sequence :  # 1784 BRENDA  3363 Uniprot
    #             i += 1   
    #             print(i)
    # print(len(nouniprot_file3))

    i = 0
    j = 0
    n = 0
    entries = list()
    for line in Kcat_lines :
        data = line.strip().split('\t')
        ECNumber, EnzymeType, Organism, Smiles = data[0], data[1], data[2], data[3]
        Substrate, UniprotID, Value, Unit = data[4], data[5], data[6], data[7]

        RetrievedSeq = ''
        entry = dict()
        # print(UniprotID)
        if UniprotID :
            # print(UniprotID)
            try :  # because a few (maybe four) UniprotIDs have no ID as the key 
                if ' ' not in UniprotID :
                    RetrievedSeq = [uniprot_file1[UniprotID]]
                    # print(RetrievedSeq)
                else :
                    # print(UniprotID)
                    RetrievedSeq1 = [uniprot_file1[UniprotID.split(' ')[0]]]
                    RetrievedSeq2 = [uniprot_file1[UniprotID.split(' ')[1]]]
                    if RetrievedSeq1 == RetrievedSeq2 :
                        RetrievedSeq = RetrievedSeq1
                    # if len(RetrievedSeq) == 1:
                    #     print(RetrievedSeq)
            except :
                continue

        else :
            if nouniprot_file2[ECNumber+'&'+Organism] :
                # print(nouniprot_file2[ECNumber+'&'+Organism])
                if len(nouniprot_file2[ECNumber+'&'+Organism]) == 1 :
                    RetrievedSeq = nouniprot_file2[ECNumber+'&'+Organism]
                    # print(RetrievedSeq)
                else :
                    RetrievedSeq = ''

        # print(RetrievedSeq)
        try:  # local variable 'RetrievedSeq' referenced before assignment
            if len(RetrievedSeq) == 1 and EnzymeType == 'wildtype':  # 21108 for all, 9529 wildtype, 11579 mutant (EnzymeType != 'wildtype')
                sequence = RetrievedSeq
                i += 1
                # print(str(i) + '---------------------------')
                # print(ECNumber)
                # print(Organism)
                # print(sequence)

                entry = {
                    'ECNumber': ECNumber,
                    'Organism': Organism,
                    'Smiles': Smiles,
                    'Substrate': Substrate,
                    'Sequence': sequence[0],
                    'Type': 'wildtype',
                    'Value': Value,
                    'Unit': Unit,
                }

                entries.append(entry)

            if len(RetrievedSeq) == 1 and EnzymeType != 'wildtype':
                sequence = RetrievedSeq[0]

                mutantSites = EnzymeType.split('/')
                # print(mutantSites)

                mutant1_1 = [mutantSite[1:-1] for mutantSite in mutantSites]
                mutant1_2 = [mutantSite for mutantSite in mutantSites]
                mutant1 = [mutant1_1, mutant1_2]
                mutant2 = set(mutant1[0])
                if len(mutant1[0]) != len(mutant2) :
                    print(mutant1)
                    n += 1
                    print(str(n) + '---------------------------')  # some are mapped, some are not mapped. R234G/R234K (60, 43 mapped, 17 not mapped)

                mutatedSeq = sequence
                for mutantSite in mutantSites :
                    # print(mutantSite)
                    # print(mutatedSeq[int(mutantSite[1:-1])-1])
                    # print(mutantSite[0])
                    # print(mutantSite[-1])
                    if mutatedSeq[int(mutantSite[1:-1])-1] == mutantSite[0] :
                        # pass
                        mutatedSeq = list(mutatedSeq)
                        mutatedSeq[int(mutantSite[1:-1])-1] = mutantSite[-1]
                        mutatedSeq = ''.join(mutatedSeq)
                        if not mutatedSeq :
                            print('-------------')
                    else :
                        # n += 1
                        # print(str(n) + '---------------------------')
                        mutatedSeq = ''

                if mutatedSeq :
                    # j += 1
                    # print(str(j) + '---------------------------')          
                    entry = {
                        'ECNumber': ECNumber,
                        'Organism': Organism,
                        'Smiles': Smiles,
                        'Substrate': Substrate,
                        'Sequence': mutatedSeq,
                        'Type': 'mutant',
                        'Value': Value,
                        'Unit': Unit,
                    }

                    entries.append(entry)


        #     if len(RetrievedSeq) == 1 :  # 21108 for all, 9529 wildtype, 11579 mutant (EnzymeType != 'wildtype')
        #         sequence = RetrievedSeq
        #         # i += 1
        #         # print(str(i) + '---------------------------')
        #         # print(ECNumber)
        #         # print(Organism)
        #         # print(sequence)

        #         entry = {
        #             'ECNumber': ECNumber,
        #             'Organism': Organism,
        #             'Smiles': Smiles,
        #             'Substrate': Substrate,
        #             'Sequence': sequence[0],
        #             'Value': Value,
        #             'Unit': Unit,
        #         }

        #         entries.append(entry)

        except:
            continue

    # mutatedSeq.replace([int(mutantSite[1:-1])-1], mutantSite[-1])
    print(i)

    print(len(entries))   # 17010  including 9529 wildtype and 7481 mutant

    # with open('../complementaryData/Kcat_combination_0918.json', 'w') as outfile :
    #     json.dump(entries, outfile, indent=4)
    with open('../Data/Kcat_combination_0918_wildtype_mutant.json', 'w') as outfile :
        json.dump(entries, outfile, indent=4)

def check_substrate_seq() :
    with open('../complementaryData/Kcat_combination_0918.json', 'r') as file :  # 0907 file didn't process the mutant sequence, so the mutant sequence is the original for carious mutants
        datasets = json.load(file)

    substrate = [data['Substrate'].lower() for data in datasets]
    sequence = [data['Sequence'] for data in datasets]
    organism = [data['Organism'].lower() for data in datasets]
    EC_number = [data['ECNumber'] for data in datasets]

    unique_substrate = len(set(substrate))
    unique_sequence = len(set(sequence))
    unique_organism = len(set(organism))
    unique_EC_number = len(set(EC_number))

    print('The number of unique substrate:', unique_substrate)
    print('The number of unique sequence:', unique_sequence)
    print('The number of unique organism:', unique_organism)
    print('The number of unique EC Number:', unique_EC_number)

    # The number of unique substrate: 2706
    # The number of unique sequence: 7857
    # The number of unique organism: 856
    # The number of unique EC Number: 1706

if __name__ == "__main__" :
    combine_sequence()
    check_substrate_seq()


