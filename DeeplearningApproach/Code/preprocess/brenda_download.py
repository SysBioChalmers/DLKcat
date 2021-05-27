#!/usr/bin/python
################################################################################
# createECfiles
# Reads all data in kinetic_data and creates all EC files.
#
# Benjamin Sanchez. Last edited: 2018-04-10
################################################################################

# Updated by:
# Author: LE YUAN
# Date: 2020-07-10

#INPUTS:
#1) Path in which all BRENDA queries are (from script retrieveBRENDA.py):
input_path = '../../Data/database/brenda_ec'
#2) Path in which you wish to store all EC files:
output_path = '../../Data/database/Kcat_brenda'

################################################################################

#Read all BRENDA file names:
import os
prev_path = os.getcwd()
os.chdir(input_path)
dir_files = os.listdir(input_path)
# print(dir_files)
dir_files.sort()
#Main loop: Adds each BRENDA file's info to the corresponding EC file.
previous  = ''
for i in dir_files:
    #Define EC number and variable name:
    sep_pos   = i.find('_')
    ec_number = i[0:sep_pos]
    var_name  = i[sep_pos+1:len(i)-4]

    #Read all data in BRENDA file:
    fid       = open(i,'r')
    data      = fid.read()
    fid.close()
    
    #Detect a change of EC number:
    if ec_number != previous:
        if previous != '':
            #Save previous ec_table in a EC file:
            os.chdir(output_path)
            fid = open(previous + '.txt','w')
            for j in ec_table:
                fid.write(j)
            
            fid.close()
            print 'Succesfully constructed ' + previous + ' file.'
            os.chdir(input_path)

        #Reset ec_table (initialize it in the first iteration):    
        ec_table = []
    
    #Define query to find in string "data", according to the variable name:

    if var_name == 'KM':
        variable  = '#kmValue*'

    elif var_name == 'MW':
        variable  = '#molecularWeight*'

    elif var_name == 'PATH':
        variable  = '#pathway*'

    elif var_name == 'SEQ':
        variable  = '#sequence*'

    elif var_name == 'SA':
        variable  = '#specificActivity*'

    elif var_name == 'KCAT':
        variable  = '#turnoverNumber*'

    #Split the string in N+1 parts, where N is the number of values for the
    #given variable:
    options = data.split(variable)
    # options = ['#kmValue*']
    for k in options:
        #Find the end of the value of interest and save it as k_value:
        value_pos = k.find('#')
        if value_pos != -1:
            k_value = k[0:value_pos]

            #If there is a substrate, split will create 2 strings and the info
            #will be at the beginning of string 2. Applies to KM & KCAT.
            k_split = k.split('#substrate*')
            if len(k_split) == 1:
                k_split     = k_split[0]
                k_substrate = '*'

            else:
                k_split     = k_split[1]
                k_substrate = k_split[0:k_split.find('#')]
                if k_substrate == '':
                    k_substrate = '*'          
            
            #If there is a commentary, split will create 2 strings and the info
            #will be at the beginning of string 2. Applies to all except PATH
            #and SEQ.
            k_split2 = k_split.split('#commentary')
            if len(k_split2) == 1:
                k_split2  = k_split2[0]
                k_comment = '*'

            else:
                k_split2  = k_split2[1]
                k_comment = k_split2[0:k_split2.find('#')]
            
            #If there is a organism, split will create 2 strings and the info
            #will be at the beginning of string 2. Applies to all except PATH.
            k_split = k.split('#organism*')
            if len(k_split) == 1:
                k_org = '*'

            else:
                k_split = k_split[1]
                k_org   = k_split[0:k_split.find('#')]
                if k_org == '':
                    k_org = '*'
            
            #Append data to ec_table in the following format:
            #[variable   organism   value]
            ec_table.append(var_name + '\t' + k_org + '\t' + k_value + '\t')
            #[substrate(if any, otherwise '*') commentary(if any, otherwise '*')
            ec_table.append(k_substrate + '\t' + k_comment + '\n')

    #Update previous ec number:
    previous = ec_number

#Write last EC file:
os.chdir(output_path)
fid = open(previous + '.txt','w')
for j in ec_table:
    fid.write(j)
            
fid.close()
print 'Succesfully constructed ' + previous + ' file.'
os.chdir(prev_path)

################################################################################
