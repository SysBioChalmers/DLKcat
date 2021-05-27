#!/usr/bin/python
################################################################################
# findMaxKcats
# Reads all EC files and finds the max value for each substrate for the chosen
# microorganism on the different enzymatic parameters [Kcat, KM, SA, MW]. 
# For each parameter Writes a table with the following columns:
#   * EC number
#   * substrate
#   * organism name//taxonomical classification (according to KEGG)
#   * Max value
#   * metabolic pathways

# Benjamin Sanchez. Last edited: 2015-08-26
# Ivan Domenzain.   Last edited: 2018-04-10
################################################################################
#INPUTS:
#1) Enzymatic parameters
features_list = ['KCAT','SA', 'MW']
#2) Path in which the EC files are stored (from script createECfiles.py):
input_path = '/Users/.../brenda_parser/EC_files'
#3) Path in which you wish to store the final table:
output_path = '/Users/.../brenda_parser/max_data'
################################################################################

#sub_max_std: Recieves a list of substrates///organism_info///values, returns
#3 lists: substrates - max - std

def sub_max_std(data):
    #Sorts list, add a last empty line and initialize variables:
    data.sort()
    data.append('')
    #for every substrate gets the index of all its appearences in the rows
    #of the EC data
    substrates, org_strings_reps, values_reps, reps_indexes = substrate_repetitions(data)
    org_strings, values =\
    find_in_substrate(substrates, org_strings_reps, values_reps, reps_indexes)
    #get maximum kvalues
    values = maximum_values(values)

    return (substrates,org_strings,values)


################################################################################
#maximum_values: Gets the maximum Kvalue for each organism related
#to each substrate
def maximum_values(values):
    for subs_values in values:
        i = values.index(subs_values)
        for org_values in subs_values:
            j = subs_values.index(org_values)
            #org_values = max(org_values)
            values[i][j] = max(values[i][j])

    return(values)
################################################################################
#find_in_substrate: Finds the organisms and K values related to each substrate
def find_in_substrate(substrates, org_strings_reps, values_reps, reps_indexes):

    org_strings = []
    values      = []
    for i in range(len(substrates)):
        subs_orgs   = []
        subs_values = []
        
        for j in reps_indexes[i]:
            
            try:
                org_index = subs_orgs.index(org_strings_reps[j])
                subs_values[org_index].append(values_reps[j])
            
            except:
                subs_orgs.append(org_strings_reps[j])
                org_index = subs_orgs.index(org_strings_reps[j])
                subs_values.append([values_reps[j]])

        org_strings.append(subs_orgs)
        values.append(subs_values)

    return(org_strings, values)

################################################################################
#substrate_repetitions: Finds from each EC file the organism and K values related
#to each substrate

def substrate_repetitions(data):
    
    substrates       = []
    org_strings_reps = []
    values_reps      = []
    reps_indexes = []

    for row in data:
        if row != '':
            row_index = data.index(row)
            #gets index if row substrate is repeated
            try:
                subs_indx = substrates.index(row[0:row.find('///')])
            #if new substrate
            except:
                substrates.append(row[0:row.find('///')])
                subs_indx = substrates.index(row[0:row.find('///')])
                reps_indexes.append([])
            #list with the indexes of the rows with repeated substrates
            reps_indexes[subs_indx].append(row_index)
            #organisms for related to each substrate
            org_strings_reps.append(row[row.find('///')+3:row.find('////')])
            #values found for each organism
            values_reps.append(float(row[row.find('////')+4:len(row)]))

    return(substrates, org_strings_reps, values_reps, reps_indexes)

################################################################################

#create_orgs_list: Finds all the organism names for which data is available in
#BRENDA database. As an output a list with all found names is created.

def brenda_orgs_list(dir):
    brenda_orgs=[]
    for ec in dir:

        fid     = open(ec,'r')
        csv_fid = csv.reader(fid,delimiter='\t')
        
        try:
            for row in csv_fid:
                if row != '' and row[0] != 'SEQ' and row[0] != '*':
                
        #Uncomment and indent properly if you want to exclude any name longer
        #two words (mutants for example, but not exclusively)
        
                    #second_blank = row[1].find(' ',row[1].find(' ')+1)

                    #if second_blank == -1:
                    org_name = row[1].lower()
                    #else:
                    #    org_name=row[1][0:second_blank]
                
                    if brenda_orgs.count(org_name)==0:
                        brenda_orgs.append(org_name)
        except:
            pass#brenda_orgs.append(org_name)
    return (brenda_orgs)


################################################################################

#KEGG_orgs_list: Creates a list with all the organisms available at KEGG, as an
#output  a table with the fields: name, KEGG code and Taxonomy is created.
def KEGG_orgs_list():

    #URL that returns available data of the gene entry on KEGG
    url = 'http://rest.kegg.jp/list/organism'
    #Try/except for avoiding timeout exceedings
    try:
        query = urllib2.urlopen(url, timeout=20).read()
    except:
        query=''

    entries   = query.split('\n')
    KEGG_list = []
    tax_kegg  = []
    codes     = []
    for row in entries:
        if row != '':
            row_list = row.split('\t')
            if len(row_list)>1:
                row_list=[row_list[2],row_list[3],row_list[1]]
                if row_list[0].find('(') != -1:
                    row_list[0]=row_list[0][0:row_list[0].find('(')-1]
       
                #Saves only organisms with specific taxonomic classifications
                taxonomy=row_list[1].lower()
                if taxonomy.find('eukaryotes')!= -1 or taxonomy.find('prokaryotes')!= -1:
                    KEGG_list.append(row_list[0].lower())
                    tax_kegg.append(taxonomy)
                    codes.append(row_list[2])
    return(KEGG_list, tax_kegg, codes)

################################################################################

#orgs_list: Two possible options 1) Merges BRENDA and KEGG organisms lists
# creates a list with only coincidences between lists.

def orgs_list(dir):

    KEGG_orgs, info_KEGG, codes = KEGG_orgs_list()

    brenda_orgs         = brenda_orgs_list(dir)
    #print brenda_orgs
    organism_list = []
    taxonomy      = []
    org_codes     = []
    #i=0
    counter=0
    
    for B_org in brenda_orgs:
        flag = False
        if B_org != '*':
            i=0
        
            while (i < len(KEGG_orgs) and flag == False):
                K_org = KEGG_orgs[i]

                if K_org.find(B_org)!= -1:
                
                    flag    = True
                    counter = counter+1
                    organism_list.append(B_org)
                    taxonomy.append(info_KEGG[i])
                    org_codes.append(codes[i])
                i = i +1
    
            if flag == False:
                organism_list.append(B_org)
                taxonomy.append('*')
                org_codes.append('*')

        #for B_org in brenda_orgs:
        #    if KEGG_orgs.count(B_org) != 0:
        #        counter=counter+1
        #        i=KEGG_orgs.index(B_org)
        #        organism_list.append(B_org)
        #        taxonomy.append(info_KEGG[i])
        #        org_codes.append(codes[i])
        #    else:
        #        organism_list.append(B_org)
        #        taxonomy.append('*')
        #        org_codes.append('*')

    return(organism_list,taxonomy,org_codes)

################################################################################

#EC_string: Receives the information in the EC file and builds a string with
#substrates, related organisms and Kvalues

def EC_string(csv_fid, feature_name):
    data_string = []
    ec_pathways = ''
    for row in csv_fid:
        if row[0] != '':
            row[4] = row[4].lower()
            mutant = max(row[4].find('mutant'),row[4].find('mutated'))
            #Ignore invalid values:
            if row[2] != '-999' and mutant == -1:
                #Only allow Kcats <= 1e7 [Bar-Even et al. 2011]
                if row[0] == feature_name and float(row[2]) <= 1e7:
                    #Looks for the organism in the organism merged list
                    #in order to include taxonomical info if available
                    try:
                        org_index  = organism_list.index(row[1].lower())
                        org_string = organism_list[org_index]+'//'+\
                                     taxonomy[org_index]+ '//'+ organism_code[org_index]
                            
                        data_string.append(row[3].lower() + '///' +\
                                            org_string + '////' + row[2])
                    except:
                        print 'Organism not found in KEGG or BRENDA'
        
            #Gets the associated not engineered pathways to the
            #EC number if present
            if row[0] == 'PATH' and row[2].lower() != 'metabolic pathways':
                if row[2].find('(engineered)') == -1:
                    ec_pathways = ec_pathways + row[2].lower() + '///'
    #If path not found an asterisk is added to the field
    if ec_pathways == '' or ec_pathways == ' ' or ec_pathways == '\0':
        ec_pathways = '*'
    if len(ec_pathways) > 3 and ec_pathways[-3] == '/':
        ec_pathways = ec_pathways[:-3]

    return(data_string, ec_pathways)

################################################################################

#Main Script

#Read all EC file names:
import os
prev_path = os.getcwd()
os.chdir(input_path)
dir_files = os.listdir(input_path)
dir_files.sort()
import urllib2
import csv
organism_list,taxonomy,organism_code = orgs_list(dir_files)

for feature_name in features_list:
    #Main loop:
    output = ''
    for ec in dir_files:
        ec_number  = ec[0:len(ec)-4]
        fid       = open(ec,'r')
        csv_fid   = csv.reader(fid,delimiter='\t')
        #Builds a string with all the information in the EC file
        data_string, ec_pathways = EC_string(csv_fid, feature_name)
        fid.close()

        import numpy
        substrates,org_strings,max_values = sub_max_std(data_string)
    
        for sub in substrates:
            i = substrates.index(sub)
            for org in org_strings[i]:
                j = org_strings[i].index(org)
                output = output+ec_number+'\t'+sub+'\t'+org+'\t'+str(max_values[i][j])+'\t'+ ec_pathways+'\n'

    
        print 'Processed file ' + ec + ' ' + feature_name
    #Write output:
    os.chdir(output_path)
    fid  = open('max_' + feature_name + '.txt','w')
    fid.write(output)
    fid.close()
    os.chdir(input_path)
os.chdir(prev_path)

################################################################################
