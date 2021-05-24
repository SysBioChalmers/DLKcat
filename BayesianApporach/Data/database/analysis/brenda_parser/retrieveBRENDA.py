#!/usr/bin/python
################################################################################
# retrieveBRENDA
# Acces the web client and retrieves all EC data from BRENDA. Creates files with
# BRENDA output for all organisms and EC numbers for which there is data.
#
# Benjamin Sanchez. Last edited: 2018-04-10
################################################################################

#INPUTS:
#1) Path in which you wish to store all BRENDA queries:
output_path = '/Users/.../brenda_parser/raw_data'
#2) Last field processed (if the program was interrupted), e.g. 'KM'. If you
#   want to start from scratch, leave empty:
last_field = ''
#3) Last EC number processed (if the program was interrupted), e.g. '1.2.3.4'.
#   If you want to start from scratch, leave empty:
last_EC = ''
#4) E-mail in BRENDA:
email = 'xxxxx'
#5) Password in BRENDA:
password = 'xxxxx'

################################################################################

#extract_field: Function that extracts all BRENDA data from a specific field.
def extract_field(field,last):

    #Construct list of EC numbers, based on the enzymes for which there is
    #data on BRENDA:
    if field == 'KM':
        ECstring = client.getEcNumbersFromKmValue(credentials)

    elif field == 'MW':
        ECstring = client.getEcNumbersFromMolecularWeight(credentials)

    elif field == 'PATH':
        ECstring = client.getEcNumbersFromPathway(credentials)

    elif field == 'SEQ':
        ECstring = client.getEcNumbersFromSequence(credentials)

    elif field == 'SA':
        ECstring = client.getEcNumbersFromSpecificActivity(credentials)

    elif field == 'KCAT':
        ECstring = client.getEcNumbersFromTurnoverNumber(credentials)
  
    EClist = ECstring.split('!')

    #Loop that retrieves data from BRENDA and saves it in txt files. Starts
    #from the last EC number queried:
    start = 0
    for ECnumber in EClist:
    
        #Detects the starting point (the last EC number queried):
        if not start and (ECnumber == last or last == ''):
            start = 1
        
        if start:
            #The code will retrieve data for all organisms:
            query  = credentials + ',ecNumber*' + ECnumber + '#organism*'
            succes = 0
            
            #The try/except block inside the while is to avoid timeout PROXY
            #and encoding errors:
            while succes < 10:
                try:
                    file_name = 'EC' + ECnumber + '_' + field
                    print file_name

                    if field == 'KM':
                        data = client.getKmValue(query)

                    elif field == 'MW':
                        data = client.getMolecularWeight(query)

                    elif field == 'PATH':
                        data = client.getPathway(query)

                    elif field == 'SEQ':
                        data = client.getSequence(query)

                    elif field == 'SA':
                        data = client.getSpecificActivity(query)

                    elif field == 'KCAT':
                        data = client.getTurnoverNumber(query)
                    
                    #Once the querie was performed succesfully, the data is
                    #copied in txt files:    
                    if data:
                        fid = open(file_name + '.txt','w')
                        fid.write(data.decode('ascii','ignore'))
                        fid.close()

                    succes = 10
                
                except:
                    #Let the server cool of for a bit. If after 10 times it
                    #still fails, the query is discarded:
                    time.sleep(1)
                    succes += 1

################################################################################

#Main script
                    
#Change path:
import os
prev_path = os.getcwd()
os.chdir(output_path)

#Construct BRENDA client:
import string
import hashlib
from SOAPpy import SOAPProxy ## for usage without WSDL file
endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
client      = SOAPProxy(endpointURL)
password    = hashlib.sha256(password).hexdigest()
credentials = email + ',' + password

#Information to retrieve: km, M.W., pathway, sequence, specific activity and kcat.
fields = ['KM','MW','PATH','SA','KCAT']
import time

#Loop that retrieves all fields. Starts by the last one queried:
start = 0
for field in fields:
    if not start and (field == last_field or last_field == ''):
        start = 1

    if start:
        extract_field(field,last_EC)

os.chdir(prev_path)

################################################################################
