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
from SOAPpy import SOAPProxy ## for usage without WSDL file
endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
client      = SOAPProxy(endpointURL)
password    = hashlib.sha256(password).hexdigest()
credentials = email + ',' + password


# parameters = "j.doe@example.edu,"+password+","+"ecNumber*1.1.1.1#organism*Homo sapiens"
# resultString = client.getSequence(parameters)

# parameters = credentials+","+"ecNumber*1.1.1.1#organism*Homo sapiens"
parameters = credentials+","+"ecNumber*4.1.1.85#organism*Escherichia coli K-12"
# parameters = credentials+","+"ecNumber*3.1.22.4#organism*Escherichia coli"
# parameters = credentials+","+"ecNumber*3.1.3.17#organism*Oryctolagus cuniculus"
sequence = client.getSequence(parameters)

# sequence = client.getSequence("ecNumber*1.1.1.1#organism*Mus musculus")

print(sequence)