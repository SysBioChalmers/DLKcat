#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-06-25

# This python script is to obtain protein sequence by uniprot protein id

from urllib import request


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
    # return IdSeq[id]
    
def main() :
    uniprot_sequence('P49384')


if __name__ == "__main__" :
    main()


