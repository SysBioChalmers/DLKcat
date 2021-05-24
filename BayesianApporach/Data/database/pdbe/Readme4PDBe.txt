Readme4PDBe
1. Download pdb id list from the website (https://www.ebi.ac.uk/pdbe/) by searching for the organism name "Saccharomyces cerevisiae". Save as a CSV file named as "PDBe_search.csv".

2. Run "pdbe.m" to extract information from raw data, which will generate the file "pdb.mat".

3. Run "pdbe_processing.m" to get the final files, i.e., "Protein_stoichiometry" 