
Required software
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Python 3.6 or Python 3.7 (Anaconda installation recommended)
- PyTorch
- scikit-learn
- Biopython
- rdkit
- seaborn
- Matplotlib
- pandas
- scipy
- numpy


Usage
~~~~~

- **For users who want to use the deep learning model for prediction, please run these command lines at the terminal:**

  - (1). Download the DLKcat package
  .. code-block:: linux

         git clone https://github.com/SysBioChalmers/DLKcat

  - (2). Download required Python package
  .. code-block:: linux

         pip install numpy requests torch torchvision rdkit-pypi sklearn

  - (3). Change directory to ``DeeplearningApproach`` under the DLKcat package
  .. code-block:: linux

         cd [directory to DLKcat/DeeplearningApproach, where the DeeplearningApproach directory is]

  - (4). Unzip the ``input.zip`` file under the ``Data`` directory
  .. code-block:: linux

         unzip Data/input.zip

  - (5). Change directory to the ``Code/example`` under the DLKcat package
  .. code-block:: linux

         cd Code/example 

  - (6). Now you can use the trained deep learning model for your prediction via one command line. Here, one input file is needed to be prepared, please check the ``Code/example/input.tsv``. For the input file, protein sequence should be provided, and users also need to provide substrate (compound) name or substrate (compound) SMILES, but substrate SMILES is recommended. If it is difficult to find the substrate SMILES, please provide the substrate name and leave the substrate SMILES blank
  .. code-block:: linux

         python prediction_for_input.py input.tsv

  - Then the prediction results (``output.tsv`` file) will be output under the ``Code/example`` directory

- **For running analysis and regenerating all figures:**
  
  - To regenerate all of the figures, unzip the ``input.zip`` file in ``Data/input.zip`` and run the corresponding figure functions in the ``Code/analysis`` directory


Preprocess
~~~~~

- **For data collection and cleaning from the BRENDA database:**
  
  - run the ``brenda_retrieve.py`` to get access to the web client and retrieve dataset from the BRENDA database
  
  - run the ``brenda_download.py`` to read all data in the retrieved files and output all EC files
  
  - run the ``findMaxKvalues_AllOrgs.py`` to read all EC files and find the max value for each substrate for the chosen microorganism

  - run the ``brenda_kcat_preprocess.py`` to generate Kcat data from all EC files into one file
  
  - run the ``brenda_kcat_clean.py`` to clean the dataset from the BRENDA database

  - run the ``brenda_sequence.py`` to get the protein sequence from BRENDA database by one example 

  - run the ``brenda_sequence_organism.py`` to obtain the protein sequences for all data based on EC number and organism and output into one file for further use
  
  - run the ``brenda_get_smiles.py`` to get canonical SMILES just by substrate name for the BRENDA data using PubChem API
  
- **For data collection and cleaning from the SABIO-RK database:**
  
  - run the ``sabio_download.py`` to get access to the web client and download the dataset from the SABIO-RK database

  - run the ``sabio_kcat_unisubstrate.py`` to read all data from the downloaded files and output into one file for further use
  
  - run the ``sabio_kcat_clean_unisubstrate.py`` to clean the data by unifying all entries

  - run the ``sabio_kcat_clean.py`` to used to clean the data for the SABIO-RK data
  
  - run the ``sabio_kcat_unisubstrate_mutant.py`` to annotate the enzyme type information, i.e., wildtype or mutant

  - run the ``uniprot_sequence.py`` to to obtain protein sequence by uniprot protein id

  - run the ``sabio_get_smiles.py`` to get canonical SMILES just by substrate name for the SABIO-RK data and output one file for use

- **For data combination based on the obtained dataset from the BRENDA and the SABIO-RK database:**
  
  - run the ``combination_brenda_sabio.py`` to preliminarily combine the Kcat data from the BRENDA and the SABIO-RK database
  
  - run the ``combination_database_data.py`` to generate all the combined data into one file for deep learning and further analysis


Note
~~~~~

- **For construction and evaluation of the deep learning model:**
  
  - To see how the deep learning pipeline is constructed, check the corresponding functions in the ``Code/model`` directory

- **For prediction of 343 yeast/fungi species via the deep learning model:**
  
  - To obtain prediction results for 343 yeast/fungi species based on the trained deep learning model, unzip the ``input.zip`` file in ``Data/input.zip`` and run the corresponding function in the ``Code/prediction`` directory


