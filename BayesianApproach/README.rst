
Required software 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MATLAB <http://www.mathworks.com/>`_ 9.1 (R2016b) or higher + Optimization Toolbox.
- The `COBRA toolbox for MATLAB <https://github.com/opencobra/cobratoolbox>`_.
- The `RAVEN toolbox for MATLAB <https://github.com/SysBioChalmers/RAVEN>`_.
- The `libSBML MATLAB API <https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface>`_ (version 5.17.0 is recommended).
- For ploting violin plots, The toolbox `Violinplot-Matlab <https://github.com/bastibe/Violinplot-Matlab>`_
- For calculating Adjusted P value, function `PVAL_ADJUST <https://github.com/fakenmc/pval_adjust>`_
- In the Bayesian process, `IBM CPLEX 12.10 <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_ is used, `gurobi <https://www.gurobi.com>`_  can also be used after adapting the function ``abc_matlab_max.m`` with replacing ``ibm_cplex`` to ``gurobi`` when calling th function ``solveModel``.

- For other process, default solver in the cobra toolbox would be used, either `IBM CPLEX 12.10 <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_ or `gurobi <https://www.gurobi.com>`_  is ok.

Usage
~~~~~
- Add the path of ``DLKcat`` to the MATLAB path. 
- If you just want to run part of the process, please download the raw results from the `zenode link <https://doi.org/10.5281/zenodo.5164210>`_ and combine that with the ``Results`` folder.


- **For extracting infos for kcat prediction using deep-learning model:**

  - The function is run on cluster, please adapt the unction ``initcluster.m`` for setting up the cluster path.
  
  - Adapt the function ``ForKcatPrediction/WriteFile_pre_cluster.m``, which calls the function ``ForKcatPrediction/writeFileForKcatPrediction.m`` for extracting information required for deep learning model prediction and generating a txt file for each model. Model files are from the `Yeast-Species-GEMs <https://github.com/SysBioChalmers/Yeast-Species-GEMs/tree/master/Reconstruction_script/ModelFiles/xml>`_ 

 
- **For reconstruction enzyme-constrained model:** Enzyme-constrained models can be used as any metabolic model, with toolboxes such as COBRA or RAVEN. 

  
  - Preliminary step is to run the deep learning step for the prediction of the kcat. Check the `deeplearningApproach folder  <https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach>`_ for how th prediction was done. All predicted kcat files can be found in the `zenode link <https://doi.org/10.5281/zenodo.5164210>`_ folder, please download the raw results from the `zenode link <https://doi.org/10.5281/zenodo.5164210>`_  and combine that with the ``Results`` folder. 
 
  - This step is run on cluster, please adapt the unction ``initcluster.m`` for setting up the cluster path. create your own bash file for runing jobs on th cluster.
  - run the ``classicDLModelGeneration_cluster.m`` to collect enzymedata which contains kcat and protein info and generate Classical-ecGEM and DL-ecGEM. Species should be input as the index order in the `Strain.txt <https://github.com/SysBioChalmers/DLKcat/blob/master/BayesianApproach/Code/ecGEMconstruction/Strain.txt>`_.
  
  - run the ``BayesianModelGeneration_cluster.m`` to get the Posterior_mean-ecGEM. Species should be input as the index order in the `Strain.txt <https://github.com/SysBioChalmers/DLKcat/blob/master/BayesianApproach/Code/ecGEMconstruction/Strain.txt>`_.
  
  - run the ``getEmodel.m`` to get GECKO version ecGEMs from the enzymedata and GEM. After that, all GECKO functions can be used for further analysis. all ecGEMs are stored in the `zenode link <https://doi.org/10.5281/zenodo.5164210>`_.

- **For running analysis of enzyme-constrained model and reproduce all figures:**
  
  - To reproduce all figures, run the corresponding figure functions in the ``code/Analysis``, raw results are required from: `zenode link <https://doi.org/10.5281/zenodo.5164210>`_, please download the raw results from the `zenode link <https://doi.org/10.5281/zenodo.5164210>`_ and combine that with the ``Results`` folder.
  - All anlysis can be done in local computer.
  
  
Note
~~~~~

- If all you need is ecGEMs for 343 yeast/fungi species, please find the ecGEMs in the folder: ``Results/ecGEMs`` from th `zenode link <https://doi.org/10.5281/zenodo.5164210>`_. There are three version of ecGEMs available for each species. ``Posterior_mean-ecGEMs`` version is recommended.

