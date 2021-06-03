
Required software 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `MATLAB <http://www.mathworks.com/>`_ 9.1 (R2016b) or higher + Optimization Toolbox.
- The `COBRA toolbox for MATLAB <https://github.com/opencobra/cobratoolbox>`_.
- The `RAVEN toolbox for MATLAB <https://github.com/SysBioChalmers/RAVEN>`_.
- The `libSBML MATLAB API <https://sourceforge.net/projects/sbml/files/libsbml/MATLAB%20Interface>`_ (version 5.17.0 is recommended).
- For ploting violin plots, The toolbox `Violinplot-Matlab <https://github.com/bastibe/Violinplot-Matlab>`_
- For calculating Adjusted P value, function `PVAL_ADJUST <https://github.com/fakenmc/pval_adjust>`_
- In the Bayesian process, `IBM CPLEX 12.10 <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_ is used, `gurobi <https://www.gurobi.com>`_  can also be used after adapting the function ``abc_matlab_max.m`` with replacing ``ibm_cplex`` to ``gurobi`` when calling th function ``solveModel``.

- For other process, default cobra solver would be used, either `IBM CPLEX 12.10 <https://www.ibm.com/products/ilog-cplex-optimization-studio>`_ or `gurobi <https://www.gurobi.com>`_  is ok.

Usage
~~~~~
- Add the path of ``DLKcat`` to the MATLAB path. 
- If you just want to run part of the process, please download the raw results from the figshare link and combine that with the ``Results`` folder.


- **For extracting infos for kcat prediction using deep-learning model:**

  - The function is run on cluster, please adapt the unction ``initcluster.m`` for setting up the cluster path.
  
  - Adapt the function ``ForKcatPrediction/WriteFile_pre_cluster.m``, which calls the function ``ForKcatPrediction/writeFileForKcatPrediction.m`` for extracting information required for deep learning model prediction and generating a txt file for each model. Model files are from the `Yeast-Species-GEMs <https://github.com/SysBioChalmers/Yeast-Species-GEMs/tree/master/Reconstruction_script/ModelFiles/xml>`_ 

 
- **For reconstruction enzyme-constrained model:** Enzyme-constrained models can be used as any other metabolic model, with toolboxes such as COBRA or RAVEN. 

  - Preliminary step is to run the deep learning step for the prediction of the kcat. All predicted kcat files can be found in the figshare folder, please download the raw results from the figshare link and combine that with the ``Results`` folder
  
  - Adapt the function ``initcluster.m`` for setting up the cluster path
  
  - run the ``autoDLecModelGeneration.m`` to collect enzymedata which contains kcat and protein info and generate auto-ecGEM and DL-ecGEM.
  
  - run the ``BayesianModelGeneration_cluster.m`` to get the Bayesian-DL-ecGEM, species can be input by the index order in the 
  
  - run the ``getEmodel.m`` to get GECKO version ecGEMs

- **For running analysis of enzyme-constrained model and regenerate all figures:**
  
  - To regenerate all figures, run the corresponding figure functions in the ``Analysis``
  
  
Note
~~~~~

- If all you need is ecGEMs, please find the ecGEMs in the folder: ``Results/ecGEMs``. There are three version of ecGEMs available for each species. ``Bayesian-DL-ecGEMs`` version is recommended.

