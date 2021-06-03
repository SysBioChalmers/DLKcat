
DLKcat
============

Introduction
------------

The **DLkcat** toolbox is a Matlab/Python package for prediction of kcats and generation of the ecGEMs. The repo is divided into two parts: ``DeeplearningApproach`` and ``BayesianApproach``. ``DeeplearningApproach`` supplies a deep-learning based prediction tool for kcat prediction, while ``BayesianApproach`` supplies an automatic Bayesian based pipeline to construct ecModels using the predicted kcats.


Usage
------------
* Please check the instruction ``README`` file under these two section ``Bayesianapproach`` and ``Deeplearning`` for reporducing all figures in the paper.

* For people who are just interested in using the trained deep-learning model to predict the kcat, we supplied an example for users to direcly use the tool for their own kcat prediction. please also refer to th `DeeplearningApproach/README <https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach>`_  file under the ``Deeplearning``.
- ``input`` for the prediction is the  ``Protein sequence`` and ``Substrate SMILES structure/Substrate name``, please check the file in `DeeplearningApproach/Code/example/input.tsv <https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach/Code/example>`_ 
- ``output`` is the correponding ``kcat`` value 



Contact
-------------------------------

* Feiran Li (`@feiranl <https://github.com/feiranl>`_), Chalmers University of Technology, Gothenburg, Sweden
* Le Yuan (`@le-yuan <https://github.com/le-yuan>`_), Chalmers University of Technology, Gothenburg, Sweden


Last update: 2021-06-02
