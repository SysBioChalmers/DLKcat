DLKcat
======

<p align="center">
  <img  src="doc/logo.png" width = "400">
</p>


Introduction
------------
[中文版](https://github.com/SysBioChalmers/DLKcat/blob/master/README_ZH.md)

The **DLKcat** toolbox is a Matlab/Python package for prediction of
kcats and generation of the ecGEMs. The repo is divided into two parts:
`DeeplearningApproach` and `BayesianApproach`. `DeeplearningApproach`
supplies a deep-learning based prediction tool for kcat prediction,
while `BayesianApproach` supplies an automatic Bayesian based pipeline
to construct ecModels using the predicted kcats.

Usage
-----

-   Please check the instruction `README` file under these two section
    `Bayesianapproach` and `DeeplearningApproach` for reporducing all figures in
    the paper.
-   For people who are interested in using the trained deep-learning
    model for their own kcat prediction, we supplied an example. please
    check usage for **detailed information** in the file
    [DeeplearningApproach/README](https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach)
    under the `DeeplearningApproach`.

    > -   `input` for the prediction is the `Protein sequence` and
    >     `Substrate SMILES structure/Substrate name`, please check the
    >     file in
    >     [DeeplearningApproach/Code/example/input.tsv](https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach/Code/example)
    > -   `output` is the correponding `kcat` value

Citation
-----

- Please cite the paper [Deep learning-based kcat prediction enables improved enzyme-constrained model reconstruction](https://www.nature.com/articles/s41929-022-00798-z)""


Notes
-------
We noticed there is a mismatch of reference list in Supplementary Table 2 of the publication, therefore we made an update for that. New supplementary Tables can be found [here](https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach/Results/figures)

Contact
-------

-   Feiran Li ([@feiranl](https://github.com/feiranl)), Chalmers
    University of Technology, Gothenburg, Sweden
-   Le Yuan ([@le-yuan](https://github.com/le-yuan)), Chalmers
    University of Technology, Gothenburg, Sweden

Last update: 2023-07-04
