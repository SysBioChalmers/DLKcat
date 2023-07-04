DLKcat
======

<p align="center">
  <img  src="doc/logo.png" width = "400">
</p>


关于DLKcat
------------
 [English version](https://github.com/SysBioChalmers/DLKcat/blob/master/README.md)
 
**DLKcat**工具箱是一个基于Matlab/Python开发的程序，用于预测 kcats 和生成 ecGEMs。工具箱分为两部分：`DeeplearningApproach` 和 `BayesianApproach`。`DeeplearningApproach` 是一种基于深度学习技术来预测 kcat 的预测工具，而 `BayesianApproach`是一个通过贝叶斯方法构建 ec模型组 并预测 kcats 的自动程序。

用法
-----
-   请参考 `Bayesianapproach` 和 `DeeplearningApproach` 下的 `README` 进行使用.
-   我们为有兴趣使用 经过训练的深度学习模型 来预测 kcats 的人提供了一个**例子**，在`DeeplearningApproach`下的[DeeplearningApproach/README](https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach)。
    > -   `输入`预测的 `Protein sequence` 和
    >     `Substrate SMILES structure/Substrate name`, 请参考
    >     [DeeplearningApproach/Code/example/input.tsv](https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach/Code/example)
    > -   `输出`预测的`kcat` 

引用
-----

- 请引用论文 [Deep learning-based kcat prediction enables improved enzyme-constrained model reconstruction](https://www.nature.com/articles/s41929-022-00798-z)


注意事项
-------
我们注意到该文件补充表2中的参考文献列表不匹配，因此我们对此进行了更新。可以在[此处](https://github.com/SysBioChalmers/DLKcat/tree/master/DeeplearningApproach/Results/figures)找到新的补充表格。

联系方式
-------

-   Feiran Li ([@feiranl](https://github.com/feiranl)), Chalmers
    University of Technology, Gothenburg, Sweden
-   Le Yuan ([@le-yuan](https://github.com/le-yuan)), Chalmers
    University of Technology, Gothenburg, Sweden

名词解释
-------
- `kcats`：
最有优条件下单位酶催化生成底物的速率，定义上等于最大反应速率除以酶浓度
-  `ecGEM`：基因组规模代谢模型
    
更新时间: 2023-07-04
