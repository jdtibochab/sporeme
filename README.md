# BACILLUSme

BACILLUSme contains all the information and building scripts required to construct a ME-model for *B. subtilis*, *i*JT964-ME using
[COBRAme](https://github.com/sbrg/cobrame). See the COBRAme
[github](https://github.com/sbrg/cobrame) and 
[documentation](https://cobrame.readthedocs.io) for further instructions on 
how to build and solve ME-models. This package is based on and is homologous to [ECOLIme](https://github.com/sbrg/ecolime). Minor changes were necessary for COBRAme to function with this organism, so find the modified COBRAme we used for *i*JT964-ME [here](https://github.com/jdtibochab/cobrame).

If you are using BACILLUSme or *i*JT964-ME, please cite [doi:10.1038/s41540-022-00259-0](https://www.nature.com/articles/s41540-022-00259-0).


## Installation

1. clone the repository
2. run ```python setup.py develop --user```
3. Run build_me_model.ipynb to construct and save *i*JT964-ME.
4. The saved model can be solved using solve_demo.ipynb
