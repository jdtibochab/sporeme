# SPOREme

SPOREme contains all the information and building scripts required to construct a ME-model for *B. subtilis* sporulation (SporeME2) based on [BACILLUSme](https://github.com/jdtibochab/bacillusme) ([*i*JT964-ME](https://www.nature.com/articles/s41540-022-00259-0)), which was based on [COBRAme](https://github.com/sbrg/cobrame). See the COBRAme
[github](https://github.com/sbrg/cobrame) and 
[documentation](https://cobrame.readthedocs.io) for further instructions on 
how to build and solve ME-models. Minor changes were necessary for COBRAme to function with this organism, so find the modified COBRAme we used for *i*JT964-ME [here](https://github.com/jdtibochab/cobrame).

SporeME2 files and model are available in bacillusme/analysis/spore. The notebooks are listed in numerical order and reproduce the reconstruction and findings of the manuscript.

If you are using BACILLUSme or *i*JT964-ME, please cite [doi:10.1038/s41540-022-00259-0](https://www.nature.com/articles/s41540-022-00259-0).

## Installation

1. clone the repository
2. run ```python setup.py develop --user```
3. Run build_me_model.ipynb to construct and save *i*JT964-ME.
4. The saved model can be solved using solve_demo.ipynb