# SPOREme

SPOREme contains all the information and building scripts required to construct a ME-model for *B. subtilis* sporulation (SporeME2) based on [BACILLUSme](https://github.com/jdtibochab/bacillusme) ([*i*JT964-ME](https://www.nature.com/articles/s41540-022-00259-0)), which was based on [COBRAme](https://github.com/sbrg/cobrame). See the COBRAme
[github](https://github.com/sbrg/cobrame) and 
[documentation](https://cobrame.readthedocs.io) for further instructions on 
how to build and solve ME-models. Minor changes were necessary for COBRAme to function with this organism, so find the modified COBRAme we used for *i*JT964-ME [here](https://github.com/jdtibochab/cobrame).

SporeME2 files and model are available in bacillusme/analysis/spore. The notebooks are listed in numerical order and reproduce the reconstruction and findings of the manuscript.

If you are using BACILLUSme or *i*JT964-ME, please cite [doi:10.1038/s41540-022-00259-0](https://www.nature.com/articles/s41540-022-00259-0).


Install locally
====================
1. clone the repository
2. run ```python setup.py develop --user```

Note: You will need to manually download and install some necessary requirements: [COBRApy](https://github.com/opencobra/cobrapy/releases/tag/0.5.11) 0.5.11, [COBRAme (modified)](https://github.com/sbrg/cobrame), and [solveME](https://github.com/sbrg/solvemepy). solveME requires qMINOS, a quad-precision solver that is available upon request to Prof. Michael A. Saunders at Stanford University

Install using docker
====================
1. Clone repository and navigate to sporeme/
2. ``docker build --file "./Dockerfile-Python3.7" . -t "python3.7-cobrame"``
3. ``docker run --detach -p 10000:8888 -v USER/PATH/TO/coralme/:/opt/notebooks/ python3.10-coralme``
4. In your browser, go to ``localhost:10000``

Understanding the layout of this repository
===========================================
There are two main folder locations that correspond to the two main steps of this project: (1) Reconstruction of ME-models, including updating *i*JT964-ME, and reconstructing the forespore and mother cell ME-models, and (2) Reconstruction and analysis of the ME2-model of _B. subtilis_ sporulation, SporeME2. The layout is as follows:

1. Reconstruction of ME-models: **Location**: [sporeme/bacillusme/](https://github.com/jdtibochab/sporeme/tree/main/bacillusme/)

   - Notebooks 1, 2, and 3 reconstruct the updated ME-model of _B. subtilis_, and individual the ME-models of the mother cell and forespore, respectively.
   - The data used for ME-model building (following previously used formatting and structure, see ECOLIme and BACILLUSme), is in [building_data](https://github.com/jdtibochab/sporeme/tree/main/bacillusme/building_data). Notebooks 1 and 2 generate input files used for ME-model building (already provided in building_data/, so there is no need to re-run).
   - The output models are stored in [me_models](https://github.com/jdtibochab/sporeme/tree/main/bacillusme/me_models), as ``solution.pickle``, ``solution_mother.pickle``, and ``olution_spore.pickle``.

2. Reconstruction of SporeME2: **Location**: [sporeme/bacillusme/analysis/spore/](https://github.com/jdtibochab/sporeme/tree/main/bacillusme/analysis/spore/)

   - Notebooks are numbered in logical order. The first steps (Notebooks 1.) reconstruct SporeME2, while the subsequent ones load it to analyze metabolic mechanisms (Notebooks 2.) as well as proteome-wide protein essentiality (Notebooks 3.).
   - **Notebook 1.1**: Estimates biomass composition corrections to the _B. subtilis_ biomass reaction from various reports.
   - **Notebook 1.2**: Merges the ME-models of the mother cell and the forespore, implements special constraining from literature, and generates SporeME2.
   - **Notebook 1.3**: Visualizes the properties of SporeME2.
   - **Notebook 2.1.1**: Estimates ATP supply mechanisms to the forespore.
   - **Notebook 2.1.2**: Estimates ATP supply mechanisms to the forespore if direct ATP transport is allowed.
   - **Notebook 2.2**: Analyzes the amino acid supply mechanisms to the forespore.
   - **Notebook 2.3**: Analyzes NTP supply mechanisms to the forespore.
   - **Notebook 2.4**: Visualizes the activity and cell-specific expression of the QA channel.
   - **Notebook 3.1.1.**: Runs the calculation of the cell-specific protein essentiality analysis.
   - **Notebook 3.1.2.1.**: Reads and analyzes the results of the protein essentiality analysis.
   - **Notebook 3.1.2.2.**: Prepares GO enrichment analysis ([The Gene Ontology Resource](https://geneontology.org/)) input files and analyzes the output.
   - **Notebook 3.1.3.**: Runs conditional activity/inactivity analysis to identify proteome-wide effects of protein depletions.
   - **Notebook 3.1.4.**: Visualizes results from the conditional activity/inactivity analysis and prepares inputs for visualization in Cytoscape.
   - **Notebook 3.2.1.**: Runs the calculation of the cell-specific protein essentiality analysis in the **naive** SporeME2 (containing no cell-specific depletion constraining). 
   - **Notebook 3.2.2.**: Analyzes the results of the protein essentiality analysis on the **naive** model.
   - **Notebook 3.2.3.**: Runs conditional activity/inactivity analysis on the **naive** model.
   - **Notebook 3.2.4.**: Visualizes results from the conditional activity/inactivity analysis in the **naive** model and prepares inputs for visualization in Cytoscape.
