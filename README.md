# iHil3966
We reconstructed the genome-scale metabolic network of Hermetia illucens using the [REVAN tool](https://github.com/SysBioChalmers/RAVEN/tree/main) and conducted several metabolic analyses on it.

## Prerequisites
The main programs were developed and tested in MATLAB R2018a [COnstraint-Based Reconstruction and Analysis (The COBRA Toolbox)](https://opencobra.github.io/cobratoolbox/stable/) and [COBRApy](https://cobrapy.readthedocs.io/en/latest/) are required to perform the analysis. Check [here](https://opencobra.github.io/cobratoolbox/stable/installation.html) for the installation guidance of COBRA Toolbox. The programs were tested for COBRA Toolbox - 2024 version, but should be compatible with an earlier version.
The Linear Program (LP) and Mixed-Integer Linear Problem (MILP) solver used in the study was [gurobi 9.0.3](https://www.gurobi.com/). The built-in solver interface of COBRA Toolbox was used.

## Identification the essential exogenous nutrients
Identified the essential exogenous nutrients required by BSF using this model (`exogenous_nutrients.ipynb`). 

## Context metabolic networks
Utilizing the [iMAT++ algorithm](https://www.embopress.org/doi/full/10.15252/msb.20209649), we generated specific context metabolic networks for different life stages of BSF (`imatplus_iHil3966.m`). Before running imatplus_iHil3966.m, please install iMAT++ from [here](https://github.com/WalhoutLab/MERGE).

## Identification of nutrient interventions
Using [a nutritional algorithm](https://www.sciencedirect.com/science/article/pii/S1096717623000186), we optimized the nutritional additives for enhancing BSF's fatty acid production, and calculated the effects of individual and various combinations of nutritional additives on BSF biomass and fatty acid production (`nutrition_iHil3966.m`).  Before running nutrition_iHil3966.m, please install nutritionAlgorithm from [here](https://github.com/opencobra/cobratoolbox/tree/master/src/analysis/nutritionAlgorithm).
