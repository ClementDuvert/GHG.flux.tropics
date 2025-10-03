# GHG.flux.tropics

The GHG.flux.tropics package consists of MATLAB scripts designed to calculate greenhouse gas (GHG) fluxes from river and lake datasets in tropical and subtropical regions. The package enables the following:
- Uploading and preprocessing of GHG data for rivers and/or lakes.
- Standardisation of units.
- Generation of mean or median values per site and system.
- Estimation of median fluxes per climate region and system size (e.g. Strahler order) using a bootstrapping approach.
- Scaling of these fluxes to the broader (sub)tropics.

The package includes the necessary raw data files and is related to the following article: Duvert C. et al. (2025) Hydroclimate and landscape diversity drive highly variable greenhouse gas emissions from tropical and subtropical inland waters, Nature Water 3, https://doi.org/10.1038/s44221-025-00522-8.

### System requirements
- MATLAB version: R2022a Update 1 or later.
- Required toolboxes: Statistics and Machine Learning Toolbox, Optimization Toolbox.
- Hardware: No specific hardware requirements beyond standard MATLAB-compatible devices.

### Instructions
1. Set up directories. Ensure that all input datasets and scripts are either in the working directory or add their location as follows:
```MATLAB
addpath('/Users/...') % Add the directory where functions and/or datasets are stored
```
2. Run the analysis. Start by running the call_functions.m script, which manages the execution of all subsequent functions in the package. To begin, type the following command in MATLAB's command window:
```MATLAB
call_functions
```
A pop-up will prompt you to select whether to analyse river or lake data.

### Expected outputs
Each function generates specific outputs, as detailed in the header of each script. Key outputs are displayed in MATLAB’s command window as the functions run.
The expected run time for the river data is between 45 seconds and 2 minutes, while the expected run time for the lake data is between 20 seconds and 1 minute.
