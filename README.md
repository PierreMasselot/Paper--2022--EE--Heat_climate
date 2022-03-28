# Fully reproducible R code and data from Masselot et al. (2022) Environmental Epidemiology

Data and R code to reproduce the analysis and results of the article :

------

It implements the two-stage methodology used to predict annual heat-related mortality using low-frequency climate indices. Specifically, it runs the time-varying DLNM and extracts annual attributable fractions (AF) in the first stage, and the functional model regressing these AFs on monthly climate index values in the second stage. The code then produces output Tables and Figures. At each stage, it loops through the two cities analysed in the paper: Montréal and Québec.

## Details

It is composed of 6 scripts.
0. *0_Packages_Parameters.R*: Load necessary packages and describes the different parameters used in the analysis.
1. *1_PrepData.R*: Load and prepare data, including daily mortality and temperature time series, as well as monthly climate values.
2. *2_AttrPred_timeVaryingDLNM.R*: Runs the first-stage on each city, i.e. the time-varying DLNM, computation of AFs, and empirical CIs.
3. *3_functionalRegression.R*: Runs the second stage, namely the functional regression of annual AF on monthly climate indices, and bootstrap CIs. **Computationally intensive**, might take few hours.
4. *4_Results.R*: Extract results and produce tables.
5. *5_Plots.R*: Produces plots.

Note that additional shapefiles are needed to produce Figure 1. These can be downloaded here: https://www.donneesquebec.ca/recherche/dataset/decoupages-administratifs#
