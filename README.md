# Neuromodulation-Induced Normalization of Cortical Metastable Dynamics in Parkinson's Disease
This repository contains code in support of my project, "Neuromodulation-Induced Normalization of Cortical Metastable Dynamics in Parkinson's Disease".

## Code
The [MS-DBS-TI]([MS-DBS-TI/) folder contains the scripts required to conduct the main analyses:
- [StateSpace_modelling.m](MS-DBS-TI/StateSpace_modelling.m) is the Weighted Eigenvector Dynamics Analysis (WEiDA) pipeline identifying 4 MS states based on the DBS cohort. Input: preprocessed Bold time series. Output: cluster-centroid WEiDA eigenvectors for MS states; state fractional occupancy, dwell time, transition probabilities and state time series for each subject. 
- [WEiDA_pred.m](MS-DBS-TI/WEiDA_pred.m) infers state switching of the tTIS cohort based on the WEiDA state space model built in DBS. nput: preprocessed Bold time series. Output: state fractional occupancy, dwell time, transition probabilities and state time series for each subject.
- [overlapRSN.m](MS-DBS-TI/overlapRSN.m) compares the spatial profiles (cluster-centroid WEiDA eigenvectors) of the 4 MS states and the seven canonical functional networks (RSNs) delineated by (Yeo et al., 2011).
- [MSanalysis.m](MS-DBS-TI/MSanalysis.m) calculates state metastability for each subject.
- [calcu_FS.py](MS-DBS-TI/calcu_FS.py) calculates functional segregation and nodal functional segregation for each subject.
- [AHBA_gene_expression.py](MS-DBS-TI/AHBA_gene_expression.py) investigates the spatial correlation between the neuromodulation-induced functional segregation alteration and the PD-related gene expression (data from AHBA) using a spatial autocorrelation-preserving null model in neuromaps (https://netneurolab.github.io/neuromaps/).
- [across_atlas.py](MS-DBS-TI/across_atlas.py) makes inferences about the topographic correlations between the WEiDA-derived brain maps of the four MS states based on AAL and Schaefer atlases through a spatial autocorrelation-preserving null model in neuromaps toolbox.
- [dynamic_features.py](MS-DBS-TI/dynamic_features.py) organizes neuromodulation-induced symptom improvements, all dynamic features changes (including state fractional occupancy, dwell time, transition probabilities, metastability and functional segregation), and their relationships.
- [group_comparison.py](MS-DBS-TI/group_comparison.py) makes group comparisons of dynamic features.
- [plot_fig.ipynb](MS-DBS-TI/plot_fig.ipynb) plots figures in this project.




