# CACTI_Proximity_Test
Source code for CACTI and Proximity tests for use in spatial transcriptomics 
Code for individual function used for CACTI and the spatial proximity test can be found in CACTI_PROXIMITY_FUNCTIONS.R
lambda_calc computes the maximum lambda value such that the classification error with respect to some low resolution classification is below a threshold alpha
res_calc determines a resolution of louvain clustering such that the resulting number of clusters produced is within some pre-specified range
prune_clust tests to see if the clusters found by CACTI are significantly dissimilar from one another 
proximity_test performs the similarity test outlined in the materials and methods  
