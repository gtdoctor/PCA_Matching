# PCA Ancestry matching algorithm 

GOAL: to find up to "nc" unique controls within a certain "threshold" "distance" of each case, based on principal components of a genetic relationship matrix. 

## INPUTS: 
Requires:
plink.eigenvec 
plink.eigenval
plink.fam or plink.psam

## OUTPUTS:
Part A: 
distance matrix .Rdata

Part B: 
IIDs of kept and excluded cases and controls
Scree plot
PC plots of PC1 vs PC2-5 with all cases vs controls; a plot splitting out cases and controls kept and removed; and a plot showing just kept cases and controls  

## Requirements:
library(ggplot2)
library(grid)
library(gridExtra)
library(parallel)

## PART A: Distance Matrix
Creates a matrix of distances between each case and all controls.  

### Parameters: 
Number of principal components to include (max 10)
Euclidean or manhattan distances
  - Manhattan = sum of absolute distances from each principal component axis included)
  - Euclidean = Euclidean distance through n PC dimenensions
    
Weighted (by correspoding eigenvector) or unweighted distances:
  - the distance contributed by each PC can be weighted by the eigenvalue of that PC
    
N_chunks for parallel processing of large datasets.
Small sample (enable for testing)


## Part B: Matching
Based on a distance matrix between cases and controls (choosing one of the above methods, finds up to "nc" unique controls within a certain "threshold" distance of each case.
Cases are removed if they match to fewer than "ncmin" controls (along with the controls that were being assigned to a removed case).
When no further case matches are possible, or all cases have  "nc" matches, all unmatched controls are removed.

### Further details of matching: 
Ignoring all distances greater than the threshold, for each case the nearest-neighour control is identified
The case with the furthest near-neighbour is assigned that control. 
The control is then removed; then that case and any case assigned to the same control, have their next nearest neighbour identified.
The case with the furthest near-neighbour is assigned that control. And so on. 

### Parameters: 
"Threshold" is defined as some proportion ("threshmod") of the median value in the  distance matrix.
"nc" - number of controls sought per case (i.e max number of controls to be matched.
"ncmin" - minimum acceptable number of controls - cases excluded if this cannot be met. 


## COMMENTS
There are no default settings as such and results are sensitive to the method chosen, especially weighting and the number of controls sought.   
I recommend using Euclidean weighted distances with 10 PCs (based on some empirical testing of data from 1000 Genomes Project and UKB). 
In UKB, seeking 7 cases with threshold distance of 0.2 of median distance seems optimal. 
Using 40 cores on a linux server, on UK Biobank scale data, part A takes approx 10 minutes. Part B quicker. 

## Developers
Ancestry matching algorithm first developedby AP Levine. https://github.com/APLevine
Updated 2024 by GT Doctor UCL Centre for Genetics and Genomics ttps://github.com/gtdoctor
