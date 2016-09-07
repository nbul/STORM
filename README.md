# STORM
This script analyse ring structures visualised in STORM

The Input data should be in a single folder in .csv format 
The head script is "STORMcsv.m", which automatically will run all 
other scripts.

Then, select a folder with .csv files to analyse from a
dialog box, they can be at any location as long as all .csv files are
in a single folder and named by sequential numbers (1,2,3...).

Change the cut-off values for minimum goodness-of-fit of radial profiles
and the column which is used to get intensity of individual events.

Currently, the script uses three approaches to get 9-fold symmetry. One is
pair-wise distance between clusters following kmeans clustering of the data. 
Another is "circular" profile calculations: angle from centre is measured for all
events, they are binned by angle and all intensities in each bin are
summed up, and align to start with max peak. Then either all profiles
are summed up to give a summarised profile, or all peaks in each
profile are identified to give a distribution of peaks. It does not
produce very conclusive results, when using only 10 rings that I have,
but it would be interesting to run it on a larger dataset.

The output data:

aligned_distribution.tif - summarised circular distribution of signal
from center normalised to signal intensity and shifted to start with
highest peak

angles_clusters.csv - distributions of pair-wise angles between clusters
identified with k-means after and before removing 2 weakest clusters

distribution_circular.csv - individual circular distributions from center
normalised to signal intensity and shifted to start with
highest peak (first column is the bin, last column - summarised values)

distribution_radial - individual radial distributions from center
normalised to signal intensity (first column is the bin, last column -
summarised values) 

peaks_distribution.csv - destribution of all peaks found in individual
circular destributions.

peaks_distribution.tif - visualisation of peaks_distribution.csv.

radiuses_curve.csv - data about individual distributions fitted with the
equation from the paper (last raw is the summarised ring)

radiuses.csv - data about individual distributions fitted with gaussian
(last raw is the summarised ring)

summarized_distributions_modified.tif - visualisation of
radiuses_curve.csv 

summarized_distributions.tif - visualisation of radiuses.csv

summarized_ring_fitted.tif - visualisation of summarised ring fitted with
the equation from the paper

summarized_ring.tif - visualisation of non-fitted summarised ring
