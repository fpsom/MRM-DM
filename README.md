# MRM-DM
This is a companion repository for the Data Mining course [1] of the MSc in Medical Research Methodology provided by the Medical School of the Aristotle University of Thessaloniki.

The repository is organized as follows:

## Data
The *data* folder contains both the data file used in the examples (Data_Cortex_Nuclear in xls and csv format) and the pdf files of two papers that have been published on the analysis of the data. The original data file was retrieved from the UC Irvine Machine Learning Repository [2].


## Scripts

1. MouseClusteringMethods.R
    This script comprises three steps
    * Data restructuring: the data is cleaned and re-annotated for the clustering process
    * k-means clustering: setting a value for k, it produces a visual representation of the clustering using selected attributes/columns
    * hiearchical clustering: selecting a method for the process, it creates an image of the produced tree

2. MouseClusteringHeatmap.R
    This script aims to provide a real-world approach to how clustering is used in research. It produces a heatmap of the protein expression, through the following steps:
    * Data restructuring: the data is cleaned and re-annotated for the clustering process
    * Heatmap initializaion: the expression ranges should be identified for better representation
    * Definition of the clustering methods: values can be clustered on both axis (samples and proteins).
    * Create the final heatmap image


## Student contributions
The students are expected to investigate and improve the clustering process by changing the parameters and functions in both scripts.

The specific fields that are should be edited are highlighted with the following comment-type:

> +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
>  *Specific instructions here*
> ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## References
[1] http://mrm.med.auth.gr/data-mining/
[2] https://archive.ics.uci.edu/ml/datasets/Mice+Protein+Expression