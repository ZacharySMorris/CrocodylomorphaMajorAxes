# Crocodylomorpha Major Axes
This project contains the data and code used to perform the exploratory and statistical analyses of Morris and colleagues (YEAR). Using the original functions to calculate and compare the major axes of cranial variation across pan-crocodylian evolution.  The scripts within this package will walk you through installing necessary packages, loading datasets, and performing analyses.

#Install GMM-Major-Axes Package
The functions used to perform Major Axes of Subgroups Analysis have been packaged for easy installation and use.
```
devtools::install_github("ZacharySMorris/GMM-Major-Axes", force = T)
library(GMMMajorAxes)
```
#Analytical Outline
1. Load data
2. Perform GPA and PCA
3. Subset shape data and specimen covariate data
4. Time calibrate phylogeny
5. 
