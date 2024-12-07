# Crocodylomorpha Major Axes
This project contains the data and code used to perform the exploratory and statistical analyses of Morris and colleagues (YEAR). Using the original functions to calculate and compare the major axes of cranial variation across pan-crocodylian evolution.  The scripts within this package will walk you through installing necessary packages, loading datasets, and performing analyses.

## Install GMM-Major-Axes Package
The functions used to perform Major Axes of Subgroups Analysis have been packaged for easy installation.
```
devtools::install_github("ZacharySMorris/GMM-Major-Axes", force = T)
library(GMMMajorAxes)
```
## Analytical Outline
1. Load raw landmark data, specimen covariate data, species level phylogeny
2. Perform GPA and PCA
3. Subset shape data and specimen covariate data
4. Estimate morphological disparity & test for differences
5. Time calibrate phylogeny & calculate phylogenetic signal
6. Estimate & Compare MAs of Extant Crocodylian Ontogeny
7. Estimate MA of Extant Crocodylian Adult Ecology
8. Compare MAs of Ecomorph Ontogeny to Extant Ecology
9. Estimate MAs of Extinct Pseudosuchian Subgroups
10. Compare MAs of Extinct Pseudosuchians with Extant Ecology & Ontogeny
11. Estimate & Compare MAs of Geologic Time Periods
12. Prepare dataset of mean adult shape for all species and mean embryo shapes
13. Perform embryonic 'convergence' test
14. Perform extant ontogeny trajectory analysis

|Script Name|Steps|
|---|---|
|Loading & Subsetting|1,2,3|

