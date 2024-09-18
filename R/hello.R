### load dependencies, datasets, and set up values for plotting ###

## dependencies
devtools::install_github("ZacharySMorris/GMM-Major-Axes", force = T)
library(GMMMajorAxes)
library(geomorph)

## load datasets
load("data/Combined_CrocDorsal.Rdata") # raw landmark data
load("data/skull_wireframe.Rdata") # landmark links for wireframe
load("data/classifier.Rdata") # specimen covariate data
load("data/Pseudosuchia.phy.Rdata") # phylogeny
# load("Chapter2Code/data/classifier.Rdata") # species occurance datas

Combined_CrocDorsal.pca <- gm.prcomp(Combined_CrocDorsal)

Combined_CrocDorsal.pca$scale

  
  
SubsettingGMM()