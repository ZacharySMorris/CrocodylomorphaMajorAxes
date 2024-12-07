### load dependencies, datasets, and set up values for plotting ###

## dependencies
devtools::install_github("ZacharySMorris/GMM-Major-Axes", auth_token = "abc", force = T)
library(GMMMajorAxes)

detach("package:GMMMajorAxes", unload=TRUE)

library(geomorph)
library(convevol)

detach("package:convevol", unload=TRUE)


## load datasets
load("data/Combined_CrocDorsal.Rdata") # raw landmark data
load("data/skull_wireframe.Rdata") # landmark links for wireframe
#load("data/classifier.Rdata") # specimen covariate data
load("data/revised_classifier.Rdata") # revised specimen covariate data
#load("data/Pseudosuchia.phy.Rdata") # phylogeny
load("data/revised_Pseudosuchia.phy.Rdata") # revised phylogeny
#load("data/Species_Durations.Rdata") # species occurrence data
#load("data/Species_Duration_Table.Rdata") # species occurrence data
load("data/revised_Species_Duration_Table.Rdata") # revised species occurrence data

#revised_classifier <- read.csv("/Users/zach/Dropbox/NHM Summer 2024/MajorAxes_NatureEE_Revisions/MajorAxes_Covariates_Revised.csv")
#classifier <- read.csv("/Users/zmorris/Library/CloudStorage/Dropbox/NHM Summer 2024/MajorAxes_NatureEE_Revisions/MajorAxes_Covariates_Revised.csv")

## perform gpa & pca
Combined_CrocDorsal.gpa <- gpagen(Combined_CrocDorsal)

#rownames(Species_Duration_Table) <- gsub("Metriorhynchus_superciliosum","Thalattosuchus_superciliosus",rownames(Species_Duration_Table))
#rownames(Species_Duration_Table) <- gsub("Steneosaurus_leedsi","Charitomenosuchus_leedsi",rownames(Species_Duration_Table))
#rownames(Species_Duration_Table) <- gsub("Hyposaurus_natator","Hyposaurus_rogersii",rownames(Species_Duration_Table))

Combined_CrocDorsal.pca <- gm.prcomp(Combined_CrocDorsal.gpa$coords)
##

## make separate extant, 'fossil', adult, and embryo versions of the data
# extant subset
Extant_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
Extant_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Extant",classifier$Extant)]
Extant_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Extant",classifier$Extant)]

Extant_CrocDorsal.pca <- Combined_CrocDorsal.pca
Extant_CrocDorsal.pca$x <- Combined_CrocDorsal.pca$x[grep("Extant",classifier$Extant),]

Extant_classifier <- classifier[grep("Extant",classifier$Extant),]
Extant_classifier$Shape <- factor(Extant_classifier$Shape)
Extant_classifier$Shape_G.T <- factor(Extant_classifier$Shape_G.T)
Extant_classifier$Age <- factor(Extant_classifier$Age)
#
# fossil subset
Fossil_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
Fossil_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Fossil",classifier$Extant)]
Fossil_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Fossil",classifier$Extant)]

Fossil_CrocDorsal.pca <- Combined_CrocDorsal.pca
Fossil_CrocDorsal.pca$x <- Combined_CrocDorsal.pca$x[grep("Fossil",classifier$Extant),]

Fossil_classifier <- classifier[grep("Fossil",classifier$Extant),]
#
# adult & subadult subset
Adult_CrocDorsal.gpa <- Combined_CrocDorsal.gpa
Adult_CrocDorsal.gpa$coords <- Combined_CrocDorsal.gpa$coords[,,grep("Adult|Subadult",classifier$Age)]
Adult_CrocDorsal.gpa$Csize <- Combined_CrocDorsal.gpa$Csize[grep("Adult|Subadult",classifier$Age)]

Adult_CrocDorsal.pca <- Combined_CrocDorsal.pca
Adult_CrocDorsal.pca$x <- Combined_CrocDorsal.pca$x[grep("Adult|Subadult",classifier$Age),]

Adult_CrocDorsal.classifier <- classifier[grep("Adult|Subadult",classifier$Age),]
#

# Embryo subset
Embryo_CrocDorsal.gpa <- Extant_CrocDorsal.gpa
Embryo_CrocDorsal.gpa$coords <- Extant_CrocDorsal.gpa$coords[,,grep("Late-Stage Embryo|Mid-Stage Embryo",Extant_classifier$Age)]
Embryo_CrocDorsal.gpa$Csize <- Extant_CrocDorsal.gpa$Csize[grep("Late-Stage Embryo|Mid-Stage Embryo",Extant_classifier$Age)]

Embryo_CrocDorsal.pca <- Extant_CrocDorsal.pca
Embryo_CrocDorsal.pca$x <- Extant_CrocDorsal.pca$x[grep("Late-Stage Embryo|Mid-Stage Embryo",Extant_classifier$Age),]

Embryo_CrocDorsal.classifier <- Extant_classifier[grep("Late-Stage Embryo|Mid-Stage Embryo", Extant_classifier$Age),]
#
##