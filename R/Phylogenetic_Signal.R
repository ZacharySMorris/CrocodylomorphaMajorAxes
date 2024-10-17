##Calculate phylogenetic signal##

##Load in phylogeny##
revised_Pseudosuchia_dichotomous.phy <- multi2di(collapse.singles(revised_Pseudosuchia.phy), random = FALSE)
revised_Pseudosuchia_calibrated.phy <- paleotree::timePaleoPhy(revised_Pseudosuchia_dichotomous.phy, Revised_Species_Duration_Table, type="equal", vartime = 0.11)

temp <- Adult_mean_coords[,,match(revised_Pseudosuchia_calibrated.phy$tip.label,gsub(" ", "_", dimnames(Adult_mean_coords)[[3]]))]
dimnames(temp)[[3]] <- revised_Pseudosuchia_calibrated.phy$tip.label

physignal(A=temp, phy=revised_Pseudosuchia_calibrated.phy)

##

##Create mean adult/subadult shapes for all species##
## list of species with adults/subadults
Adult_coords <- Adult_CrocDorsal.gpa$coords
dimnames(Adult_coords)[[3]] <- Adult_CrocDorsal.classifier$Species
Unique_Adult_Species_list <- unique(Adult_CrocDorsal.classifier$Species)
##

# create an array for the mean shapes of each species
Adult_mean_coords <- array(dim = c(14,2,length(Unique_Adult_Species_list)), dimnames = list(c(1:14),c("X","Y"),Unique_Adult_Species_list))

for (i in 1:length(Unique_Adult_Species_list)){
  current_sp <- Unique_Adult_Species_list[[i]]
  current_specimen <- grep(current_sp, dimnames(Adult_coords)[[3]])
  current_meanshape <- mshape(Adult_coords[,,current_specimen])

  Adult_mean_coords[,,i] <- current_meanshape
}
##

## Perform phylogenetic signal calculation for all trees ##
physig_list <- list()
physig_matrix <- matrix(nrow=length(final_Pseudosuchia.phy),ncol=3,dimnames = list(c(1:length(final_Pseudosuchia.phy)),c("PhySig(K)","p-value","Effect Size")))
for (i in 1:length(final_Pseudosuchia.phy)){
  tree <- final_Pseudosuchia.phy[[i]]
    # Reorder coords to match phylogeny ##
    temp <- Adult_mean_coords[,,match(tree$tip.label,gsub(" ", "_", dimnames(Adult_mean_coords)[[3]]))]
    dimnames(temp)[[3]] <- tree$tip.label
    #
  physig_temp  <- physignal(A=temp, phy=tree)
  physig_matrix[i,1] <- physig_temp$phy.signal
  physig_matrix[i,2] <- physig_temp$pvalue
  physig_matrix[i,3] <- physig_temp$Z
}
#
