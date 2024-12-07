##Calculate phylogenetic signal##

##Load in phylogeny##
revised_Pseudosuchia_dichotomous.phy <- multi2di(collapse.singles(revised_Pseudosuchia.phy), random = FALSE)
revised_Pseudosuchia_calibrated.phy <- paleotree::timePaleoPhy(revised_Pseudosuchia_dichotomous.phy, Revised_Species_Duration_Table, type="equal", vartime = 0.11)

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

# create and array of mean PC scores for each species
Adult_mean_PCs <- matrix(data=NA,nrow=length(Unique_Adult_Species_list), ncol=24, dimnames = list(Unique_Adult_Species_list,colnames(Adult_CrocDorsal.pca$x)))

for (i in 1:length(Unique_Adult_Species_list)){
  current_sp <- Unique_Adult_Species_list[[i]]
  current_specimen <- grep(current_sp, Adult_CrocDorsal.classifier$Species)
  if (length(current_specimen) == 1){
    current_meanshape <- Adult_CrocDorsal.pca$x[current_specimen,]
  } else{
    current_meanshape <- mshape(Adult_CrocDorsal.pca$x[current_specimen,])
  }
  Adult_mean_PCs[i,] <- current_meanshape
}
##

## Perform phylogenetic signal calculation for basic tree ##
temp <- Adult_mean_coords[,,match(revised_Pseudosuchia_calibrated.phy$tip.label,gsub(" ", "_", dimnames(Adult_mean_coords)[[3]]))]
dimnames(temp)[[3]] <- revised_Pseudosuchia_calibrated.phy$tip.label

physignal(A=temp, phy=revised_Pseudosuchia_calibrated.phy)

temp_pc <- Adult_mean_PCs[match(revised_Pseudosuchia_calibrated.phy$tip.label,gsub(" ", "_", rownames(Adult_mean_PCs))),]
rownames(temp_pc) <- revised_Pseudosuchia_calibrated.phy$tip.label

physignal(A=temp_pc[,1:4], phy=revised_Pseudosuchia_calibrated.phy)
physignal(A=temp_pc, phy=revised_Pseudosuchia_calibrated.phy)

PC_physignal <- apply(temp_pc, 2, FUN = physignal, phy=revised_Pseudosuchia_calibrated.phy)
PC_physignal_matrix <- matrix(data=NA,nrow=length(PC_physignal), ncol=3, dimnames = list(names(PC_physignal),c("phy.signal","pvalue","Z")))

for (i in 1:nrow(PC_physignal_matrix)){
  PC_physignal_matrix[i,"phy.signal"] <- PC_physignal[[i]][["phy.signal"]]
  PC_physignal_matrix[i,"pvalue"] <- PC_physignal[[i]][["pvalue"]]
  PC_physignal_matrix[i,"Z"] <- PC_physignal[[i]][["Z"]]
}

##

## Perform phylogenetic signal calculation for popluation of cal3TimePaleoPhy trees ##
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

apply(physig_matrix, MARGIN = 2, FUN = mean)




###Code to cut?

final_MeanEmbryo.phy <- revised_Pseudosuchia_calibrated.phy
toMatch <- c("Gavialis_gangeticus","Alligator_mississippiensis","Crocodylus_porosus")

Crown_tip_list <- grep(paste(toMatch, collapse = "|"),final_MeanEmbryo.phy$tip.label)
Crown_node <- phangorn::mrca.phylo(final_MeanEmbryo.phy,Crown_tip_list)
final_MeanEmbryo.phy <- bind.tip(final_MeanEmbryo.phy,
                                 tip.label = "Mean_NonGavialid_Embryo",
                                 edge.length = 0,
                                 where = Crown_node,
                                 position = 0)

final_MeanEmbryo.phy <- bind.tip(final_MeanEmbryo.phy,
                                 tip.label = "Mean_Gavialid_Embryo",
                                 edge.length = 0,
                                 where = grep("Mean_NonGavialid_Embryo", final_MeanEmbryo.phy$tip.label),
                                 position = 0)
final_MeanEmbryo.phy <- multi2di(final_MeanEmbryo.phy, random = FALSE)
fix_zeros <- 0.01*min(final_MeanEmbryo.phy$edge.length[final_MeanEmbryo.phy$edge.length>0])
final_MeanEmbryo.phy$edge.length[final_MeanEmbryo.phy$edge.length == 0] <- fix_zeros


