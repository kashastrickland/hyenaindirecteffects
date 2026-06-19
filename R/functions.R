## function to make pedigree into indexes
factorisePed <- function(pedigree, unknown=0){
  new_ped <- data.frame(
    1:nrow(pedigree), 
    ifelse(is.na(pedigree[,2]),unknown,match(pedigree[,2], pedigree[,1])), 
    ifelse(is.na(pedigree[,3]),unknown,match(pedigree[,3], pedigree[,1]))
    )
  colnames(new_ped) <- colnames(pedigree)[1:3]

  return(new_ped)
}

## function to make different categories of individual based on whether they have parents and offspring
divide_ped <- function(ped){
	out <- list(
		NoParentsOffspring = which(ped[,2]==0 & ped[,3]==0 & ped[,1] %in% c(ped[,2],ped[,3])),
		ParentsOffspring = which( (ped[,2]!=0 | ped[,3]!=0) & ped[,1] %in% c(ped[,2],ped[,3])),
		NoParentsNoOffspring = which(ped[,2]==0 & ped[,3]==0 & !ped[,1] %in% c(ped[,2],ped[,3])),
		ParentsNoOffspring = which((ped[,2]!=0 | ped[,3]!=0) & !ped[,1] %in% c(ped[,2],ped[,3])))
	
		append(out, list(
			NoParents = c(out$NoParentsOffspring,out$NoParentsNoOffspring),
			NoOffspring = c(out$NoParentsNoOffspring,out$ParentsNoOffspring)
		)
	)
}


make_stan_data <- function(ped,social_DM){


## calculate Mendelian sampling variance
MSV <- MCMCglmm::inverseA(ped[,1:3])$dii

## make pedigree into indexes rather than character strings
stan_ped <- factorisePed(ped)

##match animal ids in data to new id in ped
animal_id <- stan_ped[match(as.character(data2$ID), ped$ID),"ID"]
# mother_id <- stan_ped[match(as.character(data$dam), ped$dam),"dam"]

indirect_id <- stan_ped[match(as.character(colnames(social_DM)), ped$ID),"ID"]

## animals in pedigree that have no parents (with or without offspring), have parents and have offspring, and have parents but do not have offspring
parents <- divide_ped(stan_ped)

## make the first row in the pedigree the 'unknown parent' - bit of a hack to make looping over the pedigree faster
stan_ped2 <- rbind(c(0,-1,-1), stan_ped)+1


list(
	N_animal = nrow(stan_ped2), 
	animal_id =animal_id+1, 
	dam =stan_ped2$dam ,
	sire =stan_ped2$sire, 
	MSV=c(1,MSV), 
	N_NoParents=length(parents$NoParents), 
	N_ParentsOffspring=length(parents$ParentsOffspring), 
	N_ParentsNoOffspring=length(parents$ParentsNoOffspring), 
	NoParents=parents$NoParents+1, 
	ParentsOffspring=parents$ParentsOffspring+1,
	ParentsNoOffspring=parents$ParentsNoOffspring+1,

	N_indirect = ncol(social_DM),# width of dm
	indirect_id = indirect_id,
	indirect_DM = social_DM
)
}

make_stan_data_ks <- function(ped, social_DM, hro) {
  
  # 1. Mendelian sampling variance
  MSV <- MCMCglmm::inverseA(ped[,1:3])$dii
  
  # 2. Factorize pedigree
  stan_ped <- factorisePed(ped)
  
  # 3. Match animal IDs in data
  animal_id <- stan_ped[match(as.character(data$ID), ped$ID), "ID"]
  
  # 4. Match indirect IDs
  indirect_id <- stan_ped[match(as.character(colnames(social_DM)), ped$ID), "ID"]
  
  # 5. Match HRO IDs
  hro_id <- stan_ped[match(as.character(colnames(hro)), ped$ID), "ID"]
  
  # 6. Divide pedigree
  parents <- divide_ped(stan_ped)
  
  # 7. Prepend dummy parent (same style as before)
  stan_ped2 <- rbind(c(0, -1, -1), stan_ped) + 1
  
  # 8. Return list including HRO
  list(
    N_animal = nrow(stan_ped2),
    animal_id = animal_id + 1,
    dam = stan_ped2$dam,
    sire = stan_ped2$sire,
    MSV = c(1, MSV),
    N_NoParents = length(parents$NoParents),
    N_ParentsOffspring = length(parents$ParentsOffspring),
    N_ParentsNoOffspring = length(parents$ParentsNoOffspring),
    NoParents = parents$NoParents + 1,
    ParentsOffspring = parents$ParentsOffspring + 1,
    ParentsNoOffspring = parents$ParentsNoOffspring + 1,
    
    # Indirect effects
    N_indirect = ncol(social_DM),
    indirect_id = indirect_id,
    indirect_DM = social_DM,
    
    # HRO
    N_hro = ncol(hro),
    hro_id = hro_id,
    hro = hro
  )
}