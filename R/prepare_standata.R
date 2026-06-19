
## ------------------------------------------------------------
## Import data
## ------------------------------------------------------------

data <- read.csv("data/lifehistory_clans.csv", stringsAsFactors = FALSE)
ped  <- read.csv("data/hyenapedigree.csv", stringsAsFactors = FALSE)

## Make sure IDs are character
data$ID     <- as.character(data$ID)
data$mother <- as.character(data$mother)

ped$ID   <- as.character(ped$ID)
ped$dam  <- as.character(ped$dam)
ped$sire <- as.character(ped$sire)

## ------------------------------------------------------------
## Load social network
## ------------------------------------------------------------

social <- as.matrix(read.csv("data/SRI_temporal_5.csv", row.names = 1))
colnames(social) <- rownames(social)

## Ensure names are character
rownames(social) <- as.character(rownames(social))
colnames(social) <- as.character(colnames(social))

## Remove self effects
diag(social) <- 0

## ------------------------------------------------------------
## Load social similarity matrix
## ------------------------------------------------------------

load("data/SRItmp_sim.jac.R")

## Assuming the loaded object is called SRItmp_sim.jac
socsim <- SRItmp_sim.jac

rownames(socsim) <- as.character(rownames(socsim))
colnames(socsim) <- as.character(colnames(socsim))

## ------------------------------------------------------------
## Load spatial overlap matrix
## ------------------------------------------------------------

hro <- as.matrix(read.csv("data/hro_VI_5.csv", row.names = 1))
colnames(hro) <- rownames(hro)

rownames(hro) <- as.character(rownames(hro))
colnames(hro) <- as.character(colnames(hro))

## ------------------------------------------------------------
## Initial filtering of phenotype data
## ------------------------------------------------------------

## Keep individuals present in both the social network and HRO matrix
data2 <- data[
  data$ID %in% rownames(social) &
    data$ID %in% rownames(socsim) &
    data$ID %in% rownames(hro),
]

## ------------------------------------------------------------
## Preliminary social design matrix
## Used only to identify indirect-effect individuals for pruning pedigree
## ------------------------------------------------------------

social_DM_pre <- social[as.character(data2$ID), , drop = FALSE]

## Drop empty columns
social_DM_pre <- social_DM_pre[, colSums(social_DM_pre) != 0, drop = FALSE]

dim(social_DM_pre)

## ------------------------------------------------------------
## Prune pedigree
## ------------------------------------------------------------

## Keep all individuals who appear as indirect social partners.
keep_ids <- unique(c(
  colnames(social_DM_pre),
  as.character(data2$ID),
  as.character(data2$mother)
))

ped_IGE <- MCMCglmm::prunePed(
  ped,
  keep = keep_ids,
  make.base = TRUE
)

ped_IGE$ID   <- as.character(ped_IGE$ID)
ped_IGE$dam  <- as.character(ped_IGE$dam)
ped_IGE$sire <- as.character(ped_IGE$sire)

## ------------------------------------------------------------
## Final filtering after pedigree pruning
## ------------------------------------------------------------

## Keep observations whose focal ID and mother are in the pruned pedigree
data2 <- data2[
  data2$ID %in% ped_IGE$ID &
    data2$mother %in% ped_IGE$ID,
]

## reset rownames
rownames(data2) <- NULL

## ------------------------------------------------------------
## Final social design matrix
## ------------------------------------------------------------

## match order of rows of social_DM to match data2 rows
social_DM <- social[as.character(data2$ID), , drop = FALSE]

## Drop columns that are now empty after final data filtering.
social_DM <- social_DM[, colSums(social_DM) != 0, drop = FALSE]

## Keep only columns for individuals present in pruned pedigree
social_DM <- social_DM[, colnames(social_DM) %in% ped_IGE$ID, drop = FALSE]

## Drop empty columns again, in case the pedigree filter removed some
social_DM <- social_DM[, colSums(social_DM) != 0, drop = FALSE]

## Check alignment
stopifnot(nrow(social_DM) == nrow(data2))
stopifnot(identical(rownames(social_DM), as.character(data2$ID)))

## ------------------------------------------------------------
## HRO matrix and indexing
## ------------------------------------------------------------

## Use unique ID levels in case hro is a covariance matrix over individuals,
## not over observations.
hro_levels <- unique(as.character(data2$ID))

hro2 <- hro[hro_levels, hro_levels, drop = FALSE]

data2$hro_id <- match(as.character(data2$ID), hro_levels)

## If needed, add a small jitter to diagonal if not positive definite.
if (!corpcor::is.positive.definite(hro2)) {
  message("hro2 is not positive definite; adding small diagonal jitter.")
  diag(hro2) <- diag(hro2) + 1e-6
}

stopifnot(corpcor::is.positive.definite(hro2))
stopifnot(!anyNA(data2$hro_id))
stopifnot(all(data2$hro_id >= 1))
stopifnot(all(data2$hro_id <= nrow(hro2)))

## ------------------------------------------------------------
## Social similarity matrix and indexing
## ------------------------------------------------------------

## If social similarity is an individual-level covariance matrix, use unique
## focal IDs, same logic as HRO.
ss_levels <- unique(as.character(data2$ID))

socsim2 <- socsim[ss_levels, ss_levels, drop = FALSE]

data2$ss_id <- match(as.character(data2$ID), ss_levels)

## If not positive definite, add a small jitter.
if (!corpcor::is.positive.definite(socsim2)) {
  message("socsim2 is not positive definite; adding small diagonal jitter.")
  diag(socsim2) <- diag(socsim2) + 1e-6
}

stopifnot(corpcor::is.positive.definite(socsim2))
stopifnot(!anyNA(data2$ss_id))
stopifnot(all(data2$ss_id >= 1))
stopifnot(all(data2$ss_id <= nrow(socsim2)))

## ------------------------------------------------------------
## Fixed-effect design matrix
## ------------------------------------------------------------

##add in an effect of total strength
data2$social_strength <- rowSums(social_DM)
data2$degree <- rowSums(social_DM>0)

##create fixed effects matrix
X <- model.matrix(
  LRS ~ sex + mother_rank_ord + clan_size + f + X + degree,
  data = data2
)

## Check response
Y <- as.integer(data2$LRS)
stopifnot(all(Y >= 0))
stopifnot(nrow(X) == nrow(data2))

## ------------------------------------------------------------
## Random-effect index variables
## ------------------------------------------------------------

data2$mother_id <- as.integer(factor(data2$mother))
data2$by_id     <- as.integer(factor(data2$birth_year))
data2$clan_id   <- as.integer(factor(data2$birth_clan))

stopifnot(!anyNA(data2$mother_id))
stopifnot(!anyNA(data2$by_id))
stopifnot(!anyNA(data2$clan_id))

## ------------------------------------------------------------
## Pedigree indexing
## ------------------------------------------------------------

## Mendelian sampling variance.
MSV <- MCMCglmm::inverseA(ped_IGE[, 1:3])$dii

## Convert pedigree IDs to integer indexing.
stan_ped <- factorisePed(ped_IGE)

stopifnot(all(c("ID", "dam", "sire") %in% colnames(stan_ped)))

## Divide pedigree into groups needed by Stan
parents <- divide_ped(stan_ped)

## Add dummy parent row as row 1.
stan_ped2 <- rbind(c(0, -1, -1), stan_ped) + 1
stan_ped2 <- as.data.frame(stan_ped2)
colnames(stan_ped2) <- colnames(stan_ped)

stopifnot(all(stan_ped2$ID == seq_len(nrow(stan_ped2))))
stopifnot(all(stan_ped2$dam >= 0))
stopifnot(all(stan_ped2$sire >= 0))
stopifnot(all(stan_ped2$dam <= nrow(stan_ped2)))
stopifnot(all(stan_ped2$sire <= nrow(stan_ped2)))

## Animal IDs for focal observations.
animal_id_unshifted <- stan_ped[
  match(as.character(data2$ID), ped_IGE$ID),
  "ID"
]

## Animal IDs for indirect effect columns.
indirect_id_unshifted <- stan_ped[
  match(as.character(colnames(social_DM)), ped_IGE$ID),
  "ID"
]

stopifnot(!anyNA(animal_id_unshifted))
stopifnot(!anyNA(indirect_id_unshifted))

## Shift by +1 for Stan, because of dummy row.
animal_id   <- animal_id_unshifted + 1
indirect_id <- indirect_id_unshifted + 1

stopifnot(all(animal_id >= 1))
stopifnot(all(animal_id <= nrow(stan_ped2)))
stopifnot(all(indirect_id >= 1))
stopifnot(all(indirect_id <= nrow(stan_ped2)))

## Pedigree categories shifted by +1.
NoParents            <- parents$NoParents + 1
ParentsOffspring     <- parents$ParentsOffspring + 1
ParentsNoOffspring   <- parents$ParentsNoOffspring + 1

## Check that the three pedigree groups cover all non-dummy animals exactly once.
all_ped_groups <- c(NoParents, ParentsOffspring, ParentsNoOffspring)

stopifnot(!anyDuplicated(all_ped_groups))
stopifnot(setequal(all_ped_groups, 2:nrow(stan_ped2)))

## Full MSV including dummy parent.
MSV_stan <- c(1, MSV)

stopifnot(length(MSV_stan) == nrow(stan_ped2))
stopifnot(all(MSV_stan >= 0))

## ------------------------------------------------------------
## Final Stan data list
## ------------------------------------------------------------

stan_dat_IGE_NB_HRO <- list(
  ## Observation data
  N = nrow(data2),
  Y = Y,
  K = ncol(X),
  X = X,
  
  ## Pedigree/animal effect data
  N_animal = nrow(stan_ped2),
  animal_id = as.integer(animal_id),
  dam = as.integer(stan_ped2$dam),
  sire = as.integer(stan_ped2$sire),
  MSV = as.vector(MSV_stan),
  
  N_NoParents = length(NoParents),
  N_ParentsOffspring = length(ParentsOffspring),
  N_ParentsNoOffspring = length(ParentsNoOffspring),
  
  NoParents = as.integer(NoParents),
  ParentsOffspring = as.integer(ParentsOffspring),
  ParentsNoOffspring = as.integer(ParentsNoOffspring),
  
  ## Indirect effects via social network design matrix
  N_indirect = ncol(social_DM),
  indirect_id = as.integer(indirect_id),
  indirect_DM = as.matrix(social_DM),
  
  ## Mother random effect
  N_mother = length(unique(data2$mother_id)),
  mother_id = as.integer(data2$mother_id),
  
  ## Birth-year random effect
  N_by = length(unique(data2$by_id)),
  by_id = as.integer(data2$by_id),
  
  ## Birth-clan random effect
  N_clan = length(unique(data2$clan_id)),
  clan_id = as.integer(data2$clan_id),
  
  ## Spatial overlap matrix
  N_hro = nrow(hro2),
  hro_id = as.integer(data2$hro_id),
  hro = as.matrix(hro2),
  
  ## Social similarity matrix
  N_ss = nrow(socsim2),
  ss_id = as.integer(data2$ss_id),
  socsim = as.matrix(socsim2)
)

## ------------------------------------------------------------
## Final checks
## ------------------------------------------------------------

with(stan_dat_IGE_NB_HRO, {
  ## Observation data
  stopifnot(length(Y) == N)
  stopifnot(nrow(X) == N)
  stopifnot(ncol(X) == K)
  
  ## Animal IDs
  stopifnot(length(animal_id) == N)
  stopifnot(all(animal_id >= 1))
  stopifnot(all(animal_id <= N_animal))
  
  ## Pedigree arrays
  stopifnot(length(dam) == N_animal)
  stopifnot(length(sire) == N_animal)
  stopifnot(length(MSV) == N_animal)
  
  stopifnot(length(NoParents) == N_NoParents)
  stopifnot(length(ParentsOffspring) == N_ParentsOffspring)
  stopifnot(length(ParentsNoOffspring) == N_ParentsNoOffspring)
  
  ## Indirect design matrix
  stopifnot(nrow(indirect_DM) == N)
  stopifnot(ncol(indirect_DM) == N_indirect)
  stopifnot(length(indirect_id) == N_indirect)
  stopifnot(all(indirect_id >= 1))
  stopifnot(all(indirect_id <= N_animal))
  
  ## Mother
  stopifnot(length(mother_id) == N)
  stopifnot(all(mother_id >= 1))
  stopifnot(all(mother_id <= N_mother))
  
  ## Birth year
  stopifnot(length(by_id) == N)
  stopifnot(all(by_id >= 1))
  stopifnot(all(by_id <= N_by))
  
  ## Clan
  stopifnot(length(clan_id) == N)
  stopifnot(all(clan_id >= 1))
  stopifnot(all(clan_id <= N_clan))
  
  ## HRO
  stopifnot(length(hro_id) == N)
  stopifnot(all(hro_id >= 1))
  stopifnot(all(hro_id <= N_hro))
  stopifnot(nrow(hro) == N_hro)
  stopifnot(ncol(hro) == N_hro)
  stopifnot(isTRUE(all.equal(hro, t(hro))))
  
  ## Social similarity
  stopifnot(length(ss_id) == N)
  stopifnot(all(ss_id >= 1))
  stopifnot(all(ss_id <= N_ss))
  stopifnot(nrow(socsim) == N_ss)
  stopifnot(ncol(socsim) == N_ss)
  stopifnot(isTRUE(all.equal(socsim, t(socsim))))
})

str(stan_dat_IGE_NB_HRO)
