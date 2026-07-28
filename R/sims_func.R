
run_one_sim_and_fit <- function(i) {
  
  tryCatch({
    
    cat("Starting simulation", i, "\n")
    
    set.seed(10000 + i)
    
    ## ----------------------------
    ## Fixed effects
    ## ----------------------------
    
    beta <- as.numeric(random_draw[i, 1:6])
    names(beta) <- colnames(X)

    ## ----------------------------
    ## HRO effects
    ## ----------------------------
    hro_effects<-t(chol(hro2)) %*% rnorm(nrow(hro2),0,random_draw[i,"sd_hro"])
    ## name by rownames
    names(hro_effects) <- rownames(stan_dat_IGE_NB_HRO$hro)

    ## ----------------------------
    ## Genetic effects
    ## ----------------------------
    Sigma_A <- matrix(c(random_draw[i, "Sigma_A[1,1]"],random_draw[i, "Sigma_A[1,2]"],
                        random_draw[i, "Sigma_A[2,1]"],random_draw[i, "Sigma_A[2,2]"]),nrow = 2,byrow = TRUE)
    genetic_effects <- MCMCglmm::rbv(ped_IGE, Sigma_A)
    genetic_effects <- as.matrix(genetic_effects)
    colnames(genetic_effects) <- c("direct_genetic", "indirect_genetic")
    rownames(genetic_effects) <- as.character(rownames(genetic_effects))
    
    ## ----------------------------
    ## Mother effects
    ## ----------------------------
    mother_ids <- unique(as.character(data2$mother))
    mother_effects <- rnorm(length(mother_ids),mean = 0,sd = sqrt(random_draw[i, "var_mother"]))
    names(mother_effects) <- mother_ids

    ## ----------------------------
    ## Cohort effects
    ## ----------------------------
    cohort_ids <- unique(as.character(data2$birth_year))
    cohort_effects <- rnorm(length(cohort_ids),mean = 0,sd = random_draw[i, "sd_by"])
    names(cohort_effects) <- cohort_ids

    ## ----------------------------
    ## Clan effects
    ## ----------------------------
    clan_ids <- unique(as.character(data2$birth_clan))
    clan_effects <- rnorm(length(clan_ids),mean = 0,sd = random_draw[i, "sd_clan"])
    names(clan_effects) <- clan_ids
    
    ## ----------------------------
    ## Indirect environmental effects
    ## ----------------------------
    indirect_effects <- rnorm(ncol(social_DM),mean = 0,sd = random_draw[i, "sd_indirectE"])
    names(indirect_effects) <- colnames(social_DM)
    
    ## Check indexing
    ## Check indexing
    id <- as.character(data2$ID)
    mother <- as.character(data2$mother)
    birth_year <- as.character(data2$birth_year)
    birth_clan <- as.character(data2$birth_clan)
    social_ids <- as.character(colnames(social_DM))
    hro_id <- as.character(data2$hro_id)
    
    ## stop if indexing not working
    stopifnot(nrow(X) == nrow(data2))
    stopifnot(nrow(social_DM) == nrow(data2))
    stopifnot(length(beta) == ncol(X))
    
    stopifnot(all(id %in% rownames(genetic_effects)))
    stopifnot(all(mother %in% names(mother_effects)))
    stopifnot(all(birth_year %in% names(cohort_effects)))
    stopifnot(all(birth_clan %in% names(clan_effects)))
    stopifnot(all(social_ids %in% rownames(genetic_effects)))
    stopifnot(all(social_ids %in% names(indirect_effects)))
    
    ## Indirect effects
    indirect_total <- indirect_effects[social_ids] +
      genetic_effects[social_ids, "indirect_genetic"]
    
    indirect_eta <- as.numeric(social_DM %*% indirect_total)
    hro_eta <- hro_effects[id]
    
    ## ----------------------------
    ## Linear predictor
    ## ----------------------------
    
    eta <- as.numeric(
      (X %*% beta) +
        genetic_effects[id, "direct_genetic"] +
        mother_effects[mother] +
        cohort_effects[birth_year] +
        clan_effects[birth_clan] +
        indirect_eta +
        hro_eta
    )
    
    ## ----------------------------
    ## Negative binomial simulation
    ## ----------------------------
    
    mu <- exp(eta)
    
    theta <- random_draw[i, "shape"]
    
    y <- rnbinom(n = nrow(data2),
      mu = mu,
      size = theta)
    
    ## ----------------------------
    ## Local Stan data object
    ## ----------------------------
    
    stan_dat_i <- stan_dat_IGE_NB_HRO
    stan_dat_i$Y <- y
    
    ## ----------------------------
    ## Fit Stan model
    ## ----------------------------
    
    fit <- rstan::sampling(
      object = stan_model_IGE_NB_VA_fullIGE,
      data = stan_dat_i,
      cores = 3,
      chains = 3,
      iter = 3500,
      warmup = 1000,
      pars = c(
        "b",
        "shape",
        "var_mother",
        "var_by",
        "var_clan",
        "var_indirectE",
        "var_hro",
        "sd_by",
        "sd_clan",
        "sd_indirectE",
        "sd_hro",
        "Sigma_A",
        "corr_A"
      ),
      seed = 20000 + i,
      refresh = 500
    )
    
    ## ----------------------------
    ## Save full model
    ## ----------------------------
    
    out <- list(
      sim_id = i,
      true_parameters = random_draw[i, ],
      y = y,
      fit = fit,
      error = FALSE
    )
    
    out_file <- file.path(
      "sim_fits_full",
      paste0("IGE_sim_full_", sprintf("%03d", i), ".rds")
    )
    
    saveRDS(out, out_file, compress = "xz")
    
    rm(fit, stan_dat_i, out)
    gc()
    
    cat("Finished simulation", i, "\n")
    
    return(out_file)
    
  }, error = function(e) {
    
    err_file <- file.path(
      "sim_fits_full",
      paste0("IGE_sim_full_", sprintf("%03d", i), "_ERROR.rds")
    )
    
    err_out <- list(
      sim_id = i,
      true_parameters = random_draw[i, ],
      error = TRUE,
      message = conditionMessage(e),
      time = Sys.time()
    )
    
    saveRDS(err_out, err_file)
    
    cat("Simulation", i, "failed:", conditionMessage(e), "\n")
    
    return(err_file)
  })
}
