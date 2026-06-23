fetch_id_inbreeding <- function(ID) {
  
  inputID <- hyenaR::check_function_arg.ID(ID, argument.name = "ID", strict = TRUE, .fill = FALSE, arg.max.length = Inf)
  
  input.tbl <- data.frame(ID = inputID)
  
  if (!requireNamespace("ggroups", quietly = TRUE)) {
    stop("Please install R package 'ggroups' to use this function")
  }
  ## build a pedigree to start from
  pedigree <- hyenaR::create_pop_pedigree(pedigree.format = "both")
  
  ## select columns needed for the package ggroups
  pedigree[, c("ID_num", "sire_num", "dam_num")] -> pedigree_ggroups
  
  ## compute inbreeding
  ggroups::inbreed(ped = as.data.frame(pedigree_ggroups)) -> inbr_coeff
  
  ## format to typical hyenaR tibble::tibble
  tibble::tibble(ID = pedigree$ID, inbr_coeff) -> inbr_tbl
  
  hyenaR::check_function_output(input.tbl = input.tbl, output.tbl = inbr_tbl, join.by = "ID",
                                duplicates = "input", output.IDcolumn = "inbr_coeff")
}


fetch_id_paternity.checked <- function(ID) {
  inputID <- hyenaR::check_function_arg.ID(ID, argument.name = "ID", strict = TRUE, .fill = FALSE, arg.max.length = Inf)
  ID_microsat <- extract_database_table("alleles") |> pull(ID)
  ID %in% ID_microsat
}

## Right margin calculation based on lifespan quantiles:

calculate_lost_IDs_right <- function(comparison_data, right.censoring = NULL, sex = NULL) {
  
  if (sex == "both") {
    sex <- c("male", "female", NA)
  } else if (sex == "males") {
    sex <- "male"
  } else if (sex == "females") {
    sex <- "female"
  } else 
    stop("sex missing")
  
  ref_data <- hyenaR::create_id_starting.table(clan.birth = hyenaR::find_clan_name.all(main.clans = TRUE),
                                               sex = sex,
                                               lifestage = "cub", lifestage.overlap = "start",
                                               from = hyenaR::find_pop_date.observation.first()) |>
    dplyr::semi_join(comparison_data, by = c("ID" = "FID"))
  
  extended_data <- hyenaR::create_id_starting.table(clan.birth = hyenaR::find_clan_name.all(main.clans = TRUE),
                                                    sex = sex,
                                                    lifestage = "cub", lifestage.overlap = "start",
                                                    from = hyenaR::find_pop_date.observation.first(),
                                                    to = right.censoring) |>
    dplyr::semi_join(comparison_data, by = c("ID" = "FID"))
  
  lost_data <- dplyr::anti_join(ref_data, extended_data, by = join_by(ID)) |>
    left_join(comparison_data, by = c("ID" = "FID"))
  tibble::tibble(n_lost = nrow(lost_data), pct_lost = round(100*nrow(lost_data)/nrow(ref_data)), Fhat3_max = max(lost_data$Fhat3), ID = list(lost_data$ID)) 
  
}
