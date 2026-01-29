
library(dplyr)


load("FINAL_sample_ids.RData")


totals <- readRDS("<path>/totals_merged.rds")
totals <- totals[totals$TE_type == 'ALU',]
totals <- totals[!is.na(totals$pop),]

totals<- totals[totals$study_code != '1KG',]



get_matched_pops <- function(cancer_type,totals) {
  # 1) Subset and stack treated vs controls
  
  
  if(cancer_type == 'all') {
    df_cases    <- totals[totals$sample_id %in% all_cancer_ids, ]
    df_controls <- totals[totals$sample_id %in% all_unrelated_names_noncancer, ]
  } else if(cancer_type == 'solid') {
    df_cases    <- totals[totals$sample_id %in% solid_ids, ]
    df_controls <- totals[totals$sample_id %in% all_unrelated_names_noncancer, ]
  } else if(cancer_type == 'liquid') {
    df_cases    <- totals[totals$sample_id %in% liquid_ids, ]
    df_controls <- totals[totals$sample_id %in% all_unrelated_names_noncancer, ]
  }
  
  tol <- 0.01  # 5% tolerance
  
  # compute case counts & control capacities by pop
  case_counts    <- table(df_cases$group)
  control_counts <- table(df_controls$group)
  pop_levels     <- names(case_counts)
  
  # target case proportions
  p_case <- case_counts / sum(case_counts)
  
  #––– 1) WATER-FILLING –––
  # find the max scale factor f so that floor(f * case_counts) ≤ control_counts
  factors <- control_counts[pop_levels] / case_counts[pop_levels]
  f       <- min(factors)
  
  # how many to take in each stratum initially
  init_counts <- floor(f * case_counts)
  
  # sample that many from each stratum
  set.seed(42)
  selected <- unlist(lapply(pop_levels, function(pop) {
    ids <- df_controls$sample_id[df_controls$group == pop]
    sample(ids, init_counts[pop])
  }))
  
  # remaining pool
  remaining <- setdiff(df_controls$sample_id, selected)
  
  
  #––– 2) FIRST-FIT GREEDY ADDITION UP TO ±5% –––
  # shuffle remaining to avoid order bias
  remaining <- sample(remaining)
  
  for (cid in remaining) {
    cand   <- c(selected, cid)
    # compute new pop proportions
    pops_sel <- df_controls$group[match(cand, df_controls$sample_id)]
    ctrl_tab  <- table(factor(pops_sel, levels = pop_levels))
    p_ctrl    <- ctrl_tab / sum(ctrl_tab)
    
    # if still within ±tol, keep it
    if (max(abs(p_ctrl - p_case)) <= tol) {
      selected <- cand
    }
    # else skip and move on
  }
  selected
}


all_unrelated_names_noncancer <- all_unrelated_names_noncancer[!all_unrelated_names_noncancer %in% all_cancer_ids]
totals$group <- paste0(totals$pop,totals$sex)

all_noncancer_pop_sex <- get_matched_pops('all',totals)
solid_noncancer_pop_sex <- get_matched_pops('solid',totals)
liquid_noncancer_pop_sex <- get_matched_pops('liquid',totals)


totals$group <- totals$pop

all_noncancer_pop <- get_matched_pops('all',totals)
solid_noncancer_pop <- get_matched_pops('solid',totals)
liquid_noncancer_pop <- get_matched_pops('liquid',totals)

table(totals[totals$sample_id %in% all_cancer_ids,]$pop)/length(all_cancer_ids)
table(totals[totals$sample_id %in% all_noncancer_pop_sex,]$pop)/length(all_noncancer_pop_sex)

table(totals[totals$sample_id %in% solid_ids,]$pop)/length(solid_ids)
table(totals[totals$sample_id %in% solid_noncancer_pop_sex,]$pop)/length(solid_noncancer_pop_sex)


table(totals[totals$sample_id %in% liquid_ids,]$pop)/length(liquid_ids)
table(totals[totals$sample_id %in% liquid_noncancer_pop_sex,]$pop)/length(liquid_noncancer_pop_sex)

save(all_noncancer_pop_sex,all_noncancer_pop,
     solid_noncancer_pop_sex,solid_noncancer_pop,
     liquid_noncancer_pop_sex,liquid_noncancer_pop,file = "NEW_ancestry_matched_nokg.RData")

