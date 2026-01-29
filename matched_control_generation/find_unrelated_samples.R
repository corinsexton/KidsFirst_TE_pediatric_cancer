# plot distribution of insertions

library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(igraph)


# output of vcftools --relatedness2
relatedness <- read_tsv("related_individuals_gt0.tsv")

# boxplot(relatedness$RELATEDNESS_PHI)

relatedness_cutoff <- relatedness[relatedness$RELATEDNESS_PHI >= 0.08,]
relatedness_cutoff <- relatedness_cutoff %>%
  rowwise() %>%
  mutate(
    # Sort the two IDs
    first_sorted  = sort(c(INDV1, INDV2))[1],
    second_sorted = sort(c(INDV1, INDV2))[2]
  ) %>%
  ungroup() %>%
  transmute(
    INDV1 = first_sorted,
    INDV2 = second_sorted
  ) %>%
  distinct(INDV1, INDV2) %>%
  filter(INDV1 != INDV2)  # optionally remove self-pairs if they exist



# merged insertion file
totals_merged <- readRDS("<path>/totals_merged_newOBP.rds")
totals_merged <- totals_merged[which(totals_merged$filter == 'pass'),]

cancer_studies = c("KF-ESGR","KF-FALL","KF-IGCT","KF-MMC",
                   "KF-NBL","KF-NCSF","KF-TALL","KF-GNINT")
cancer_parent_controls_ids <-  totals_merged[totals_merged$study_code %in% cancer_studies &
                                               totals_merged$proband_status != 'proband',]$sample_id %>% unique

totals_merged = totals_merged[!totals_merged$sample_id %in% c("H1009CLAb1",
                                                              "H1009CLAc1",
                                                              "H1014HOSb1",
                                                              "H1014HOSc1"),]


liquid_ids <- totals_merged$sample_id[grepl('liquid',totals_merged$cancer_type)] %>% unique()
solid_ids <- totals_merged$sample_id[grepl('solid',totals_merged$cancer_type)] %>% unique()

length(totals_merged$sample_id[grepl('id',totals_merged$cancer_type)] %>% unique())


all_cancer_ids <- totals_merged$sample_id[grepl('id',totals_merged$cancer_type)] %>% unique()
control_non_cancer_ids <- totals_merged[totals_merged$proband_status == 'control' &
                                          !totals_merged$sample_id %in% cancer_parent_controls_ids,]$sample_id %>% unique()
control_ids <- totals_merged[totals_merged$proband_status == 'control',]$sample_id %>% unique()
all_samples <- totals_merged$sample_id %>% unique()


relatedness_cutoff <- relatedness_cutoff[relatedness_cutoff$INDV1 %in% all_samples & relatedness_cutoff$INDV2 %in% all_samples,]


# Assume your data frame is named df and has columns "INDV1" and "INDV2".
# Each row represents a relationship between INDV1 and INDV2.
# For example:
# df <- data.frame(INDV1 = c("A", "B", "C", "D"),
#                  INDV2 = c("B", "C", "D", "E"))


approxMIS_smallest_degree_with_forced <- function(g, forced_vertices) {
  # 1. Ensure forced vertices actually exist in the graph
  #    (use intersection in case 'forced_vertices' includes IDs not in g)
  forced_vertices <- intersect(forced_vertices, V(g)$name)
  print(length(forced_vertices))
  # Start the independent set with the forced vertices
  independent_set <- forced_vertices
  

  # 2. Remove forced vertices and their neighbors from the graph
  forced_neighbors <- unique(unlist(
    lapply(forced_vertices, function(fv) {
      # Convert igraph vertex sequence to character vector of names
      as_ids(neighbors(g, fv))
    })
  ))
  
  # Build a union of forced vertices + their neighbor *names*
  vertices_to_remove <- union(forced_vertices, forced_neighbors)
  
  # Now delete by vertex name
  g_temp <- delete_vertices(g, vertices_to_remove)
  
  # 3. Proceed with the standard greedy approach on the subgraph
  while (vcount(g_temp) > 0) {
    # Pick the vertex with the smallest degree in the subgraph
    v <- V(g_temp)[which.min(degree(g_temp))]
    # Add it to the independent set
    independent_set <- c(independent_set, as_ids(v))
    
    # Remove that vertex and its neighbors
    nbrs <- neighbors(g_temp, v)
    g_temp <- delete_vertices(g_temp, c(v, nbrs))
  }
  # # 4. The final set is 'independent_set'. 
  # #    To return the *excluded* set, take all vertex names minus those in the independent set.
  # all_vertices <- V(g)$name
  # excluded_set <- setdiff(all_vertices, independent_set)
  
  return(independent_set)
}


get_sample_names <- function(df,required_ids,all_samples) {
  # Create an undirected graph from the data frame.
  g <- graph_from_data_frame(df, directed = FALSE,
                             vertices = data.frame(name = all_samples, stringsAsFactors = FALSE))
  
  # Find the largest clique in the complement graph.
  # A clique in the complement graph corresponds to an independent set in the original graph.
  set_approx <- approxMIS_smallest_degree_with_forced(g,required_ids)
  
  set_approx
}




liquid_unrelated_names <- get_sample_names(relatedness_cutoff,liquid_ids,all_samples)
solid_unrelated_names <- get_sample_names(relatedness_cutoff,solid_ids,all_samples)
all_unrelated_names <- get_sample_names(relatedness_cutoff,all_cancer_ids,all_samples)

all_unrelated_names_noncancer <- all_unrelated_names[all_unrelated_names %in% c(control_non_cancer_ids,all_cancer_ids)]
liquid_unrelated_names_noncancer <- liquid_unrelated_names[liquid_unrelated_names %in% c(control_non_cancer_ids,solid_ids)]
solid_unrelated_names_noncancer <- solid_unrelated_names[solid_unrelated_names %in% c(control_non_cancer_ids,liquid_ids)]

save(liquid_unrelated_names_noncancer,
     solid_unrelated_names_noncancer,
     all_unrelated_names_noncancer,
     liquid_ids,
     solid_ids,
     all_cancer_ids,
     cancer_parent_controls_ids,
     control_ids,
     control_non_cancer_ids,
     file = "FINAL_sample_ids.RData")

