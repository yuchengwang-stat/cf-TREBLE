library(ggtree)
library(ggplot2)
library(ape)
library(impute)        # for impute.knn
library(rtracklayer)
library(GenomicRanges)
#' BuildClusterTree
#' @description Estimate the heirachical tree
#' @export
BuildClusterTree <- function(bed,
                             k = 5,
                             dist_function = "euclidean",
                             linkage_function = "ward.D2")
{
  # Split metadata and score matrix
  metadata <- bed[, 1:3, drop = FALSE]
  score_data <- bed[, -(1:3), drop = FALSE]

  # Filter rows with too many NAs (<= 50% missing kept)
  na_percentage <- rowMeans(is.na(score_data))
  rows_to_keep <- na_percentage <= 0.5
  metadata_filtered <- metadata[rows_to_keep, , drop = FALSE]
  score_data_filtered <- score_data[rows_to_keep, , drop = FALSE]

  # impute.knn expects features in rows, samples in cols -> transpose
  score_matrix <- t(as.matrix(score_data_filtered))

  # k-NN imputation (use the function argument k)
  imputed <- impute.knn(score_matrix, k = k)
  imputed_matrix <- imputed$data

  # Back to original orientation: rows = loci, cols = samples
  imputed_scores <- t(imputed_matrix)
  imputed_scores_df <- as.data.frame(imputed_scores)
  colnames(imputed_scores_df) <- colnames(score_data)

  # Recombine with metadata
  final_df <- cbind(metadata_filtered, imputed_scores_df)

  # Distance & clustering (honor function arguments)
  score_mat_for_dist <- t(as.matrix(final_df[, -(1:3), drop = FALSE]))  # samples x features
  dist_matrix <- dist(score_mat_for_dist, method = dist_function)
  hc <- hclust(dist_matrix, method = linkage_function)

  # Convert to phylo and plot a fan tree
  phylo_tree <- as.phylo(hc)
  op <- par(cex = 1.4, mar = c(10, 10, 10, 1))
  on.exit(par(op), add = TRUE)

  p <- ggtree(phylo_tree, layout = "fan") +
    geom_tiplab(size = 3) +
    ggtitle("Clustering tree")
  print(p)

  return(hc)
}
