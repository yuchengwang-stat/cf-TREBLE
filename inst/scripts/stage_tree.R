#!/usr/bin/env Rscript
## Stage 1 -- build and validate the tree bundle every later stage reads.
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a   <- args_parse(list(config = cft_default_config(), from = "lists"))
cfg <- load_config(a$config)
dir <- run_dir(cfg, "tree")

tr <- if (identical(a$from, "lists")) {
  fl <- cfg$tree$from_lists
  log_msg("building tree from ", basename(fl$children))
  tree_from_lists(fl$children, fl$parent, fl$celltypes)
} else if (identical(a$from, "ref")) {
  fr <- cfg$tree$from_reference
  sl <- readRDS(cfg$data$sample_list)
  ch <- as.integer(unlist(fr$chunks %||% cfg$data$chunks))
  log_msg("building tree from the reference counts, chunks ", paste(ch, collapse = ","))
  h <- tree_from_reference(cfg$data$m_chunks, cfg$data$n_chunks, ch, sl,
                           n_cpg = as.integer(fr$n_cpg %||% 20000L),
                           min_depth = as.numeric(fr$min_depth %||% 10),
                           min_sd = as.numeric(fr$min_sd %||% 0.05),
                           dist_method = fr$dist %||% "euclidean",
                           linkage = fr$linkage %||% "ward.D2",
                           seed = as.integer(fr$seed %||% 1L))
  saveRDS(h$hclust, file.path(dir, "hclust.rds"))
  saveRDS(h$profile, file.path(dir, "celltype_profile.rds"))
  hclust_to_tree(h$hclust, names(sl))
} else if (identical(a$from, "bed")) {
  fb <- cfg$tree$from_bed
  log_msg("building tree from ", fb$bed)
  bed <- readRDS(fb$bed)
  ct  <- names(readRDS(cfg$data$sample_list))
  h   <- tree_from_bed(bed, ct, fb$knn_k, fb$dist, fb$linkage)
  saveRDS(h$hclust, file.path(dir, "hclust.rds"))
  hclust_to_tree(h$hclust, ct)
} else die("--from must be one of: lists, ref, bed")

save_rds(tr, file.path(dir, "tree.rds"))
write_manifest(dir, "tree",
  entries = list(list(path = file.path(dir, "tree.rds"))),
  extra = list(n_leaf = tr$n_leaf, n_node = tr$n_node, n_col = tr$n_col,
               n_layer = length(tr$layer), celltypes = tr$celltypes))
cat(sprintf("tree: %d leaves, %d nodes, %d likelihood columns, %d layers\n",
            tr$n_leaf, tr$n_node, tr$n_col, length(tr$layer)))
