#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "boot.R"))
a <- args_parse(list(config = cft_default_config(),
                     max_per_celltype = "", max_per_node = ""))
cfg  <- load_config(a$config)
tree <- readRDS(file.path(run_dir(cfg, "tree", FALSE), "tree.rds"))
mdir <- run_dir(cfg, "markers", FALSE)
dir  <- run_dir(cfg, "signature")

fs <- sort(list.files(mdir, pattern = "^markers_c[0-9]+\\.rds$", full.names = TRUE))
if (!length(fs)) die("no marker output in ", mdir)
parts <- lapply(fs, readRDS)
sig <- do.call(rbind, Filter(Negate(is.null), lapply(parts, `[[`, "markers")))
if (is.null(sig) || !nrow(sig)) die("no markers survived selection")
log_msg("merged ", length(fs), " chunks -> ", nrow(sig), " rows over ",
        length(unique(sig$index)), " CpGs")

cap_ovr  <- as.integer(if (nzchar(a$max_per_celltype)) a$max_per_celltype
                       else (cfg$markers$max_per_celltype %||% 0L))
cap_node <- as.integer(if (nzchar(a$max_per_node)) a$max_per_node
                       else (cfg$markers$max_per_node %||% 0L))
apply_cap <- function(d, cap, what) {
  if (!cap || !nrow(d)) return(d)
  before <- table(d$celltype)
  d <- d[order(d$celltype, d$median_iqr, d$index), , drop = FALSE]
  keep <- unlist(lapply(split(seq_len(nrow(d)), d$celltype), head, cap), use.names = FALSE)
  d <- d[sort(keep), , drop = FALSE]
  hit <- names(before)[before > cap]
  log_msg(sprintf("  %s cap %d: %d of %d labels were over it (%s), %d -> %d rows",
                  what, cap, length(hit), length(before),
                  paste(head(hit, 4), collapse = ","), sum(before), nrow(d)))
  d
}
sig <- rbind(apply_cap(sig[sig$onevsrest == 1L, , drop = FALSE], cap_ovr,  "celltype/class"),
             apply_cap(sig[sig$onevsrest == 0L, , drop = FALSE], cap_node, "semi"))
sig <- sig[order(sig$index, sig$node), , drop = FALSE]

by_cpg <- split(seq_len(nrow(sig)), sig$index)
first  <- vapply(by_cpg, `[`, 0L, 1L)
label_of <- function(nd)
  ifelse(nd <= tree$n_leaf, tree$celltypes[nd], paste0("node", nd))
targets <- vapply(by_cpg, function(r)
  paste(label_of(sig$node[r]), collapse = ";"), "")
sp <- if (is.null(sig$semi_pair)) rep(NA, nrow(sig)) else sig$semi_pair
kind_of <- ifelse(sig$onevsrest == 0L,
                  ifelse(is.na(sp), "semi", ifelse(sp, "semi_pair", "semi_single")),
                  ifelse(sig$node <= tree$n_leaf, "celltype", "class"))
kinds <- vapply(by_cpg, function(r) paste(kind_of[r], collapse = ";"), "")
wide <- sig[first, , drop = FALSE]
wide$target   <- unname(targets)
wide$kind     <- unname(kinds)
wide$n_target <- lengths(by_cpg)
wide$celltype <- NULL; wide$node <- NULL; wide$onevsrest <- NULL; wide$score <- NULL
wide$semi_pair <- NULL
log_msg(nrow(wide), " CpGs; ", sum(wide$n_target > 1L), " resolve more than one block")

beta <- signature_beta(wide, tree$celltypes)
meta <- wide[, c("index", "target", "kind", "n_target")]

is_ovr <- meta$kind %in% c("celltype", "class")
for (k in c("celltype", "class", "semi", "semi_pair", "semi_single"))
  cat(sprintf("  %-11s %6d CpGs\n", k, sum(grepl(k, meta$kind, fixed = TRUE))))

save_rds(list(table = wide, long = sig, beta = beta, celltypes = tree$celltypes,
              caps = c(celltype_class = cap_ovr, semi = cap_node),
              n_scanned = sum(vapply(parts, `[[`, 0L, "n_scanned"))),
         file.path(dir, "signature.rds"))
gz <- file.path(dir, "signature.txt.gz")
data.table::fwrite(cbind(meta, beta), gz, sep = "\t", compress = "gzip")
for (nm in c("one_vs_rest", "semi")) {
  sel <- if (nm == "one_vs_rest") is_ovr else !is_ovr
  f <- file.path(dir, sprintf("signature_%s.txt.gz", nm))
  data.table::fwrite(cbind(meta[sel, , drop = FALSE], beta[sel, , drop = FALSE]),
                     f, sep = "\t", compress = "gzip")
  log_msg("wrote ", basename(f), " (", sum(sel), " CpGs)")
}
log_msg("wrote ", gz)
write_manifest(dir, "signature",
  entries = list(list(path = file.path(dir, "signature.rds")), list(path = gz)),
  extra = list(n_cpg = nrow(wide), n_rows = nrow(sig),
               n_celltype = sum(grepl("celltype", wide$kind, fixed = TRUE)),
               n_class    = sum(grepl("class", wide$kind, fixed = TRUE)),
               n_semi     = sum(grepl("semi", wide$kind, fixed = TRUE)),
               n_semi_pair   = sum(grepl("semi_pair", wide$kind, fixed = TRUE)),
               n_semi_single = sum(grepl("semi_single", wide$kind, fixed = TRUE)),
               n_multi_target = sum(wide$n_target > 1L),
               caps = list(celltype_class = cap_ovr, semi = cap_node),
               n_scanned = sum(vapply(parts, `[[`, 0L, "n_scanned")),
               celltypes = tree$celltypes))
cat(sprintf("signature: %d CpGs (%d rows) x %d cell types -> %s\n",
            nrow(wide), nrow(sig), length(tree$celltypes), gz))
