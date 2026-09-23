#!/usr/bin/env Rscript
## Stage 5 -- merge per-chunk markers into one signature.
##
## Two things happen here that cannot happen per chunk:
##
##  * Per-label caps.  Marker counts are very uneven across cell types, and a cap
##    has to be applied after every chunk is in, or "the top N" would depend on
##    how much of the genome was scanned.
##  * The target list.  A CpG can resolve several blocks at once, which marker
##    selection records as one row per block.
##    The exported signature carries one row per CpG with every block it
##    resolves in `target`, so "this CpG separates 34 from 35" is legible
##    instead of being spread over two rows.
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

## ---- per-label caps ---------------------------------------------------------
cap_ovr  <- as.integer(if (nzchar(a$max_per_celltype)) a$max_per_celltype
                       else (cfg$markers$max_per_celltype %||% 0L))
cap_node <- as.integer(if (nzchar(a$max_per_node)) a$max_per_node
                       else (cfg$markers$max_per_node %||% 0L))
apply_cap <- function(d, cap, what) {
  if (!cap || !nrow(d)) return(d)
  before <- table(d$celltype)
  ## Rank by cross-cell-type median IQR, smallest first -- keep the markers whose
  ## methylation varies least within a cell type.  Ties broken by CpG index so
  ## the choice is reproducible.
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

## ---- one row per CpG, every block it resolves listed -------------------------
by_cpg <- split(seq_len(nrow(sig)), sig$index)
first  <- vapply(by_cpg, `[`, 0L, 1L)
## Label every target the same way regardless of which scan found it: a leaf by
## its cell type name, an internal node as node<k>.  A bare number would not say
## which of the two it is.  `sig$celltype` keeps the raw label the scans produced.
label_of <- function(nd)
  ifelse(nd <= tree$n_leaf, tree$celltypes[nd], paste0("node", nd))
targets <- vapply(by_cpg, function(r)
  paste(label_of(sig$node[r]), collapse = ";"), "")
## Four kinds, and every marker is exactly one of them:
##   celltype     one cell type (a leaf) against all the others
##   class        one internal node against all the others
##   semi_pair    a node and its sibling are each resolved as a block
##   semi_single  a node is resolved as a block; its sibling need not be
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

## A one-vs-rest marker has a single target; a semi one can have several, so the
## two are written to separate files as well as to the combined one.
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
