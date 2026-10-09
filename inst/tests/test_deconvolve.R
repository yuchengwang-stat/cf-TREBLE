#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])), "..", "scripts", "boot.R"))
set.seed(4)
tmp <- file.path(tempdir(), paste0("cft_dec_", Sys.getpid()))
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)

G   <- 3000L
ct  <- sprintf("CT%02d", 1:6)
beta <- matrix(rbeta(G * length(ct), 0.5, 0.5), G, length(ct), dimnames = list(NULL, ct))
truth <- c(0.40, 0.25, 0.15, 0.10, 0.07, 0.03); names(truth) <- ct
idx <- sort(sample.int(20000L, G))

nBeta <- 20000L
depth <- rep(200L, nBeta); meth <- integer(nBeta)
p <- as.vector(beta %*% truth)
meth[idx] <- rbinom(G, depth[idx], p)
bf <- file.path(tmp, "plasma.beta")
writeBin(as.integer(t(cbind(pmin(meth, 255L), pmin(depth, 255L)))), bf, size = 1L)

sigfile <- file.path(tmp, "external_signature.tsv")
write.table(cbind(index = idx, as.data.frame(beta)), sigfile,
            sep = "\t", row.names = FALSE, quote = FALSE)
S <- read_signature(sigfile)
stopifnot(identical(S$celltypes, ct), identical(S$index, idx))
est <- deconvolve_beta(bf, NULL, S$index, signature_phi(S$beta))
cat(sprintf("  external tsv    max |est - truth| = %.4f   %s\n", max(abs(est - truth)),
            if (max(abs(est - truth)) < 0.02) "PASS" else "FAIL"))
stopifnot(max(abs(est - truth)) < 0.02, abs(sum(est) - 1) < 1e-6)

withmeta <- file.path(tmp, "with_metadata.tsv")
write.table(cbind(index = idx, target = ct[1], kind = "celltype", n_target = 1L,
                  as.data.frame(beta)), withmeta, sep = "\t", row.names = FALSE, quote = FALSE)
stopifnot(identical(read_signature(withmeta)$celltypes, ct))
est2 <- deconvolve_beta(bf, NULL, S$index, signature_phi(read_signature(withmeta)$beta))
stopifnot(max(abs(est2 - est)) < 1e-10)
cat("  metadata columns ignored                        PASS\n")

sub <- read_signature(sigfile, celltypes = ct[1:4])
stopifnot(identical(sub$celltypes, ct[1:4]), ncol(sub$beta) == 4L)
cat("  --celltypes subsets the reference                PASS\n")

rep_phi <- lapply(seq_along(ct), function(k) beta[, rep(k, 3), drop = FALSE])
names(rep_phi) <- ct
dat <- read_beta(bf)[idx, , drop = FALSE]
a <- cf_BLEND(dat[, 1], dat[, 2], signature_phi(beta))$cellular_frac
b <- cf_BLEND(dat[, 1], dat[, 2], rep_phi)$cellular_frac
cat(sprintf("  replicated reference columns: max diff %.2e   %s\n", max(abs(a - b)),
            if (max(abs(a - b)) < 1e-10) "PASS" else "FAIL"))
stopifnot(max(abs(a - b)) < 1e-10)

noidx <- file.path(tmp, "no_index.tsv")
write.table(as.data.frame(beta), noidx, sep = "\t", row.names = FALSE, quote = FALSE)
e1 <- tryCatch(read_signature(noidx), error = conditionMessage)
stopifnot(grepl("index", e1))
oor <- file.path(tmp, "out_of_range.tsv")
b2 <- beta; b2[1, 1] <- 1.4
write.table(cbind(index = idx, as.data.frame(b2)), oor, sep = "\t", row.names = FALSE, quote = FALSE)
e2 <- tryCatch(read_signature(oor), error = conditionMessage)
stopifnot(grepl("\\[0, 1\\]", e2))
e3 <- tryCatch(read_signature(sigfile, celltypes = "Nope"), error = conditionMessage)
stopifnot(grepl("not in the signature", e3))
cat("  missing index / beta out of range / unknown type refused  PASS\n")

e4 <- tryCatch(deconvolve_beta(bf, NULL, c(idx, 999999L), signature_phi(beta)),
               error = conditionMessage)
stopifnot(grepl("different CpG sets", e4))
cat("  index beyond the .beta is refused                PASS\n")

sl <- setNames(split(sprintf("s%02d", seq_len(12)), rep(seq_len(6), each = 2)), ct)
subdir <- file.path(tmp, "subject")
dir.create(subdir, showWarnings = FALSE)
P <- beta[, rep(seq_along(ct), each = 2), drop = FALSE] +
     matrix(rnorm(G * 12, 0, 0.01), G, 12)
P <- pmin(pmax(P, 1e-4), 1 - 1e-4)
colnames(P) <- unlist(sl, use.names = FALSE)
half <- seq_len(G %/% 2L)
saveRDS(list(chunk = 1L, index = c(idx[half], idx[half[1]]),
             p_subj = P[c(half, half[1]), , drop = FALSE],
             samples = colnames(P)), file.path(subdir, "subject_c001.rds"))
saveRDS(list(chunk = 2L, index = idx[-half], p_subj = P[-half, , drop = FALSE],
             samples = colnames(P)), file.path(subdir, "subject_c002.rds"))

sub <- read_subject(subdir, idx)
stopifnot(nrow(sub$beta) == G, identical(sub$index, idx))
stopifnot(max(abs(sub$beta - P)) < 1e-12)
sphi <- subject_phi(sub$beta, sl)
stopifnot(identical(names(sphi), ct), all(vapply(sphi, ncol, 0L) == 2L))
cat("  subject matrix deduped, aligned, split 2 columns per cell type   PASS\n")

est3 <- cf_BLEND(dat[, 1], dat[, 2], sphi)
stopifnot(max(abs(est3$cellular_frac - truth)) < 0.03)
stopifnot(all(vapply(est3$ref_mixing_prop, length, 0L) == 2L),
          all(abs(vapply(est3$ref_mixing_prop, sum, 0) - 1) < 1e-8))
cat(sprintf("  subject-specific recovers truth: max |err| %.4f            PASS\n",
            max(abs(est3$cellular_frac - truth))))

P2 <- P; P2[5, 3] <- NA_real_
subdir2 <- file.path(tmp, "subject_na")
dir.create(subdir2, showWarnings = FALSE)
saveRDS(list(chunk = 1L, index = idx, p_subj = P2, samples = colnames(P2)),
        file.path(subdir2, "subject_c001.rds"))
s_na <- read_subject(subdir2, idx)
stopifnot(sum(!complete.cases(s_na$beta)) == 1L)
cat("  a no-read subject/CpG survives as NA for the caller to handle     PASS\n")

e5 <- tryCatch(read_subject(subdir, c(idx, 987654L)), error = conditionMessage)
stopifnot(grepl("different marker set", e5))
e6 <- tryCatch(read_subject(file.path(tmp, "nothing-here"), idx), error = conditionMessage)
stopifnot(grepl("no subject_c", e6))
cat("  a mismatched or missing subject directory is refused             PASS\n")

unlink(tmp, recursive = TRUE)
cat("\nDECONVOLVE TESTS PASSED\n")
