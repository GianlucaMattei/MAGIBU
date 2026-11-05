# ==============================================================================
# Typicality p-value of a projected sample for each methylation family (Magibu 0.8.0)
# ------------------------------------------------------------------------------
# Conformal p-value that accounts for the spread of each family on the UMAP map
# (same method as Magibu/Leukemia/Scripts/31_typicality.r, see Strategia_p_tipicita.md):
#   spread of family F: for every member i, s_i = mean distance to its k nearest members of F (i excluded)
#   new sample x:       s(x) = mean distance to the k members of F closest to x
#   p(F) = (1 + #{s_i >= s(x)}) / (n_F + 1)
# p high = x is as close to F as its own members; p <= alpha = x is farther than almost all members.
# Plausible families = p > alpha (one, several = ambiguous, none = out of distribution).
# Family = text before the first comma of the class ("LGG, PA PF" -> "LGG"), so it works for any model.
# Families with fewer than min_n reference samples are not assessed (p = NA): with few members the
# spread is unreliable and small scattered families look compatible with everything.
#
# Validation (external cohorts Sturm/Hardin/CBTN, n = 1735, Capper model FinalTest.v2, k = 5, min_n = 20):
#   true family with p > 0.05 in 82% of samples; kNN5 family accuracy 70% when its p <= 0.05, 94% when p > 0.5.
#   At class level the p does NOT separate subclasses (they overlap on the map): not reported.
#   Results: Magibu/Results/Typicality/
# ==============================================================================

typicality_family <- function(x) trimws(sub(",.*", "", x))

# Reference spread: per family the sorted s_i (families with >= min_n members)
typicality_reference <- function(ref_xy, ref_group, k = 5, min_n = 20) {
  tab <- table(ref_group)
  out <- list()
  for (g in names(tab)[tab >= min_n]) {
    X <- ref_xy[ref_group == g, , drop = FALSE]
    D <- as.matrix(stats::dist(X)); diag(D) <- Inf
    kk <- min(k, nrow(X) - 1)
    out[[g]] <- sort(apply(D, 1, function(d) mean(sort(d)[seq_len(kk)])))
  }
  out
}

#' @param metadata.in projected metadata (UMAP1, UMAP2, sample.title, Batch == "New_Projected" for new samples)
#' @param sample_ids  Sample_ID_Full of the samples to assess (default: all New_Projected)
#' @return list(summary = one row per sample, all_p = long table Sample x Family)
compute_typicality <- function(metadata.in, sample_ids = NULL, group_fun = typicality_family,
                               k = 5, min_n = 20, alpha = 0.05) {
  meta <- data.table::as.data.table(metadata.in)
  is_new <- !is.na(meta$Batch) & meta$Batch == "New_Projected"
  ref <- meta[!is_new & !is.na(sample.title)]
  if (is.null(sample_ids)) sample_ids <- meta$Sample_ID_Full[is_new]
  ref_xy <- as.matrix(ref[, .(UMAP1, UMAP2)])
  grp <- group_fun(ref$sample.title)
  rs <- typicality_reference(ref_xy, grp, k = k, min_n = min_n)
  n_ref <- table(grp)

  all_p <- list(); summ <- list()
  for (s in sample_ids) {
    t <- meta[Sample_ID_Full == s][1]
    d <- sqrt((ref_xy[, 1] - t$UMAP1)^2 + (ref_xy[, 2] - t$UMAP2)^2)
    p <- vapply(names(rs), function(g) {
      dg <- sort(d[grp == g]); kk <- min(k, length(dg) - 1)
      sx <- mean(dg[seq_len(kk)]); n <- length(rs[[g]])
      (1 + n - findInterval(sx, rs[[g]], left.open = TRUE)) / (n + 1)
    }, 0)
    all_p[[s]] <- data.table::data.table(Sample = s, Family = names(p), N_ref = as.integer(n_ref[names(p)]),
                                         p = round(p, 3))[order(-p)]
    # family voted by the 5 nearest reference samples (k-NN majority, as in the Summary)
    tb <- sort(table(grp[order(d)][seq_len(min(k, length(d)))]), decreasing = TRUE)
    fam <- names(tb)[1]
    pf <- if (fam %in% names(p)) unname(p[fam]) else NA_real_
    plaus <- names(sort(p[p > alpha], decreasing = TRUE))
    flag <- if (length(plaus) == 0) "OUT OF DISTRIBUTION: no family with p > alpha"
            else if (is.na(pf)) paste0("NOT ASSESSED: family with < ", min_n, " reference samples")
            else if (pf <= alpha) "ATYPICAL: farther from its family than almost all members"
            else "OK"
    summ[[s]] <- data.table::data.table(Sample = s, Family_kNN5 = fam, Family_p = round(pf, 3), Typicality = flag,
                                        Plausible_families = paste(plaus, collapse = " | "),
                                        N_plausible = length(plaus))
  }
  list(summary = data.table::rbindlist(summ), all_p = data.table::rbindlist(all_p),
       params = list(k = k, min_n = min_n, alpha = alpha))
}
