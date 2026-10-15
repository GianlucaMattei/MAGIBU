#!/usr/bin/env Rscript
# Shrinks a Magibu model: keeps only the parts of prcomp that predict() needs (rotation[,1:k], center, scale).
# Usage: Rscript slim_model.R <components.rds> <out.rds>
args <- commandArgs(trailingOnly = TRUE)
m <- readRDS(args[1])
k <- if (!is.null(m$pca_components)) m$pca_components else 50L
p <- m$pca_model
slim <- list(sdev = p$sdev[seq_len(k)], rotation = p$rotation[, seq_len(k), drop = FALSE],
             center = p$center, scale = p$scale)
class(slim) <- "prcomp"
m$pca_model <- slim
saveRDS(m, args[2], compress = "xz")
cat("slim model:", round(file.size(args[2]) / 1e6, 1), "MB (from", round(file.size(args[1]) / 1e6, 1), "MB)\n")
