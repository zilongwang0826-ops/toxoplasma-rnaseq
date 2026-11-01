#!/usr/bin/env Rscript
# Render labelled QC figures from frozen plotting tables; no transformations or fits.
source("workflow/scripts/analyze.R")

place_labels <- function(points, widths, height) {
  placed <- matrix(NA_real_, nrow(points), 2)
  angles <- seq(0, 2 * pi, length.out = 17)[-17]
  offsets <- do.call(rbind, lapply(c(0.035, 0.06, 0.09, 0.13, 0.18, 0.24),
                                  function(r) cbind(cos(angles) * r, sin(angles) * r)))
  for (i in seq_len(nrow(points))) {
    candidates <- sweep(offsets, 2, points[i, ], "+")
    cost <- apply(candidates, 1, function(candidate) {
      w <- widths[i] / 2
      outside <- candidate[1] < w + 0.01 || candidate[1] > 0.99 - w || candidate[2] < height || candidate[2] > 1 - height
      prior <- if (i > 1) seq_len(i - 1) else integer()
      collisions <- sum(abs(placed[prior, 1] - candidate[1]) < (widths[prior] + widths[i]) / 2 + 0.008 &
                          abs(placed[prior, 2] - candidate[2]) < height + 0.008)
      point_overlap <- sum(abs(points[, 1] - candidate[1]) < w + 0.01 & abs(points[, 2] - candidate[2]) < height / 2 + 0.01)
      1000 * outside + 100 * collisions + 100 * point_overlap + sum((candidate - points[i, ])^2)
    })
    placed[i, ] <- candidates[which.min(cost), ]
  }
  placed
}

draw_labelled_pca <- function(coords, variance, title) {
  colours <- c(control = "#31688E", infected = "#C95D31")[coords$condition]
  shapes <- c(`11` = 16, `33` = 17)[as.character(coords$day)]
  padded <- function(x) { r <- range(x); r + c(-1, 1) * max(diff(r) * 0.18, 1) }
  graphics::par(mar = c(7, 5, 4, 1.5))
  graphics::plot(coords$PC1, coords$PC2, xlim = padded(coords$PC1), ylim = padded(coords$PC2),
    col = colours, pch = shapes, cex = 1.2, main = title,
    xlab = sprintf("PC1 (%.1f%%)", 100 * variance$fraction_variance[1]),
    ylab = sprintf("PC2 (%.1f%%)", 100 * variance$fraction_variance[2]))
  bounds <- graphics::par("usr"); spans <- c(diff(bounds[1:2]), diff(bounds[3:4]))
  points <- sweep(sweep(as.matrix(coords[, c("PC1", "PC2")]), 2, bounds[c(1, 3)], "-"), 2, spans, "/")
  labels <- place_labels(points, graphics::strwidth(coords$sample_id, cex = 0.8) / spans[1], graphics::strheight("M", cex = 0.8) / spans[2])
  labels <- sweep(sweep(labels, 2, spans, "*"), 2, bounds[c(1, 3)], "+")
  graphics::segments(coords$PC1, coords$PC2, labels[, 1], labels[, 2], col = "grey65", lwd = 0.6)
  graphics::text(labels[, 1], labels[, 2], labels = coords$sample_id, cex = 0.8, col = colours)
  groups <- unique(coords[, c("condition", "day")])
  graphics::legend("bottom", inset = c(0, -0.27), xpd = NA, bty = "n", ncol = 2,
    legend = paste(groups$condition, paste0(groups$day, " dpi")),
    col = c(control = "#31688E", infected = "#C95D31")[groups$condition],
    pch = c(`11` = 16, `33` = 17)[as.character(groups$day)], cex = 0.85)
}

draw_distance <- function(distance, title) {
  order <- stats::hclust(stats::as.dist(distance), method = "complete")$order
  shown <- distance[rev(order), order, drop = FALSE]
  palette <- grDevices::hcl.colors(100, "Blues 3", rev = TRUE)
  limit <- max(distance)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(7, 1.3))
  on.exit(graphics::layout(1))
  graphics::par(mar = c(6, 6, 4, 1))
  graphics::image(seq_len(ncol(shown)), seq_len(nrow(shown)), t(shown), zlim = c(0, limit),
    col = palette, axes = FALSE, xlab = "", ylab = "", main = title)
  graphics::axis(1, seq_len(ncol(shown)), colnames(shown), las = 2, tick = FALSE)
  graphics::axis(2, seq_len(nrow(shown)), rownames(shown), las = 2, tick = FALSE)
  graphics::mtext("Complete-linkage ordering; Euclidean distances on variance-stabilized values", 1, line = 4.5, cex = 0.65)
  graphics::par(mar = c(6, 1, 4, 4))
  scale <- seq(0, limit, length.out = 100)
  graphics::image(1, scale, matrix(scale, nrow = 1), zlim = c(0, limit), col = palette,
                  axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, las = 1, cex.axis = 0.8)
  graphics::mtext("Distance (VST)", 4, line = 2.5, cex = 0.8)
}

draw_go <- function(result, title, alpha) {
  selected <- result[!is.na(result$padj) & result$padj < alpha, , drop = FALSE]
  selected <- head(selected[order(selected$padj, selected$pvalue, selected$term_id), ], 15)
  if (!nrow(selected)) {
    graphics::plot.new(); graphics::title(main = title)
    graphics::text(0.5, 0.5, "No enriched terms at the prespecified BH threshold")
    return(invisible(NULL))
  }
  selected <- selected[nrow(selected):1, ]
  labels <- vapply(paste(selected$term_id, selected$term_name), function(x)
    paste(strwrap(x, width = 48), collapse = "\n"), character(1))
  lines <- lengths(strsplit(labels, "\n", fixed = TRUE))
  edges <- c(0, cumsum(lines + 0.8))
  positions <- (head(edges, -1) + tail(edges, -1)) / 2
  graphics::par(mar = c(6.5, 24, 4, 2))
  graphics::plot(selected$enrichment_ratio, positions, pch = 21, bg = "#31688E", col = "white",
    cex = 0.7 + 2 * sqrt(selected$overlap / max(selected$overlap)),
    ylim = range(edges), xlim = c(0, max(selected$enrichment_ratio) * 1.08),
    yaxt = "n", ylab = "", xlab = "Observed / expected gene proportion", main = title)
  graphics::abline(v = 1, lty = 3, col = "grey60")
  graphics::axis(2, at = positions, labels = labels, las = 1, cex.axis = 0.72)
  graphics::mtext(sprintf("BH-adjusted P < %s; up to 15 terms. Point size increases with overlapping gene count.", alpha),
                  side = 1, line = 5, cex = 0.7)
}

publication_figures_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_options(args)
  if (!setequal(names(opt), c("analysis", "enrichment", "sensitivity", "outdir"))) stop("Required: --analysis DIR --enrichment DIR --sensitivity DIR --outdir NEW_DIR")
  if (dir.exists(opt$outdir) && length(list.files(opt$outdir, all.files = TRUE, no.. = TRUE))) stop("Refusing to overwrite existing figure outputs")
  inputs <- character()
  parameter_path <- file.path(opt$enrichment, "parameters.tsv")
  parameters <- utils::read.delim(parameter_path, colClasses = "character")
  alpha <- as.numeric(parameters$value[parameters$parameter == "alpha"])
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("Invalid enrichment alpha")
  inputs <- c(inputs, parameter_path)
  for (level in c("gene", "transcript")) for (scope in c("qc", "day11/figures", "day33/figures")) {
    base <- file.path(opt$analysis, level, scope)
    files <- file.path(base, c("pca-coordinates.tsv", "pca-variance.tsv", "sample-distances.tsv"))
    inputs <- c(inputs, files)
    coords <- utils::read.delim(files[1]); variance <- utils::read.delim(files[2])
    frame <- utils::read.delim(files[3], check.names = FALSE)
    distance <- as.matrix(frame[, -1]); rownames(distance) <- frame$sample_id
    if (anyDuplicated(coords$sample_id) || !identical(coords$sample_id, rownames(distance)) ||
        !identical(rownames(distance), colnames(distance)) || any(!is.finite(distance)) ||
        !isTRUE(all.equal(distance, t(distance))) || any(diag(distance) != 0)) stop("Invalid QC table identities/distances")
    time <- if (scope == "qc") "all samples" else if (grepl("day11", scope)) "11 dpi" else "33 dpi"
    title <- paste(tools::toTitleCase(level), "expression:", time)
    dest <- file.path(opt$outdir, level, scope)
    save_plot(file.path(dest, "pca"), function() draw_labelled_pca(coords, variance, paste(title, "PCA")), 9, 7)
    save_plot(file.path(dest, "sample-distances"), function() draw_distance(distance, paste(title, "distances")), 9, 8)
  }
  for (day in c("11", "33")) for (direction in c("up", "down")) {
    primary <- file.path(opt$enrichment, paste0("day", day), "GO", paste0(direction, ".tsv"))
    sensitivity <- file.path(opt$sensitivity, paste0("day", day), direction, "collapsed.tsv")
    for (kind in c("primary", "sensitivity")) {
      path <- if (kind == "primary") primary else sensitivity
      result <- utils::read.delim(path)
      inputs <- c(inputs, path)
      stem <- if (kind == "primary") file.path(opt$outdir, "enrichment", paste0("day", day), "GO", paste0(direction, "-enrichment")) else
        file.path(opt$outdir, "sensitivity", "lilrb4a", paste0("day", day), direction, "enrichment")
      title <- paste(day, "dpi", direction, "GO", if (kind == "sensitivity") "- locus sensitivity" else "")
      save_plot(stem, function() draw_go(result, title, alpha), 13, 10)
    }
  }
  write_tsv(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))), file.path(opt$outdir, "input-checksums.tsv"))
  writeLines("Presentation-only rendering from saved PCA/distance and GO tables. No count transformation, model fit, sample exclusion or change to coordinates, distances or enrichment statistics.", file.path(opt$outdir, "interpretation.txt"))
}
if (sys.nframe() == 0L) publication_figures_main()
