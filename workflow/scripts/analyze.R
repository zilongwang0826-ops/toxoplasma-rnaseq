#!/usr/bin/env Rscript
# DESeq2 inference from Salmon estimates; TPM is never used directly as counts.
parse_options <- function(args, defaults = list()) {
  if (!length(args)) return(defaults)
  if (length(args) %% 2L || any(!grepl("^--[a-z][a-z0-9-]*$", args[seq(1L, length(args), 2L)]))) stop("Expected --name value pairs")
  keys <- sub("^--", "", args[seq(1L, length(args), 2L)])
  if (anyDuplicated(keys)) stop("Duplicate command-line option")
  for (i in seq_along(keys)) defaults[[keys[i]]] <- args[2L * i]
  defaults
}
need_packages <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing R packages: ", paste(absent, collapse = ", "), ". Use the declared analysis container.")
}
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
strip_version <- function(x) sub("\\.[0-9]+$", "", x)
validate_samples <- function(samples) {
  required <- c("sample_id", "condition", "day")
  if (!all(required %in% names(samples))) stop("Missing sample columns")
  if (anyNA(samples[required]) || any(!nzchar(samples$sample_id)) || anyDuplicated(samples$sample_id)) stop("Missing or duplicate sample identity")
  if (any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", samples$sample_id))) stop("Unsafe sample ID")
  samples$day <- as.character(samples$day)
  if (!all(samples$condition %in% c("control", "infected")) || !all(samples$day %in% c("11", "33"))) stop("Invalid condition or time point")
  groups <- table(factor(samples$day, levels = c("11", "33")), factor(samples$condition, levels = c("control", "infected")))
  if (!all(groups == 3L) || nrow(samples) != 12L) stop("Study requires 12 libraries, 3 in each day-by-condition group")
  samples$condition <- factor(samples$condition, levels = c("control", "infected"))
  samples$day <- factor(samples$day, levels = c("11", "33"))
  rownames(samples) <- samples$sample_id
  samples
}
validate_txi <- function(txi, samples) {
  required <- c("counts", "abundance", "length", "countsFromAbundance")
  if (!all(required %in% names(txi))) stop("Incomplete tximport object")
  for (field in c("counts", "abundance", "length")) {
    m <- txi[[field]]
    if (!is.matrix(m) || any(!is.finite(m)) || any(m < 0)) stop("Invalid matrix: ", field)
    if (!identical(dimnames(m), dimnames(txi$counts))) stop("Inconsistent tximport matrix identities")
  }
  if (!identical(colnames(txi$counts), as.character(samples$sample_id))) stop("Count/sample order mismatch")
  if (is.null(rownames(txi$counts)) || anyDuplicated(rownames(txi$counts))) stop("Invalid feature IDs")
  if (any(colSums(txi$counts) <= 0)) stop("Empty library")
}
subset_txi <- function(txi, columns) {
  for (field in c("counts", "abundance", "length")) txi[[field]] <- txi[[field]][, columns, drop = FALSE]
  txi
}
new_dataset <- function(txi, samples, design) {
  validate_txi(txi, samples)
  model <- stats::model.matrix(design, data = samples)
  if (qr(model)$rank != ncol(model)) stop("Design is not full rank")
  DESeq2::DESeqDataSetFromTximport(txi, colData = samples, design = design)
}
save_plot <- function(stem, draw, width = 7, height = 6) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(paste0(stem, ".pdf"), width = width, height = height, useDingbats = FALSE)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 300)
  tryCatch(draw(), finally = grDevices::dev.off())
}
draw_expression_heatmap <- function(z, samples, title) {
  row_order <- stats::hclust(stats::dist(z))$order
  col_order <- stats::hclust(stats::dist(t(z)))$order
  shown <- z[rev(row_order), col_order, drop = FALSE]
  limit <- max(abs(shown))
  palette <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(7, 1))
  on.exit(graphics::layout(1))
  graphics::par(mar = c(7, 1, 4, 9))
  graphics::image(seq_len(ncol(shown)), seq_len(nrow(shown)), t(shown), col = palette,
                  zlim = c(-limit, limit), axes = FALSE, xlab = "", ylab = "", main = title)
  graphics::axis(4, at = seq_len(nrow(shown)), labels = rownames(shown), las = 1, tick = FALSE, cex.axis = 0.55)
  label_colours <- c(control = "#31688E", infected = "#C95D31")[as.character(samples[colnames(shown), "condition"])]
  graphics::text(seq_len(ncol(shown)), graphics::par("usr")[3] - 0.7, labels = colnames(shown), srt = 90,
                 adj = 1, xpd = NA, cex = 0.8, col = label_colours)
  graphics::mtext("Columns: control (blue), infected (orange)", side = 1, line = 5, cex = 0.65)
  graphics::par(mar = c(7, 1, 4, 3))
  scale <- seq(-limit, limit, length.out = 101)
  graphics::image(1, scale, matrix(scale, nrow = 1), col = palette, axes = FALSE, xlab = "", ylab = "", zlim = c(-limit, limit))
  graphics::axis(4, las = 1, cex.axis = 0.7)
  graphics::mtext("Row Z score", side = 4, line = 2, cex = 0.8)
}
plot_qc <- function(dds, samples, outdir, title_prefix) {
  # VST uses all count-filtered features, not a significance-selected gene list.
  transformed <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  values <- SummarizedExperiment::assay(transformed)
  variance <- apply(values, 1L, stats::var)
  if (sum(variance > 0) < 3L) stop("Too few variable features for QC")
  selected <- order(variance, decreasing = TRUE)[seq_len(min(500L, sum(variance > 0)))]
  pca <- stats::prcomp(t(values[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
  explained <- pca$sdev^2 / sum(pca$sdev^2)
  write_tsv(data.frame(sample_id = samples$sample_id, condition = samples$condition, day = samples$day,
                       PC1 = pca$x[, 1], PC2 = pca$x[, 2]), file.path(outdir, "pca-coordinates.tsv"))
  write_tsv(data.frame(component = seq_along(explained), fraction_variance = explained), file.path(outdir, "pca-variance.tsv"))
  colours <- c(control = "#31688E", infected = "#C95D31")[as.character(samples$condition)]
  shapes <- c(`11` = 16, `33` = 17)[as.character(samples$day)]
  save_plot(file.path(outdir, "pca"), function() {
    graphics::par(mar = c(4.5, 4.5, 3, 9))
    graphics::plot(pca$x[, 1:2], pch = shapes, col = colours, cex = 1.3,
                   xlab = sprintf("PC1 (%.1f%%)", explained[1] * 100), ylab = sprintf("PC2 (%.1f%%)", explained[2] * 100), main = paste(title_prefix, "PCA"))
    bounds <- graphics::par("usr")
    graphics::legend(bounds[2] + 0.05 * diff(bounds[1:2]), bounds[4],
                     paste(samples$sample_id, samples$condition, paste0(samples$day, " dpi")),
                     col = colours, pch = shapes, bty = "n", cex = 0.72, xpd = NA, xjust = 0, title = "Sample / condition / day")
  }, 9, 7)
  distance <- as.matrix(stats::dist(t(values)))
  write_tsv(data.frame(sample_id = rownames(distance), distance, check.names = FALSE), file.path(outdir, "sample-distances.tsv"))
  save_plot(file.path(outdir, "sample-distances"), function() {
    stats::heatmap(distance, Rowv = as.dendrogram(stats::hclust(stats::as.dist(distance))), Colv = "Rowv", scale = "none", symm = TRUE,
                   col = grDevices::hcl.colors(64, "Blues 3", rev = TRUE), margins = c(8, 8), main = paste(title_prefix, "sample distances"))
  }, 8, 8)
  top <- order(variance, decreasing = TRUE)[seq_len(min(40L, sum(variance > 0)))]
  z <- t(scale(t(values[top, , drop = FALSE])))
  write_tsv(data.frame(feature_id = rownames(z), z, check.names = FALSE), file.path(outdir, "heatmap-row-zscores.tsv"))
  save_plot(file.path(outdir, "variable-features-heatmap"), function() {
    draw_expression_heatmap(z, samples, paste(title_prefix, "top variable features"))
  }, 9, 10)
  invisible(values)
}
annotate_result <- function(result, ids, keep, alpha, lfc) {
  frame <- as.data.frame(result)
  frame$feature_id <- rownames(frame)
  frame <- frame[match(ids, frame$feature_id), , drop = FALSE]
  frame$feature_id <- ids
  frame$filter_pass <- keep
  frame$test_status <- ifelse(!keep, "low_count", ifelse(is.na(frame$pvalue), "cooks_or_untestable", ifelse(is.na(frame$padj), "independent_filter", "tested")))
  frame$significant <- !is.na(frame$padj) & frame$padj < alpha & abs(frame$log2FoldChange) > lfc
  frame$nominal_sensitivity <- !is.na(frame$pvalue) & frame$pvalue < alpha & abs(frame$log2FoldChange) >= lfc
  frame$direction <- ifelse(frame$significant, ifelse(frame$log2FoldChange > 0, "up", "down"), "not_significant")
  rownames(frame) <- NULL
  frame[, c("feature_id", setdiff(names(frame), "feature_id"))]
}
fit_contrast <- function(txi, samples, outdir, alpha = 0.05, lfc = 1.5, min_count = 10L, min_samples = 3L, level = "gene", make_plots = TRUE) {
  need_packages("DESeq2")
  if (length(unique(samples$day)) != 1L) stop("Single time point required")
  dds <- new_dataset(txi, samples, ~ condition)
  enough <- rowSums(DESeq2::counts(dds) >= min_count)
  keep <- enough >= min_samples
  if (sum(keep) < 10L) stop("Fewer than 10 features pass count filtering")
  write_tsv(data.frame(feature_id = rownames(dds), samples_at_min_count = enough, filter_pass = keep), file.path(outdir, "feature-filter.tsv"))
  dds <- DESeq2::DESeq(dds[keep, ], parallel = FALSE, quiet = TRUE)
  result <- DESeq2::results(dds, contrast = c("condition", "infected", "control"), alpha = alpha, pAdjustMethod = "BH", independentFiltering = TRUE)
  table <- annotate_result(result, rownames(txi$counts), keep, alpha, lfc)
  write_tsv(table, file.path(outdir, "results.tsv"))
  write_tsv(table[table$significant, , drop = FALSE], file.path(outdir, "significant.tsv"))
  normalized <- DESeq2::counts(dds, normalized = TRUE)
  write_tsv(data.frame(feature_id = rownames(normalized), normalized, check.names = FALSE), file.path(outdir, "normalized-counts.tsv"))
  factors <- DESeq2::sizeFactors(dds)
  if (!is.null(factors)) write_tsv(data.frame(sample_id = colnames(dds), size_factor = factors), file.path(outdir, "size-factors.tsv"))
  nf <- DESeq2::normalizationFactors(dds)
  if (!is.null(nf)) write_tsv(data.frame(sample_id = colnames(dds), median_normalization_factor = apply(nf, 2L, stats::median)), file.path(outdir, "normalization-summary.tsv"))
  saveRDS(dds, file.path(outdir, "deseq2.rds"))
  writeLines(capture.output(S4Vectors::metadata(result)), file.path(outdir, "result-metadata.txt"))
  summary <- data.frame(level = level, day = as.character(samples$day[1]), contrast = "infected_vs_control", n_control = sum(samples$condition == "control"),
                        n_infected = sum(samples$condition == "infected"), features_input = nrow(txi$counts), features_filter_pass = sum(keep), features_eligible = sum(table$test_status == "tested"),
                        significant_up = sum(table$direction == "up"), significant_down = sum(table$direction == "down"), nominal_sensitivity = sum(table$nominal_sensitivity),
                        alpha = alpha, lfc_cutoff = lfc, counts_from_abundance = txi$countsFromAbundance)
  write_tsv(summary, file.path(outdir, "summary.tsv"))
  if (make_plots) {
    title <- paste(tools::toTitleCase(level), as.character(samples$day[1]), "dpi")
    plot_qc(dds, samples, file.path(outdir, "figures"), title)
    save_plot(file.path(outdir, "figures", "ma"), function() {
      DESeq2::plotMA(result, alpha = alpha, main = paste(title, "MA plot"), ylim = c(-8, 8))
      graphics::abline(h = c(-lfc, lfc), lty = 3, col = "grey40")
    })
    save_plot(file.path(outdir, "figures", "volcano"), function() {
      ok <- is.finite(table$log2FoldChange) & !is.na(table$padj)
      if (!any(ok)) { graphics::plot.new(); graphics::title(main = paste(title, "no FDR-eligible features")); return(invisible(NULL)) }
      graphics::plot(table$log2FoldChange[ok], -log10(pmax(table$padj[ok], .Machine$double.xmin)), pch = 16, cex = 0.45,
                     col = c(up = "#C95D31", down = "#31688E", not_significant = "grey75")[table$direction[ok]],
                     xlab = "log2 fold change: infected / control", ylab = "-log10 BH-adjusted P", main = paste(title, "volcano"))
      graphics::abline(v = c(-lfc, lfc), h = -log10(alpha), lty = 3, col = "grey40")
      graphics::legend("top", c("Up", "Down", "Other"), col = c("#C95D31", "#31688E", "grey75"), pch = 16, bty = "n")
    })
  }
  invisible(list(dds = dds, results = table, summary = summary))
}
fit_interaction <- function(txi, samples, outdir, alpha, min_count, min_samples) {
  dds <- new_dataset(txi, samples, ~ day + condition + day:condition)
  keep <- rowSums(DESeq2::counts(dds) >= min_count) >= min_samples
  if (sum(keep) < 10L) stop("Insufficient features for interaction model")
  dds <- DESeq2::DESeq(dds[keep, ], parallel = FALSE, quiet = TRUE)
  coefficient <- grep("day33.*conditioninfected|conditioninfected.*day33", DESeq2::resultsNames(dds), value = TRUE)
  if (length(coefficient) != 1L) stop("Cannot identify unique interaction coefficient")
  result <- DESeq2::results(dds, name = coefficient, alpha = alpha)
  frame <- as.data.frame(result)
  write_tsv(data.frame(feature_id = rownames(frame), frame, row.names = NULL), file.path(outdir, "results.tsv"))
  writeLines(c("Exploratory: (infected33 - control33) - (infected11 - control11).", "Cross-sectional design; no pairing or repeated measures are assumed.", paste("Coefficient:", coefficient)), file.path(outdir, "interpretation.txt"))
  saveRDS(dds, file.path(outdir, "deseq2.rds"))
}
validate_quant_files <- function(files) {
  expected <- NULL
  for (path in files) {
    q <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!identical(names(q), c("Name", "Length", "EffectiveLength", "TPM", "NumReads")) || !nrow(q)) stop("Malformed quant.sf: ", path)
    if (anyDuplicated(q$Name) || any(!nzchar(q$Name))) stop("Invalid transcript identities: ", path)
    if (any(!is.finite(as.matrix(q[, -1]))) || any(as.matrix(q[, -1]) < 0) || sum(q$NumReads) <= 0) stop("Invalid quantification values: ", path)
    if (is.null(expected)) expected <- q$Name
    if (!identical(q$Name, expected)) stop("Transcript identities/order differ across samples")
  }
  expected
}
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (identical(args, "--help")) {
    cat("Usage: analyze.R --samples FILE --quant-root DIR --tx2gene FILE --outdir DIR [--alpha 0.05 --lfc 1.5 --min-count 10 --min-samples 3 --interaction false]\n")
    return(invisible(NULL))
  }
  opt <- parse_options(args, list(alpha = "0.05", lfc = "1.5", `min-count` = "10", `min-samples` = "3", interaction = "false"))
  required <- c("samples", "quant-root", "tx2gene", "outdir")
  if (!all(required %in% names(opt))) stop("Missing required arguments")
  if (!all(names(opt) %in% c(required, "alpha", "lfc", "min-count", "min-samples", "interaction"))) stop("Unknown option")
  alpha <- as.numeric(opt$alpha); lfc <- as.numeric(opt$lfc); min_count <- as.integer(opt[["min-count"]]); min_samples <- as.integer(opt[["min-samples"]])
  if (anyNA(c(alpha, lfc, min_count, min_samples)) || alpha <= 0 || alpha >= 1 || lfc < 0 || min_count < 1 || min_samples < 2 || min_samples > 6) stop("Invalid analysis parameters")
  if (!opt$interaction %in% c("true", "false")) stop("--interaction must be true or false")
  need_packages(c("tximport", "DESeq2"))
  samples <- validate_samples(utils::read.delim(opt$samples, stringsAsFactors = FALSE, check.names = FALSE))
  files <- file.path(opt[["quant-root"]], samples$sample_id, "quant.sf"); names(files) <- samples$sample_id
  if (any(!file.exists(files))) stop("Missing quant.sf: ", paste(files[!file.exists(files)], collapse = ", "))
  ids <- validate_quant_files(files)
  mapping <- utils::read.delim(opt$tx2gene, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("transcript_id", "gene_id") %in% names(mapping))) stop("tx2gene requires transcript_id,gene_id")
  mapping$transcript_id <- strip_version(mapping$transcript_id); mapping$gene_id <- strip_version(mapping$gene_id)
  tx2gene <- unique(mapping[, c("transcript_id", "gene_id")])
  if (anyNA(tx2gene) || anyDuplicated(tx2gene$transcript_id) || any(!nzchar(tx2gene$gene_id))) stop("Ambiguous transcript-to-gene mapping")
  if (anyDuplicated(strip_version(ids)) || !all(strip_version(ids) %in% tx2gene$transcript_id)) stop("Quantified IDs mismatch reference annotation")
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(samples, file.path(opt$outdir, "samples-used.tsv"))
  write_tsv(data.frame(parameter = names(opt), value = unlist(opt, use.names = FALSE)), file.path(opt$outdir, "parameters.tsv"))
  paths <- c(opt$samples, opt$tx2gene, unname(files))
  write_tsv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))), file.path(opt$outdir, "input-checksums.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "session-info.txt"))
  imported <- list(gene = tximport::tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "no"),
                   transcript = tximport::tximport(files, type = "salmon", txOut = TRUE, countsFromAbundance = "lengthScaledTPM"))
  summaries <- list()
  for (level in names(imported)) {
    txi <- imported[[level]]; validate_txi(txi, samples)
    global_dds <- new_dataset(txi, samples, ~ day + condition + day:condition)
    keep <- rowSums(DESeq2::counts(global_dds) >= min_count) >= min_samples
    if (sum(keep) < 10L) stop("Insufficient features for global QC")
    plot_qc(global_dds[keep, ], samples, file.path(opt$outdir, level, "qc"), paste(tools::toTitleCase(level), "all samples"))
    for (day in levels(samples$day)) {
      idx <- which(samples$day == day)
      fitted <- fit_contrast(subset_txi(txi, idx), samples[idx, , drop = FALSE], file.path(opt$outdir, level, paste0("day", day)), alpha, lfc, min_count, min_samples, level)
      summaries[[paste(level, day)]] <- fitted$summary
    }
    if (opt$interaction == "true") fit_interaction(txi, samples, file.path(opt$outdir, level, "interaction"), alpha, min_count, min_samples)
  }
  write_tsv(do.call(rbind, summaries), file.path(opt$outdir, "summary.tsv"))
  writeLines(c("Status: executed", paste("Finished UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
               "Per-day design: ~ condition; contrast infected versus control.", "Primary: BH-adjusted P < alpha and |unshrunken log2FC| > cutoff.",
               "Nominal-P sensitivity is exploratory. Execution alone does not establish reproduction."), file.path(opt$outdir, "execution-status.txt"))
}
if (sys.nframe() == 0L) main()
