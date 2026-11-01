#!/usr/bin/env Rscript
# Compare executed transcript results with the published significant-only lists.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!exists("parse_options", mode = "function")) {
  if (!length(script_arg)) stop("Source analyze.R first")
  source(file.path(dirname(sub("^--file=", "", script_arg[1])), "analyze.R"))
}
comparison_metrics <- function(reproduced, published, day, selection) {
  tested <- reproduced$feature_id[!is.na(reproduced$pvalue)]
  eligible <- reproduced$feature_id[!is.na(reproduced$padj)]
  selected <- reproduced$feature_id[reproduced[[selection]]]
  known <- intersect(published$feature_id, reproduced$feature_id)
  overlap <- intersect(selected, published$feature_id)
  merged <- merge(published[, c("feature_id", "log2FoldChange")], reproduced[, c("feature_id", "log2FoldChange")], by = "feature_id", suffixes = c("_published", "_reanalysis"))
  finite <- is.finite(merged$log2FoldChange_published) & is.finite(merged$log2FoldChange_reanalysis)
  merged <- merged[finite, ]
  same_sign <- if (nrow(merged)) mean(sign(merged$log2FoldChange_published) == sign(merged$log2FoldChange_reanalysis)) else NA_real_
  correlation <- if (nrow(merged) >= 3L && stats::sd(merged$log2FoldChange_published) > 0 && stats::sd(merged$log2FoldChange_reanalysis) > 0) stats::cor(merged$log2FoldChange_published, merged$log2FoldChange_reanalysis, method = "spearman") else NA_real_
  data.frame(day = day, level = "transcript", selection = selection, published_selected = nrow(published),
             published_available = length(known), published_testable = sum(published$feature_id %in% tested), published_FDR_eligible = sum(published$feature_id %in% eligible),
             reanalysis_selected = length(selected), shared_selected = length(overlap),
             published_list_recovery = length(overlap) / nrow(published),
             jaccard = if (length(union(selected, published$feature_id))) length(overlap) / length(union(selected, published$feature_id)) else NA_real_,
             finite_effect_pairs = nrow(merged), same_direction_fraction = same_sign, effect_spearman = correlation)
}
comparison_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (identical(args, "--help")) {
    cat("Usage: compare_published.R --analysis DIR --published metadata/published_dets.tsv --outdir DIR\n")
    return(invisible(NULL))
  }
  opt <- parse_options(args)
  if (!setequal(names(opt), c("analysis", "published", "outdir"))) stop("Required options: --analysis --published --outdir")
  published <- utils::read.delim(opt$published, stringsAsFactors = FALSE)
  required <- c("level", "feature_id", "day", "log2FoldChange", "pvalue", "padj", "direction")
  if (!all(required %in% names(published)) || !nrow(published)) stop("Invalid published comparison table")
  if (any(published$level != "transcript")) stop("Published primary comparison must use transcript identifiers")
  published$feature_id <- strip_version(published$feature_id)
  if (anyDuplicated(paste(published$feature_id, published$day))) stop("Duplicate published transcript-by-day identity")
  if (!setequal(unique(as.character(published$day)), c("11", "33"))) stop("Published table requires both time points")
  summaries <- list(); paths <- opt$published
  for (day in c("11", "33")) {
    path <- file.path(opt$analysis, "transcript", paste0("day", day), "results.tsv")
    if (!file.exists(path)) stop("Executed transcript analysis is required: ", path)
    paths <- c(paths, path)
    result <- utils::read.delim(path, stringsAsFactors = FALSE)
    if (!all(c("feature_id", "log2FoldChange", "pvalue", "padj", "significant", "nominal_sensitivity", "direction") %in% names(result))) stop("Invalid executed result table")
    result$feature_id <- strip_version(result$feature_id)
    if (anyDuplicated(result$feature_id)) stop("Transcript IDs collide after version removal")
    ref <- published[as.character(published$day) == day, ]
    for (selection in c("significant", "nominal_sensitivity")) summaries[[paste(day, selection)]] <- comparison_metrics(result, ref, day, selection)
    merged <- merge(ref, result, by = "feature_id", all.x = TRUE, suffixes = c("_published", "_reanalysis"))
    merged$available_in_reanalysis <- merged$feature_id %in% result$feature_id
    merged$same_direction <- with(merged, ifelse(is.finite(log2FoldChange_reanalysis), sign(log2FoldChange_published) == sign(log2FoldChange_reanalysis), NA))
    write_tsv(merged, file.path(opt$outdir, paste0("day", day, "-published-transcript-comparison.tsv")))
    save_plot(file.path(opt$outdir, paste0("day", day, "-effect-comparison")), function() {
      ok <- is.finite(merged$log2FoldChange_published) & is.finite(merged$log2FoldChange_reanalysis)
      if (!any(ok)) {
        graphics::plot.new(); graphics::title(main = paste(day, "dpi: no comparable finite effect estimates"))
      } else {
        graphics::plot(merged$log2FoldChange_published[ok], merged$log2FoldChange_reanalysis[ok], pch = 16, cex = 0.5,
                       col = grDevices::adjustcolor("#31688E", alpha.f = 0.5), xlab = "Published transcript log2 fold change",
                       ylab = "Reanalysis transcript log2 fold change", main = paste(day, "dpi: published selected transcripts"))
        graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
        graphics::abline(h = 0, v = 0, col = "grey80")
      }
    })
  }
  write_tsv(do.call(rbind, summaries), file.path(opt$outdir, "summary.tsv"))
  write_tsv(data.frame(path = paths, md5 = unname(tools::md5sum(paths))), file.path(opt$outdir, "input-checksums.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "session-info.txt"))
  writeLines(c("Comparison uses transcript IDs, with terminal numeric Ensembl versions removed.",
               "Published tables contain significant hits only. Missing entries are not known negatives.",
               "List recovery and Jaccard are descriptive agreement statistics, not sensitivity, specificity or scientific validation.",
               "Effect correlation is conditional on the published selected set and is subject to selection bias.",
               "Gene-level DE counts must not be compared directly with the published transcript totals.",
               "The nominal-P sensitivity reuses the published adjusted-P-selected list because nominal P values were not supplied."), file.path(opt$outdir, "interpretation.txt"))
}
if (sys.nframe() == 0L) comparison_main()
