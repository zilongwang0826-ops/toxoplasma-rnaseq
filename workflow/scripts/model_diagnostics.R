#!/usr/bin/env Rscript
# Inspect saved DESeq2 objects without refitting models or excluding samples.

read_tsv <- function(path) utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
finite_positive <- function(x) !is.null(x) && length(x) > 0L && all(is.finite(x) & x > 0)
numeric_summary <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(min = NA_real_, median = NA_real_, p95 = NA_real_, p99 = NA_real_, max = NA_real_))
  stats::setNames(as.numeric(stats::quantile(x, c(0, .5, .95, .99, 1), names = FALSE)), c("min", "median", "p95", "p99", "max"))
}

inspect_saved_model <- function(path, samples, level, comparison) {
  dds <- readRDS(path)
  stopifnot(inherits(dds, "DESeqDataSet"))
  result_path <- file.path(dirname(path), "results.tsv")
  result <- read_tsv(result_path)
  interaction <- identical(comparison, "interaction")
  expected <- if (interaction) samples else samples[as.character(samples$day) == sub("^day", "", comparison), , drop = FALSE]
  metadata <- as.data.frame(SummarizedExperiment::colData(dds))
  sample_ids <- as.character(expected$sample_id)
  identity <- identical(colnames(dds), sample_ids) && identical(rownames(metadata), sample_ids)
  metadata_match <- all(names(expected) %in% names(metadata)) && all(vapply(names(expected), function(k) identical(as.character(metadata[[k]]), as.character(expected[[k]])), logical(1)))
  design_text <- paste(deparse(DESeq2::design(dds)), collapse = " ")
  expected_design <- if (interaction) "~day + condition + day:condition" else "~condition"
  design_match <- identical(gsub(" ", "", design_text, fixed = TRUE), gsub(" ", "", expected_design, fixed = TRUE))
  condition_levels_match <- identical(levels(metadata$condition), c("control", "infected"))
  day_levels_match <- identical(levels(metadata$day), c("11", "33"))
  model_matrix <- stats::model.matrix(DESeq2::design(dds), metadata)
  full_rank <- qr(model_matrix)$rank == ncol(model_matrix)
  results_ids <- if (interaction) result$feature_id else result$feature_id[result$filter_pass]
  feature_order_match <- identical(rownames(dds), results_ids) && !anyDuplicated(rownames(dds))
  counts <- DESeq2::counts(dds)
  counts_valid <- all(is.finite(counts) & counts >= 0)
  normalization <- DESeq2::normalizationFactors(dds)
  lengths <- SummarizedExperiment::assays(dds)[["avgTxLength"]]
  size_factors <- DESeq2::sizeFactors(dds)
  normalization_valid <- if (level == "gene") finite_positive(normalization) && identical(dimnames(normalization), dimnames(dds)) else is.null(normalization) && finite_positive(size_factors)
  lengths_valid <- if (level == "gene") finite_positive(lengths) && identical(dimnames(lengths), dimnames(dds)) else is.null(lengths)
  mcols <- as.data.frame(S4Vectors::mcols(dds))
  beta <- mcols$betaConv
  beta_valid <- length(beta) == nrow(dds) && !anyNA(beta) && all(beta)
  dispersions <- DESeq2::dispersions(dds)
  dispersion_valid <- finite_positive(dispersions) && length(dispersions) == nrow(dds)
  cooks <- SummarizedExperiment::assays(dds)[["cooks"]]
  cooks_valid <- !is.null(cooks) && identical(dimnames(cooks), dimnames(dds)) && all(is.finite(cooks) & cooks >= 0)
  max_cooks <- mcols$maxCooks
  cutoff <- stats::qf(.99, ncol(model_matrix), nrow(model_matrix) - ncol(model_matrix))
  if (is.null(cooks) || length(max_cooks) != nrow(dds)) stop("Missing saved Cook diagnostics: ", path)
  flagged <- is.finite(max_cooks) & max_cooks > cutoff
  index <- match(rownames(dds), result$feature_id)
  if (anyNA(index)) stop("Saved features missing from result table: ", path)
  missing_p <- is.na(result$pvalue[index])
  model_id <- paste(level, comparison, sep = "/")
  checks <- c(sample_identity = identity, sample_metadata = metadata_match, design = design_match,
              condition_levels = condition_levels_match, day_levels = day_levels_match, full_rank = full_rank,
              feature_order = feature_order_match, counts = counts_valid, normalization = normalization_valid,
              lengths = lengths_valid, beta_convergence = beta_valid, dispersions = dispersion_valid, cooks = cooks_valid)
  model <- data.frame(model = model_id, design = design_text, n_samples = ncol(dds), n_features = nrow(dds),
    model_rank = qr(model_matrix)$rank, model_columns = ncol(model_matrix), all_checks_pass = all(checks),
    failed_checks = paste(names(checks)[!checks], collapse = ";"), beta_nonconverged = sum(!beta, na.rm = TRUE),
    beta_missing = sum(is.na(beta)), dispersion_nonfinite = sum(!is.finite(dispersions)),
    dispersion_nonpositive = sum(dispersions <= 0, na.rm = TRUE), dispersion_min = min(dispersions, na.rm = TRUE),
    dispersion_median = stats::median(dispersions, na.rm = TRUE), dispersion_max = max(dispersions, na.rm = TRUE),
    cooks_cutoff = cutoff, cooks_flagged = sum(flagged), missing_pvalue = sum(missing_p),
    missing_pvalue_without_saved_cook_flag = sum(missing_p & !flagged),
    saved_cook_flag_with_defined_pvalue = sum(flagged & !missing_p),
    replacement_assay_present = "replaceCounts" %in% SummarizedExperiment::assayNames(dds),
    stringsAsFactors = FALSE)
  sample_rows <- do.call(rbind, lapply(seq_len(ncol(dds)), function(j) {
    cs <- numeric_summary(cooks[, j])
    data.frame(model = model_id, sample_id = colnames(dds)[j],
      normalization_min = if (is.null(normalization)) NA_real_ else min(normalization[, j]),
      normalization_median = if (is.null(normalization)) NA_real_ else stats::median(normalization[, j]),
      normalization_max = if (is.null(normalization)) NA_real_ else max(normalization[, j]),
      avg_length_min = if (is.null(lengths)) NA_real_ else min(lengths[, j]),
      avg_length_max = if (is.null(lengths)) NA_real_ else max(lengths[, j]),
      size_factor = if (is.null(size_factors)) NA_real_ else size_factors[j],
      cooks_median = cs[["median"]], cooks_p95 = cs[["p95"]], cooks_p99 = cs[["p99"]], cooks_max = cs[["max"]],
      cooks_above_cutoff = sum(cooks[, j] > cutoff, na.rm = TRUE),
      flagged_features_above_cutoff = sum(flagged & cooks[, j] > cutoff, na.rm = TRUE), stringsAsFactors = FALSE)
  }))
  review <- flagged | missing_p | is.na(beta) | !beta
  review[is.na(review)] <- TRUE
  feature_rows <- data.frame(model = rep(model_id, sum(review)), feature_id = rownames(dds)[review],
    beta_converged = beta[review], saved_max_cooks = max_cooks[review], saved_cook_flag = flagged[review],
    pvalue_missing = missing_p[review], padj = result$padj[index][review], stringsAsFactors = FALSE)
  list(model = model, samples = sample_rows, features = feature_rows,
       inputs = c(path, result_path), checks = as.list(checks))
}

inspect_entrez_groups <- function(analysis) {
  # Mapping groups are an audit, not an instruction to merge biological entities.
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE) || !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(list(status = "not_executed: pinned annotation packages unavailable"))
  }
  results <- lapply(c("day11", "day33"), function(day) read_tsv(file.path(analysis, "gene", day, "results.tsv")))
  ids <- unique(unlist(lapply(results, function(x) x$feature_id)))
  db <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
  available <- intersect(ids, AnnotationDbi::keys(db, keytype = "ENSEMBL"))
  mapping <- suppressMessages(AnnotationDbi::select(db, keys = available, keytype = "ENSEMBL", columns = c("ENTREZID", "SYMBOL")))
  mapping <- unique(mapping[!is.na(mapping$ENTREZID), ])
  groups <- lapply(split(mapping$ENSEMBL, mapping$ENTREZID), unique)
  groups <- groups[lengths(groups) > 1L]
  group_ids <- unique(unlist(groups))
  if (!length(group_ids)) return(list(status = "executed", reference_multigroups = 0L))
  go <- suppressMessages(AnnotationDbi::select(db, keys = group_ids, keytype = "ENSEMBL", columns = "GOALL"))
  go_ids <- unique(go$ENSEMBL[!is.na(go$GOALL)])
  summaries <- list(); selected_groups <- list()
  for (i in seq_along(results)) {
    res <- results[[i]]; day <- c("11", "33")[i]
    eligible <- res$feature_id[res$test_status == "tested"]
    eligible_groups <- lapply(groups, intersect, y = eligible)
    go_eligible_groups <- lapply(eligible_groups, intersect, y = go_ids)
    for (direction in c("up", "down")) {
      selected <- res$feature_id[res$significant & res$direction == direction]
      selected_members <- lapply(groups, intersect, y = selected)
      go_selected <- lapply(selected_members, intersect, y = go_ids)
      key <- paste(day, direction)
      summaries[[key]] <- data.frame(day = day, direction = direction,
        entrez_groups_with_multiple_eligible = sum(lengths(eligible_groups) > 1L),
        go_groups_with_multiple_eligible = sum(lengths(go_eligible_groups) > 1L),
        entrez_groups_with_multiple_selected = sum(lengths(selected_members) > 1L),
        go_groups_with_multiple_selected = sum(lengths(go_selected) > 1L), stringsAsFactors = FALSE)
      for (entrez in names(groups)[lengths(selected_members) > 0L & lengths(eligible_groups) > 1L]) {
        members <- groups[[entrez]]
        selected_groups[[paste(key, entrez)]] <- data.frame(day = day, direction = direction, entrez_id = entrez,
          symbols = paste(sort(unique(mapping$SYMBOL[mapping$ENTREZID == entrez & !is.na(mapping$SYMBOL)])), collapse = ";"),
          reference_members = paste(sort(members), collapse = ";"),
          eligible_members = paste(sort(eligible_groups[[entrez]]), collapse = ";"),
          selected_members = paste(sort(selected_members[[entrez]]), collapse = ";"),
          go_eligible_n = length(go_eligible_groups[[entrez]]), go_selected_n = length(go_selected[[entrez]]), stringsAsFactors = FALSE)
      }
    }
  }
  list(status = "executed", org_package_version = as.character(utils::packageVersion("org.Mm.eg.db")),
    reference_multigroups = length(groups), reference_distinct_members = length(group_ids),
    go_reference_multigroups = sum(vapply(groups, function(g) sum(g %in% go_ids) > 1L, logical(1))),
    summary = do.call(rbind, summaries), selected_groups = if (length(selected_groups)) do.call(rbind, selected_groups) else data.frame(),
    interpretation = "Entrez mapping groups can overlap and do not establish locus equivalence; no IDs were collapsed.")
}

collect_model_diagnostics <- function(analysis) {
  if (!requireNamespace("DESeq2", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) stop("Use the declared statistics container")
  samples_path <- file.path(analysis, "samples-used.tsv")
  samples <- read_tsv(samples_path)
  if (nrow(samples) != 12L || anyDuplicated(samples$sample_id)) stop("Expected 12 distinct study libraries")
  output <- list(); inputs <- samples_path
  for (level in c("gene", "transcript")) {
    for (comparison in c("day11", "day33", "interaction")) {
      path <- file.path(analysis, level, comparison, "deseq2.rds")
      if (!file.exists(path)) {
        if (comparison == "interaction") next
        stop("Missing primary model: ", path)
      }
      key <- paste(level, comparison, sep = "/")
      output[[key]] <- inspect_saved_model(path, samples, level, comparison)
      inputs <- c(inputs, output[[key]]$inputs)
      invisible(gc())
    }
  }
  models <- do.call(rbind, lapply(output, `[[`, "model"))
  list(status = if (all(models$all_checks_pass)) "checks_passed" else "review_required",
       execution = "Read saved objects only; no model refits, count replacement, exclusions, or DE result changes.",
       models = models, samples = do.call(rbind, lapply(output, `[[`, "samples")),
       flagged_features = do.call(rbind, lapply(output, `[[`, "features")),
       checks = lapply(output, `[[`, "checks"), entrez_audit = inspect_entrez_groups(analysis),
       input_checksums = data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))),
       session_info = capture.output(sessionInfo()))
}

model_diagnostics_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (identical(args, "--help")) {
    cat("Usage: model_diagnostics.R --analysis DIR --outdir NEW_DIR\n")
    return(invisible(NULL))
  }
  if (length(args) != 4L || !setequal(args[c(1L, 3L)], c("--analysis", "--outdir"))) stop("Expected --analysis DIR --outdir NEW_DIR")
  opt <- stats::setNames(as.list(args[c(2L, 4L)]), sub("^--", "", args[c(1L, 3L)]))
  destination <- opt$outdir
  if (file.exists(destination) && (!dir.exists(destination) || length(list.files(destination, all.files = TRUE, no.. = TRUE)))) stop("Refusing to overwrite a nonempty destination: ", destination)
  result <- collect_model_diagnostics(opt$analysis)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(x, name) utils::write.table(x, file.path(destination, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_table(result$models, "models.tsv")
  write_table(result$samples, "samples.tsv")
  write_table(result$flagged_features, "flagged-features.tsv")
  write_table(result$input_checksums, "input-checksums.tsv")
  if (!is.null(result$entrez_audit$summary)) write_table(result$entrez_audit$summary, "entrez-group-summary.tsv")
  if (!is.null(result$entrez_audit$selected_groups)) write_table(result$entrez_audit$selected_groups, "selected-entrez-groups.tsv")
  writeLines(result$session_info, file.path(destination, "session-info.txt"))
  writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE, digits = 16, na = "null"), file.path(destination, "diagnostics.json"))
  if (result$status != "checks_passed") stop("Saved-model diagnostics require review; inspect diagnostics.json")
  invisible(result)
}

if (sys.nframe() == 0L) model_diagnostics_main()
