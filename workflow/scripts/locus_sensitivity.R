#!/usr/bin/env Rscript
# Counting sensitivity for one documented overlapping Ensembl locus, not a DE refit.
source("workflow/scripts/analyze.R")
source("workflow/scripts/enrichment.R")

collapse_ids <- function(ids, mapping) {
  hit <- match(ids, names(mapping))
  ids[!is.na(hit)] <- unname(mapping[hit[!is.na(hit)]])
  ids
}

counting_ora <- function(selected, universe, annotations, direction, min_size, max_size) {
  universe <- unique(universe)
  selected <- intersect(unique(selected), universe)
  ann <- unique(annotations[annotations$gene_id %in% universe, ])
  sets <- lapply(split(ann$gene_id, ann$term_id), unique)
  sets <- sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
  if (!length(sets) || !length(selected) || !length(universe)) return(empty_ora())
  terms <- names(sets)
  descriptions <- ann[match(terms, ann$term_id), ]
  hits <- lapply(sets, intersect, y = selected)
  result <- data.frame(term_id = terms, term_name = descriptions$term_name,
    ontology = descriptions$ontology, direction = direction, overlap = lengths(hits),
    term_size = lengths(sets), selected_size = length(selected), universe_size = length(universe),
    enrichment_ratio = (lengths(hits) / length(selected)) / (lengths(sets) / length(universe)),
    pvalue = stats::phyper(lengths(hits) - 1L, lengths(sets), length(universe) - lengths(sets),
                          length(selected), lower.tail = FALSE),
    padj = NA_real_, genes = vapply(hits, function(x) paste(sort(x), collapse = ";"), character(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  result$padj <- stats::p.adjust(result$pvalue, "BH")
  result[order(result$padj, result$pvalue, result$term_id), , drop = FALSE]
}

locus_sensitivity_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_options(args)
  if (!setequal(names(opt), c("analysis", "enrichment", "outdir"))) stop("Required: --analysis DIR --enrichment DIR --outdir NEW_DIR")
  if (dir.exists(opt$outdir) && length(list.files(opt$outdir, all.files = TRUE, no.. = TRUE))) stop("Refusing to overwrite existing sensitivity outputs")
  paths <- file.path(opt$analysis, "gene", c("day11", "day33"), "results.tsv")
  snapshot <- file.path(opt$enrichment, "go-annotation-snapshot.tsv")
  params_path <- file.path(opt$enrichment, "parameters.tsv")
  if (!all(file.exists(c(paths, snapshot, params_path, file.path(opt$enrichment, "summary.tsv"))))) stop("Completed primary analysis and enrichment are required")
  params <- utils::read.delim(params_path, colClasses = "character")
  value <- setNames(params$value, params$parameter)
  alpha <- as.numeric(value[["alpha"]]); min_size <- as.integer(value[["min-size"]]); max_size <- as.integer(value[["max-size"]])
  if (anyNA(c(alpha, min_size, max_size))) stop("Missing primary enrichment parameters")
  ann <- utils::read.delim(snapshot, colClasses = "character", check.names = FALSE)
  if (!all(c("gene_id", "term_id", "term_name", "ontology") %in% names(ann))) stop("Invalid annotation snapshot")
  pair <- c("ENSMUSG00000062593", "ENSMUSG00000112148")
  mapping <- setNames(rep("LOCUS_Lilrb4a_GRCm38_chr10", 2), pair)
  pair_terms <- lapply(pair, function(id) unique(ann$term_id[ann$gene_id == id]))
  if (any(lengths(pair_terms) == 0L)) stop("Documented pair lacks frozen GO annotation")
  equal_terms <- setequal(pair_terms[[1]], pair_terms[[2]])
  collapsed <- ann; collapsed$gene_id <- collapse_ids(collapsed$gene_id, mapping)
  collapsed <- unique(collapsed)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(data.frame(gene_id = pair, locus_id = unname(mapping), annotated_terms = lengths(pair_terms),
                      identical_term_sets = equal_terms), file.path(opt$outdir, "locus-definition.tsv"))
  write_tsv(ann[ann$gene_id %in% pair, ], file.path(opt$outdir, "locus-go-associations.tsv"))
  summaries <- list(); states <- list()
  for (i in seq_along(paths)) {
    day <- c("11", "33")[i]
    res <- utils::read.delim(paths[i], stringsAsFactors = FALSE)
    if (anyDuplicated(res$feature_id)) stop("Duplicate DE identifiers")
    states[[day]] <- cbind(day = day, res[match(pair, res$feature_id), c("feature_id", "log2FoldChange", "padj", "test_status", "significant", "direction")])
    universe <- intersect(res$feature_id[res$test_status == "tested"], unique(ann$gene_id))
    new_universe <- unique(collapse_ids(universe, mapping))
    for (direction in c("up", "down")) {
      selected <- intersect(res$feature_id[res$significant & res$direction == direction], universe)
      new_selected <- unique(collapse_ids(selected, mapping))
      before <- counting_ora(selected, universe, ann, direction, min_size, max_size)
      primary_path <- file.path(opt$enrichment, paste0("day", day), "GO", paste0(direction, ".tsv"))
      primary <- utils::read.delim(primary_path, stringsAsFactors = FALSE)
      if (!setequal(before$term_id, primary$term_id)) stop("Recomputed primary term family differs")
      primary <- primary[match(before$term_id, primary$term_id), ]
      for (column in c("overlap", "term_size", "selected_size", "universe_size", "pvalue", "padj")) {
        if (!isTRUE(all.equal(before[[column]], primary[[column]], tolerance = 1e-12, check.attributes = FALSE))) stop("Primary ORA validation failed: ", column)
      }
      after <- counting_ora(new_selected, new_universe, collapsed, direction, min_size, max_size)
      destination <- file.path(opt$outdir, paste0("day", day), direction)
      write_tsv(after, file.path(destination, "collapsed.tsv"))
      compare <- merge(before[, c("term_id", "term_name", "overlap", "term_size", "pvalue", "padj")],
                       after[, c("term_id", "overlap", "term_size", "pvalue", "padj")],
                       by = "term_id", all = TRUE, suffixes = c("_primary", "_collapsed"))
      missing_names <- is.na(compare$term_name)
      compare$term_name[missing_names] <- ann$term_name[match(compare$term_id[missing_names], ann$term_id)]
      compare$family_change <- ifelse(is.na(compare$pvalue_primary), "added", ifelse(is.na(compare$pvalue_collapsed), "removed", "retained"))
      compare$primary_significant <- !is.na(compare$padj_primary) & compare$padj_primary < alpha
      compare$collapsed_significant <- !is.na(compare$padj_collapsed) & compare$padj_collapsed < alpha
      write_tsv(compare, file.path(destination, "comparison.tsv"))
      plot_enrichment(after, file.path(destination, "enrichment"), paste(day, "dpi", direction, "GO: locus sensitivity"), alpha)
      summaries[[paste(day, direction)]] <- data.frame(day = day, direction = direction,
        primary_universe = length(universe), collapsed_universe = length(new_universe),
        primary_selected = length(selected), collapsed_selected = length(new_selected),
        primary_terms = nrow(before), collapsed_terms = nrow(after),
        terms_added = sum(compare$family_change == "added"), terms_removed = sum(compare$family_change == "removed"),
        primary_significant = sum(compare$primary_significant), collapsed_significant = sum(compare$collapsed_significant),
        significant_lost = sum(compare$primary_significant & !compare$collapsed_significant),
        significant_gained = sum(!compare$primary_significant & compare$collapsed_significant))
    }
  }
  write_tsv(do.call(rbind, states), file.path(opt$outdir, "pair-de-status.tsv"))
  write_tsv(do.call(rbind, summaries), file.path(opt$outdir, "summary.tsv"))
  inputs <- c(paths, snapshot, params_path)
  write_tsv(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))), file.path(opt$outdir, "input-checksums.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "session-info.txt"))
  writeLines(c("Completed counting sensitivity for one documented overlapping Ensembl locus.",
    "Primary DE estimates and P values are unchanged. Selected locus membership is the union of originally selected member IDs within each direction.",
    "This is not a locus-level differential-expression test; union selection does not provide new locus-level FDR control.",
    paste("Frozen GO term sets identical:", equal_terms, "; union annotations used when collapsing."),
    "Term-size eligibility and BH correction were recomputed for all four day/direction families, including zero-overlap terms.",
    "Other aliases are not automatically merged."), file.path(opt$outdir, "interpretation.txt"))
}
if (sys.nframe() == 0L) locus_sensitivity_main()
