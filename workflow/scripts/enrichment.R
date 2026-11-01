#!/usr/bin/env Rscript
# Gene-level over-representation with an explicitly tested and annotated universe.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!exists("parse_options", mode = "function")) {
  if (!length(script_arg)) stop("Source analyze.R before sourcing enrichment.R")
  source(file.path(dirname(sub("^--file=", "", script_arg[1])), "analyze.R"))
}

empty_ora <- function() data.frame(term_id = character(), term_name = character(), ontology = character(), direction = character(), overlap = integer(), term_size = integer(), selected_size = integer(), universe_size = integer(), enrichment_ratio = numeric(), pvalue = numeric(), padj = numeric(), genes = character(), stringsAsFactors = FALSE)

ora <- function(selected, universe, annotations, direction, min_size = 10L, max_size = 500L) {
  universe <- unique(universe)
  selected <- intersect(unique(selected), universe)
  annotations <- unique(annotations[annotations$gene_id %in% universe, c("term_id", "term_name", "ontology", "gene_id")])
  sets <- split(annotations$gene_id, annotations$term_id)
  sets <- lapply(sets, unique)
  sets <- sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
  if (!length(sets) || !length(selected) || !length(universe)) return(empty_ora())
  results <- lapply(names(sets), function(term) {
    genes <- sets[[term]]
    hit <- intersect(selected, genes)
    description <- annotations[match(term, annotations$term_id), ]
    data.frame(term_id = term, term_name = description$term_name, ontology = description$ontology, direction = direction,
               overlap = length(hit), term_size = length(genes), selected_size = length(selected), universe_size = length(universe),
               enrichment_ratio = (length(hit) / length(selected)) / (length(genes) / length(universe)),
               pvalue = stats::phyper(length(hit) - 1L, length(genes), length(universe) - length(genes), length(selected), lower.tail = FALSE),
               padj = NA_real_, genes = paste(sort(hit), collapse = ";"), stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, results)
  # Retain zero-overlap eligible terms so the correction family is not selected by hits.
  result$padj <- stats::p.adjust(result$pvalue, method = "BH")
  result[order(result$padj, result$pvalue, result$term_id), , drop = FALSE]
}

plot_enrichment <- function(result, path, title, alpha = 0.05) {
  selected <- result[!is.na(result$padj) & result$padj < alpha, , drop = FALSE]
  selected <- head(selected[order(selected$padj), , drop = FALSE], 15L)
  save_plot(path, function() {
    if (!nrow(selected)) {
      graphics::plot.new(); graphics::title(main = title)
      graphics::text(0.5, 0.5, "No enriched terms at the prespecified BH threshold")
    } else {
      selected <- selected[nrow(selected):1L, , drop = FALSE]
      graphics::par(mar = c(5, 17, 4, 1))
      labels <- ifelse(nchar(selected$term_name) > 65, paste0(substr(selected$term_name, 1, 62), "..."), selected$term_name)
      graphics::plot(selected$enrichment_ratio, seq_len(nrow(selected)), pch = 21, bg = "#31688E", col = "white",
                     cex = 0.7 + 2 * sqrt(selected$overlap / max(selected$overlap)), yaxt = "n", ylab = "",
                     xlab = "Observed / expected proportion", main = title)
      graphics::axis(2, at = seq_len(nrow(selected)), labels = labels, las = 1, cex.axis = 0.65)
      graphics::mtext("Point area indicates number of overlapping genes; significant terms only", side = 1, line = 3.5, cex = 0.7)
    }
  }, 11, 7)
}

fetch_go_annotation <- function(gene_ids) {
  need_packages(c("AnnotationDbi", "org.Mm.eg.db", "GO.db"))
  database <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
  available <- intersect(gene_ids, AnnotationDbi::keys(database, keytype = "ENSEMBL"))
  if (!length(available)) stop("No tested Ensembl genes map to the GO annotation database")
  table <- suppressMessages(AnnotationDbi::select(database, keys = available, keytype = "ENSEMBL", columns = c("GOALL", "ONTOLOGYALL")))
  table <- unique(table[!is.na(table$GOALL) & !is.na(table$ONTOLOGYALL), c("ENSEMBL", "GOALL", "ONTOLOGYALL")])
  names(table) <- c("gene_id", "term_id", "ontology")
  godb <- get("GO.db", envir = asNamespace("GO.db"))
  terms <- suppressMessages(AnnotationDbi::select(godb, keys = unique(table$term_id), keytype = "GOID", columns = "TERM"))
  table$term_name <- terms$TERM[match(table$term_id, terms$GOID)]
  table <- table[!is.na(table$term_name), ]
  list(annotations = table, database = database)
}

enrichment_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (identical(args, "--help")) {
    cat("Usage: enrichment.R --analysis DIR --tx2gene FILE --outdir DIR [--alpha 0.05 --min-size 10 --max-size 500 --kegg-tsv FILE]\n")
    return(invisible(NULL))
  }
  opt <- parse_options(args, list(alpha = "0.05", `min-size` = "10", `max-size` = "500", `kegg-tsv` = ""))
  required <- c("analysis", "tx2gene", "outdir")
  if (!all(required %in% names(opt)) || !all(names(opt) %in% c(required, "alpha", "min-size", "max-size", "kegg-tsv"))) stop("Missing or unknown enrichment options")
  alpha <- as.numeric(opt$alpha); min_size <- as.integer(opt[["min-size"]]); max_size <- as.integer(opt[["max-size"]])
  if (anyNA(c(alpha, min_size, max_size)) || alpha <= 0 || alpha >= 1 || min_size < 1 || max_size < min_size) stop("Invalid enrichment parameters")
  paths <- file.path(opt$analysis, "gene", c("day11", "day33"), "results.tsv")
  if (!all(file.exists(paths))) stop("Both gene-level differential-expression tables are required")
  results <- lapply(paths, utils::read.delim, stringsAsFactors = FALSE)
  for (x in results) {
    if (!all(c("feature_id", "test_status", "significant", "direction", "padj") %in% names(x)) || anyDuplicated(x$feature_id)) stop("Invalid DE result table")
  }
  ids <- unique(unlist(lapply(results, function(x) strip_version(x$feature_id[x$test_status == "tested"]))))
  go <- fetch_go_annotation(ids)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  write_tsv(go$annotations, file.path(opt$outdir, "go-annotation-snapshot.tsv"))
  write_tsv(AnnotationDbi::metadata(go$database), file.path(opt$outdir, "go-database-metadata.tsv"))
  writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "session-info.txt"))
  write_tsv(data.frame(parameter = names(opt), value = unlist(opt, use.names = FALSE)), file.path(opt$outdir, "parameters.tsv"))
  input_paths <- c(paths, opt$tx2gene)
  if (!all(file.exists(input_paths))) stop("Missing enrichment provenance input")
  families <- list(GO = go$annotations)
  if (nzchar(opt[["kegg-tsv"]])) {
    kegg <- utils::read.delim(opt[["kegg-tsv"]], stringsAsFactors = FALSE)
    if (!all(c("term_id", "term_name", "gene_id") %in% names(kegg)) || !nrow(kegg) || anyNA(kegg[, c("term_id", "term_name", "gene_id")])) stop("KEGG snapshot requires nonempty term_id,term_name,gene_id")
    kegg$gene_id <- strip_version(kegg$gene_id); kegg$ontology <- "KEGG"
    families$KEGG <- unique(kegg)
    input_paths <- c(input_paths, opt[["kegg-tsv"]])
  } else {
    writeLines("KEGG not executed: no frozen licensed gene-set snapshot was supplied. GO is the primary functional analysis.", file.path(opt$outdir, "kegg-status.txt"))
  }
  write_tsv(data.frame(path = input_paths, md5 = unname(tools::md5sum(input_paths))), file.path(opt$outdir, "input-checksums.tsv"))
  summaries <- list()
  for (i in seq_along(results)) {
    day <- c("11", "33")[i]; res <- results[[i]]; res$feature_id <- strip_version(res$feature_id)
    if (anyDuplicated(res$feature_id)) stop("Gene IDs collide after version removal")
    eligible <- res$feature_id[res$test_status == "tested"]
    for (family in names(families)) {
      annotation <- families[[family]]
      universe <- intersect(eligible, unique(annotation$gene_id))
      destination <- file.path(opt$outdir, paste0("day", day), family)
      write_tsv(data.frame(gene_id = eligible, mapped_to_family = eligible %in% universe), file.path(destination, "universe.tsv"))
      for (direction in c("up", "down")) {
        selected <- res$feature_id[res$significant & res$direction == direction]
        result <- ora(selected, universe, annotation, direction, min_size, max_size)
        write_tsv(result, file.path(destination, paste0(direction, ".tsv")))
        plot_enrichment(result, file.path(destination, paste0(direction, "-enrichment")), paste(day, "dpi", direction, family), alpha)
        summaries[[paste(day, family, direction)]] <- data.frame(day = day, family = family, direction = direction,
          eligible_genes = length(eligible), annotated_universe = length(universe), selected_genes = length(selected),
          selected_annotated = length(intersect(selected, universe)), selected_unmapped = length(setdiff(selected, universe)),
          tested_terms = nrow(result), significant_terms = sum(result$padj < alpha), alpha = alpha)
      }
    }
  }
  write_tsv(do.call(rbind, summaries), file.path(opt$outdir, "summary.tsv"))
  writeLines(c("Executed gene-level GO over-representation.", "BH correction spans all size-eligible terms across BP/CC/MF within each day and direction.",
               "GO-derived associations include ancestor propagation (GOALL); genes are counted once per term.",
               "Separate day/direction/family families do not provide study-wide FDR control.",
               "Enrichment is exploratory and does not demonstrate pathway activation or causality."), file.path(opt$outdir, "execution-status.txt"))
}
if (sys.nframe() == 0L) enrichment_main()
