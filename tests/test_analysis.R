#!/usr/bin/env Rscript
# All data here are synthetic; this test cannot establish biological reproduction.
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else "results/validation/r-analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source("workflow/scripts/analyze.R")
source("workflow/scripts/enrichment.R")
source("workflow/scripts/compare_published.R")
need_packages("DESeq2")
expect_error <- function(expression) {
  failed <- tryCatch({force(expression); FALSE}, error = function(e) TRUE)
  stopifnot(failed)
}
set.seed(570903)
samples <- expand.grid(replicate = 1:3, condition = c("control", "infected"), day = c("11", "33"), stringsAsFactors = FALSE)
samples$sample_id <- paste0("S", seq_len(nrow(samples)))
samples <- validate_samples(samples)
expect_error(validate_samples(samples[-1, ]))
expect_error(validate_samples(transform(samples, sample_id = rep("same", nrow(samples)))))
ngenes <- 1200L
mu <- exp(stats::runif(ngenes, log(80), log(1200)))
base <- matrix(rep(mu, nrow(samples)), nrow = ngenes)
base[1:70, samples$condition == "infected"] <- base[1:70, samples$condition == "infected"] * 8
base[71:140, samples$condition == "infected"] <- base[71:140, samples$condition == "infected"] / 8
base[141:200, samples$condition == "infected" & samples$day == "33"] <- base[141:200, samples$condition == "infected" & samples$day == "33"] * 10
counts <- matrix(stats::rnbinom(length(base), mu = as.vector(base), size = 35), nrow = ngenes,
                 dimnames = list(paste0("synthetic_", seq_len(ngenes)), samples$sample_id))
counts[ngenes, ] <- 0L
lengths <- matrix(stats::runif(length(base), 700, 1500), nrow = ngenes, dimnames = dimnames(counts))
abundance <- sweep(counts / lengths, 2, colSums(counts / lengths), "/") * 1e6
txi <- list(counts = counts, abundance = abundance, length = lengths, countsFromAbundance = "lengthScaledTPM")
validate_txi(txi, samples)
expect_error(validate_txi(txi, samples[rev(seq_len(nrow(samples))), ]))
bad <- txi; bad$counts[1, 1] <- -1
expect_error(validate_txi(bad, samples))
idx <- which(samples$day == "11")
fitted <- fit_contrast(subset_txi(txi, idx), samples[idx, ], file.path(outdir, "transcript", "day11"), make_plots = TRUE, level = "synthetic transcript")
r <- fitted$results
stopifnot(mean(r$significant[1:70] & r$log2FoldChange[1:70] > 0) > 0.8,
          mean(r$significant[71:140] & r$log2FoldChange[71:140] < 0) > 0.8,
          r$test_status[ngenes] == "low_count", !r$significant[ngenes],
          all(r$padj[r$test_status == "tested"] >= 0 & r$padj[r$test_status == "tested"] <= 1))
# Exercise the exact supported gene-level tximport offset path without installing tximport.
gene_txi <- txi; gene_txi$countsFromAbundance <- "no"
gene_result <- fit_contrast(subset_txi(gene_txi, idx), samples[idx, ], file.path(outdir, "gene", "day11"), make_plots = FALSE, level = "synthetic gene")
stopifnot(!is.null(DESeq2::normalizationFactors(gene_result$dds)),
          all(is.finite(DESeq2::normalizationFactors(gene_result$dds))))
fit_interaction(txi, samples, file.path(outdir, "interaction"), 0.05, 10L, 3L)
interaction <- read.delim(file.path(outdir, "interaction", "results.tsv"))
stopifnot(mean(interaction$log2FoldChange[match(paste0("synthetic_", 141:200), interaction$feature_id)] > 0) > 0.9)
plot_qc(new_dataset(txi, samples, ~ day + condition + day:condition)[1:1199, ], samples, file.path(outdir, "global-qc"), "Synthetic only")
annotation <- data.frame(term_id = rep(c("termA", "termB", "termC"), each = 20), term_name = rep(c("A", "B", "C"), each = 20),
                         ontology = "BP", gene_id = paste0("g", c(1:20, 21:40, 41:60)))
ora_result <- ora(paste0("g", 1:10), paste0("g", 1:100), annotation, "up", 10, 50)
stopifnot(nrow(ora_result) == 3L, ora_result$overlap[ora_result$term_id == "termA"] == 10,
          ora_result$pvalue[ora_result$term_id == "termB"] == 1,
          isTRUE(all.equal(ora_result$padj, p.adjust(ora_result$pvalue, "BH"))),
          nrow(ora(character(), paste0("g", 1:100), annotation, "down")) == 0L)
expected <- phyper(9, 20, 80, 10, lower.tail = FALSE)
stopifnot(isTRUE(all.equal(ora_result$pvalue[ora_result$term_id == "termA"], expected)))
plot_enrichment(ora_result, file.path(outdir, "synthetic-enrichment"), "Synthetic enrichment only")
plot_enrichment(empty_ora(), file.path(outdir, "empty-enrichment"), "Empty-result behavior")
pdf_files <- list.files(outdir, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
stopifnot(length(pdf_files) >= 10L, all(file.info(pdf_files)$size > 1000))
fake_published <- data.frame(feature_id = c("a", "b", "c"), log2FoldChange = c(2, -2, 3))
fake_reanalysis <- data.frame(feature_id = c("a", "b", "d"), log2FoldChange = c(2.2, 1.1, -4),
  pvalue = c(0.001, 0.1, 0.001), padj = c(0.003, 0.1, 0.003), significant = c(TRUE, FALSE, TRUE))
metrics <- comparison_metrics(fake_reanalysis, fake_published, "11", "significant")
stopifnot(metrics$shared_selected == 1L, metrics$published_available == 2L,
          metrics$jaccard == 0.25, metrics$same_direction_fraction == 0.5)
write_tsv(data.frame(test = c("positive_control_up", "positive_control_down", "low_count_guard", "gene_offsets", "interaction_direction", "hypergeometric", "empty_enrichment", "published_comparison"),
                     status = "PASS"), file.path(outdir, "test-summary.tsv"))
# Exercise full comparison file I/O using only synthetic identifiers and values.
comparison_fixture <- file.path(outdir, "comparison-fixture")
fake_reanalysis$nominal_sensitivity <- fake_reanalysis$significant
fake_reanalysis$direction <- c("up", "not_significant", "down")
pub_fixture <- do.call(rbind, lapply(c("11", "33"), function(day) {
  data.frame(level = "transcript", feature_id = fake_published$feature_id, day = day,
    log2FoldChange = fake_published$log2FoldChange, pvalue = NA_real_, padj = 0.001, direction = c("up", "down", "up"))
}))
for (day in c("11", "33")) write_tsv(fake_reanalysis, file.path(comparison_fixture, "transcript", paste0("day", day), "results.tsv"))
pub_path <- file.path(comparison_fixture, "published.tsv")
write_tsv(pub_fixture, pub_path)
comparison_main(c("--analysis", comparison_fixture, "--published", pub_path, "--outdir", file.path(outdir, "comparison")))
comparison_summary <- read.delim(file.path(outdir, "comparison", "summary.tsv"))
stopifnot(nrow(comparison_summary) == 4L, all(comparison_summary$shared_selected == 1L))
# Container tests additionally exercise actual Salmon-file parsing when tximport exists.
integration <- data.frame(component = c("tximport_file_integration", "GO_database_integration"), status = "NOT_RUN_missing_package")
if (requireNamespace("tximport", quietly = TRUE)) {
  fixture <- file.path(outdir, "tximport-fixture")
  for (i in seq_len(nrow(samples))) {
    q <- data.frame(Name = rownames(counts), Length = 1800, EffectiveLength = lengths[, i], TPM = abundance[, i], NumReads = counts[, i])
    write_tsv(q, file.path(fixture, "quant", samples$sample_id[i], "quant.sf"))
  }
  write_tsv(samples, file.path(fixture, "samples.tsv"))
  write_tsv(data.frame(transcript_id = rownames(counts), gene_id = paste0("synthetic_gene_", ceiling(seq_len(ngenes) / 2))), file.path(fixture, "tx2gene.tsv"))
  main(c("--samples", file.path(fixture, "samples.tsv"), "--quant-root", file.path(fixture, "quant"),
         "--tx2gene", file.path(fixture, "tx2gene.tsv"), "--outdir", file.path(outdir, "tximport-integration"), "--interaction", "true"))
  integrated <- read.delim(file.path(outdir, "tximport-integration", "summary.tsv"))
  stopifnot(nrow(integrated) == 4L, setequal(integrated$level, c("gene", "transcript")))
  integration$status[1] <- "PASS_synthetic"
}
if (all(vapply(c("AnnotationDbi", "org.Mm.eg.db", "GO.db"), requireNamespace, logical(1), quietly = TRUE))) {
  archived <- read.delim("metadata/published_dets.tsv")
  go_test <- fetch_go_annotation(head(unique(archived$gene_id), 50L))
  stopifnot(nrow(go_test$annotations) > 0L, all(go_test$annotations$ontology %in% c("BP", "CC", "MF")))
  integration$status[2] <- "PASS_database_lookup_only"
}
write_tsv(integration, file.path(outdir, "integration-status.tsv"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))
cat("PASS: synthetic DESeq2, gene offsets, directionality, interaction, plotting, input failures, and enrichment arithmetic.\n")
