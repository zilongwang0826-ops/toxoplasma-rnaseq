#!/usr/bin/env Rscript
# Synthetic classification labels on real annotation IDs test data integration only.
# These labels are not estimates from an experiment and cannot support biology.
source("workflow/scripts/analyze.R")
source("workflow/scripts/enrichment.R")
need_packages(c("AnnotationDbi", "org.Mm.eg.db", "GO.db"))
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "results/tests/go-integration"
dir.create(root, recursive = TRUE, showWarnings = FALSE)
writeLines("SYNTHETIC SOFTWARE TEST: assigned significance labels, not experimental results.", file.path(root, "DATA-TYPE.txt"))
db <- get("org.Mm.eg.db", envir = asNamespace("org.Mm.eg.db"))
ids <- head(sort(AnnotationDbi::keys(db, keytype = "ENSEMBL")), 800L)
stopifnot(length(ids) == 800L)
mock <- data.frame(feature_id = ids, test_status = c(rep("tested", 750L), rep("independent_filter", 50L)),
                   significant = c(rep(TRUE, 100L), rep(FALSE, 700L)),
                   direction = c(rep("up", 50L), rep("down", 50L), rep("not_significant", 700L)),
                   padj = c(rep(0.01, 100L), rep(1, 650L), rep(NA_real_, 50L)))
for (day in c("11", "33")) write_tsv(mock, file.path(root, "input", "gene", paste0("day", day), "results.tsv"))
map <- file.path(root, "synthetic-tx2gene.tsv")
write_tsv(data.frame(transcript_id = paste0("synthetic_tx", seq_along(ids)), gene_id = ids), map)
out <- file.path(root, "enrichment")
enrichment_main(c("--analysis", file.path(root, "input"), "--tx2gene", map, "--outdir", out, "--min-size", "5"))
summaries <- read.delim(file.path(out, "summary.tsv"))
stopifnot(nrow(summaries) == 4L, all(summaries$eligible_genes == 750L),
          all(summaries$annotated_universe <= 750L), all(summaries$selected_genes == 50L),
          all(summaries$selected_annotated + summaries$selected_unmapped == 50L))
for (day in c("11", "33")) for (direction in c("up", "down")) {
  table <- read.delim(file.path(out, paste0("day", day), "GO", paste0(direction, ".tsv")))
  stopifnot(nrow(table) > 0L, all(table$overlap <= table$selected_size),
            isTRUE(all.equal(table$padj, p.adjust(table$pvalue, "BH"))))
}
plots <- list.files(out, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
stopifnot(length(plots) == 4L, all(file.info(plots)$size > 1000))
writeLines('{"status":"PASS","data_type":"synthetic_labels_with_real_GO_annotations","checks":["annotation_join","eligible_background","mapping_loss","BH","tables","plots"]}', file.path(root, "verified.json"))
cat("PASS: full GO annotation, background, enrichment, table and plotting integration on synthetic labels.\n")
