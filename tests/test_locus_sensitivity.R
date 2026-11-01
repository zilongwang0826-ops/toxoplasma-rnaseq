#!/usr/bin/env Rscript
source("workflow/scripts/locus_sensitivity.R")
mapping <- c(a = "L", b = "L")
stopifnot(identical(collapse_ids(c("a", "b", "c"), mapping), c("L", "L", "c")))
ann <- data.frame(gene_id = c("a", "b", "c", "a", "b", "c", "d", "e", "f", "g", "h"),
                  term_id = c(rep("small", 3), rep("large", 5), rep("zero", 3)),
                  term_name = c(rep("Small", 3), rep("Large", 5), rep("Zero", 3)), ontology = "BP")
universe <- letters[1:8]; selected <- c("a", "b")
base <- counting_ora(selected, universe, ann, "up", 3, 4)
old <- ora(selected, universe, ann, "up", 3, 4)
stopifnot(isTRUE(all.equal(base, old, check.attributes = FALSE)))
new_ann <- ann; new_ann$gene_id <- collapse_ids(new_ann$gene_id, mapping)
after <- counting_ora(collapse_ids(selected, mapping), collapse_ids(universe, mapping), new_ann, "up", 3, 4)
stopifnot(setequal(base$term_id, c("small", "zero")), setequal(after$term_id, c("large", "zero")),
          all(after$selected_size == 1), all(after$universe_size == 7),
          after$overlap[after$term_id == "large"] == 1,
          after$pvalue[after$term_id == "zero"] == 1,
          isTRUE(all.equal(after$padj, p.adjust(after$pvalue, "BH"))))
stopifnot(nrow(counting_ora(character(), universe, ann, "down", 3, 4)) == 0L)
cat("PASS: exact primary ORA agreement, locus deduplication, 3-to-2/5-to-4 term boundaries, zero-hit BH family and empty selection.\n")
