# Enrichment identifier audit and executed sensitivity

The prepared Ensembl 100 reference contains 53,028 gene IDs: 52,249 represented on the primary decoy regions and 779 on alternative regions. The latter account for 2,040 transcripts. In the pinned `org.Mm.eg.db 3.20.0`, none of those 779 alternative-region gene IDs have Entrez or GO mappings. They therefore do not directly duplicate the current GO background, but their annotation loss must be reported if they become eligible or significant.

The published lists contain three alternative-region transcripts at 11 dpi and six at 33 dpi (six distinct IDs across days), including reference records labelled Tgtp1, Ifi30 and Ctla4. Retaining these targets is necessary for transcript-level comparison. Their absence from the newer annotation mapping does not imply absence of biological function.

A separate mapping ambiguity exists on the primary assembly. `ENSMUSG00000062593` (GTF symbol Gm49339; chr10:51,480,612–51,496,613, positive strand) and `ENSMUSG00000112148` (Lilrb4a; chr10:51,490,956–51,496,611, positive strand) overlap. Both map to Entrez 14728 and share 467 unique GO associations in the pinned database. Distinct Ensembl identifiers are currently counted separately; if both enter the actual eligible or significant sets, the same locus can contribute twice.

Across the complete reference, 51 GO-annotated many-to-one Entrez groups involve 94 Ensembl IDs. This mapping alone does not establish biological equivalence for every group. Do not automatically merge all such IDs or change the primary transcript comparison.

## Executed real-result assessment

Both overlapping Lilrb4a-related IDs are FDR-eligible at both days. At 11 dpi, Gm49339 has log2FC 5.0205 and adjusted P 0.06111; Lilrb4a has log2FC 3.7479 and adjusted P 2.974 × 10⁻⁷. At 33 dpi both are significantly upregulated, with log2FC 5.8960/5.4132 and adjusted P 2.280 × 10⁻¹⁹/1.730 × 10⁻³⁶. The actual frozen GO snapshot confirms identical sets of 467 terms.

The executed sensitivity maps this pair to one locus token, takes the union of original selected members and annotations, and recalculates all four families. It does not refit DE or combine P values. The actual eligible GO background loses one counting unit on each day; the selected upward list loses one unit only at day 33. Other many-to-one mappings are not collapsed: 20 GO mapping groups per day have multiple eligible members, but only this Lilrb4a group contains multiple significant members (33 dpi upward).

| Family | Tested terms, primary → collapsed | Significant terms, primary → collapsed | Significant lost / gained |
| --- | --- | --- | --- |
| 11 dpi up | 8,035 → 8,032 | 929 → 929 | 1 / 1 |
| 11 dpi down | 8,035 → 8,032 | 0 → 0 | 0 / 0 |
| 33 dpi up | 8,112 → 8,099 | 1,129 → 1,116 | 13 / 0 |
| 33 dpi down | 8,112 → 8,099 | 0 → 0 | 0 / 0 |

Term-size boundaries matter. Eleven of the 13 lost significant day-33 terms move from 10 eligible units to nine and leave the testing family; they were not simply assigned larger adjusted P values. The other two remain tested: blood coagulation changes from adjusted P 0.0460 to 0.0847, and negative regulation of MAPK cascade from 0.0326 to 0.0610. At day 11, one formerly significant term leaves the minimum-size family and one retained term crosses into significance.

The prespecified interferon response, antigen presentation, cytokine-production, T-cell-regulation and lymphocyte-proliferation themes remain strongly significant, as does adaptive immune response. Broad interpretation is therefore stable to this specific counting sensitivity; this does not establish robustness to every annotation/background choice. Individual small-term conclusions remain sensitive. The primary Ensembl-unit analysis is retained, with separate full sensitivity/comparison tables under `results/sensitivity/lilrb4a/`.

Independent review recomputed BH and hypergeometric probabilities from the sensitivity tables and found agreement within floating-point tolerance. The 69 primary analysis/enrichment/comparison TSV files retained identical SHA-256 hashes after this additional analysis.

Inputs: the comprehensive Ensembl100 GTF, prepared tx2gene mapping, selected FASTA headers, AnnotationDbi 1.68.0, org.Mm.eg.db 3.20.0 and GO.db 3.20.0. The annotation package records Ensembl source date 2024-May14 and Entrez/GO-Entrez source date 2024-Sep20, newer than the sequence annotation. The audit used read-only queries in the validated statistics container.
