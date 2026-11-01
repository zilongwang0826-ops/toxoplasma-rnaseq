rule differential_expression:
    input:
        container_runner="workflow/scripts/container_exec.py",
        qc=OUT + '/qc/quantification.json', samples=CFG['samples'],
        quant=expand(OUT + '/salmon/{sample}/quant.sf',sample=IDS),
        mapping=REF + '/prepared/tx2gene.tsv', code='workflow/scripts/analyze.R'
    output: directory(OUT + '/analysis')
    params:
        image='statistics', image_id=image_identity('statistics'), quant=OUT + '/salmon', alpha=CFG['analysis']['alpha'],
        lfc=CFG['analysis']['abs_log2fc'], minimum_count=CFG['analysis']['min_count'], n=CFG['analysis']['min_samples'],
        interaction=str(CFG['analysis']['interaction']).lower()
    threads: CFG['resources']['threads']
    resources: mem_mb=CFG['resources']['analysis_mem_mb']
    log: OUT + '/logs/differential_expression.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'Rscript --vanilla {input.code:q} --samples {input.samples:q} --quant-root {params.quant:q} '
        '--tx2gene {input.mapping:q} --outdir {output:q} --alpha {params.alpha} --lfc {params.lfc} '
        '--min-count {params.minimum_count} --min-samples {params.n} --interaction {params.interaction} > {log:q} 2>&1'

rule enrichment:
    input:
        container_runner="workflow/scripts/container_exec.py",
        analysis=OUT + '/analysis', mapping=REF + '/prepared/tx2gene.tsv',
        kegg=[CFG['enrichment']['kegg_tsv']] if CFG['enrichment']['kegg_tsv'] else [],
        code=['workflow/scripts/enrichment.R','workflow/scripts/analyze.R']
    output: directory(OUT + '/enrichment')
    params:
        image='statistics', image_id=image_identity('statistics'), alpha=CFG['analysis']['alpha'],
        kegg=('--kegg-tsv ' + shlex.quote(CFG['enrichment']['kegg_tsv'])) if CFG['enrichment']['kegg_tsv'] else ''
    threads: CFG['resources']['threads']
    resources: mem_mb=CFG['resources']['analysis_mem_mb']
    log: OUT + '/logs/enrichment.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'Rscript --vanilla workflow/scripts/enrichment.R --analysis {input.analysis:q} '
        '--tx2gene {input.mapping:q} --outdir {output:q} --alpha {params.alpha} {params.kegg} > {log:q} 2>&1'

rule compare_published:
    input:
        container_runner="workflow/scripts/container_exec.py",
        analysis=OUT + '/analysis', published='metadata/published_dets.tsv',
        code=['workflow/scripts/compare_published.R','workflow/scripts/analyze.R']
    output: directory(OUT + '/comparison')
    params: image='statistics', image_id=image_identity('statistics')
    resources: mem_mb=2000
    log: OUT + '/logs/comparison.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'Rscript --vanilla workflow/scripts/compare_published.R --analysis {input.analysis:q} '
        '--published {input.published:q} --outdir {output:q} > {log:q} 2>&1'

rule locus_sensitivity:
    input:
        container_runner="workflow/scripts/container_exec.py",
        analysis=OUT + '/analysis', enrichment=OUT + '/enrichment',
        code=['workflow/scripts/locus_sensitivity.R','workflow/scripts/analyze.R','workflow/scripts/enrichment.R']
    output: directory(OUT + '/sensitivity/lilrb4a')
    params: image='statistics', image_id=image_identity('statistics')
    threads: CFG['resources']['threads']
    resources: mem_mb=CFG['resources']['analysis_mem_mb']
    log: OUT + '/logs/locus_sensitivity.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'Rscript --vanilla workflow/scripts/locus_sensitivity.R --analysis {input.analysis:q} '
        '--enrichment {input.enrichment:q} --outdir {output:q} > {log:q} 2>&1'

rule publication_figures:
    input:
        container_runner="workflow/scripts/container_exec.py",
        analysis=OUT + '/analysis', enrichment=OUT + '/enrichment', sensitivity=OUT + '/sensitivity/lilrb4a',
        code=['workflow/scripts/publication_figures.R','workflow/scripts/analyze.R']
    output: directory(OUT + '/figures')
    params: image='statistics', image_id=image_identity('statistics')
    threads: CFG['resources']['threads']
    resources: mem_mb=2000
    log: OUT + '/logs/publication_figures.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'Rscript --vanilla workflow/scripts/publication_figures.R --analysis {input.analysis:q} '
        '--enrichment {input.enrichment:q} --sensitivity {input.sensitivity:q} --outdir {output:q} > {log:q} 2>&1'

rule model_diagnostics:
    input:
        container_runner="workflow/scripts/container_exec.py",
        analysis=OUT + '/analysis', code='workflow/scripts/model_diagnostics.R'
    output: directory(OUT + '/validation/models')
    params: image='statistics', image_id=image_identity('statistics')
    threads: CFG['resources']['threads']
    resources: mem_mb=4000
    log: OUT + '/logs/model_diagnostics.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'Rscript --vanilla {input.code:q} --analysis {input.analysis:q} --outdir {output:q} > {log:q} 2>&1'

rule provenance:
    input:
        container_runner="workflow/scripts/container_exec.py",
        analysis=OUT + '/analysis', enrichment=OUT + '/enrichment', comparison=OUT + '/comparison',
        sensitivity=OUT + '/sensitivity/lilrb4a', diagnostics=OUT + '/validation/models', figures=OUT + '/figures',
        reference=REF + '/prepared/reference.json', qc=OUT + '/qc/quantification.json',
        code=['workflow/scripts/provenance.py','workflow/scripts/common.py']
    output: OUT + '/provenance/run.json'
    params: image='tools', image_id=image_identity('tools'), effective=json.dumps(CFG,sort_keys=True)
    resources: mem_mb=1000
    log: OUT + '/logs/provenance.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'python workflow/scripts/provenance.py --output {output:q} --config-json {params.effective:q} --reference {input.reference:q} --qc {input.qc:q} > {log:q} 2>&1'

rule report:
    input:
        container_runner="workflow/scripts/container_exec.py",
        analysis=OUT + '/analysis', enrichment=OUT + '/enrichment', comparison=OUT + '/comparison',
        sensitivity=OUT + '/sensitivity/lilrb4a', diagnostics=OUT + '/validation/models', figures=OUT + '/figures',
        provenance=OUT + '/provenance/run.json', multiqc=OUT + '/qc/multiqc_report.html',
        code='workflow/scripts/build_report.py'
    output: OUT + '/report/index.html'
    params: image='tools', image_id=image_identity('tools'), root=OUT
    resources: mem_mb=1000
    log: OUT + '/logs/report.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'python {input.code:q} --results {params.root:q} --output {output:q} > {log:q} 2>&1'
