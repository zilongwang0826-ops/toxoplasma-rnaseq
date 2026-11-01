rule fastqc:
    input: container_runner='workflow/scripts/container_exec.py', reads=RAW + '/{sample}_{mate}.fastq.gz'
    output:
        html=OUT + '/qc/raw/{sample}_{mate}_fastqc.html', zip=OUT + '/qc/raw/{sample}_{mate}_fastqc.zip'
    params: image='tools', image_id=image_identity('tools'), outdir=OUT + '/qc/raw'
    threads: 2
    resources: mem_mb=2000
    log: OUT + '/logs/fastqc/{sample}_{mate}.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'fastqc --threads {threads} --outdir {params.outdir:q} {input.reads:q} > {log:q} 2>&1'

if CFG['quantification']['trim']:
    rule trim:
        input:
            container_runner="workflow/scripts/container_exec.py",
            r1=RAW + '/{sample}_1.fastq.gz', r2=RAW + '/{sample}_2.fastq.gz',
            qc=expand(OUT + '/qc/raw/{{sample}}_{mate}_fastqc.zip',mate=[1,2])
        output:
            r1=OUT + '/trimmed/{sample}_1.fastq.gz', r2=OUT + '/trimmed/{sample}_2.fastq.gz',
            json=OUT + '/qc/trim/{sample}.json', html=OUT + '/qc/trim/{sample}.html'
        params:
            image='tools', image_id=image_identity('tools'), quality=CFG['quantification']['trim_quality'], length=CFG['quantification']['trim_min_length']
        threads: CFG['resources']['threads']
        resources: mem_mb=2000
        log: OUT + '/logs/trim/{sample}.log'
        shell:
            '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
            'fastp --in1 {input.r1:q} --in2 {input.r2:q} --out1 {output.r1:q} --out2 {output.r2:q} '
            '--detect_adapter_for_pe --qualified_quality_phred {params.quality} --length_required {params.length} '
            '--thread {threads} --json {output.json:q} --html {output.html:q} > {log:q} 2>&1'

    rule fastqc_trimmed:
        input: reads=OUT + '/trimmed/{sample}_{mate}.fastq.gz', container_runner='workflow/scripts/container_exec.py'
        output:
            html=OUT + '/qc/trimmed/{sample}_{mate}_fastqc.html', zip=OUT + '/qc/trimmed/{sample}_{mate}_fastqc.zip'
        params: image='tools', image_id=image_identity('tools'), outdir=OUT + '/qc/trimmed'
        threads: 2
        resources: mem_mb=2000
        log: OUT + '/logs/fastqc_trimmed/{sample}_{mate}.log'
        shell:
            '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
            'fastqc --threads {threads} --outdir {params.outdir:q} {input.reads:q} > {log:q} 2>&1'

rule salmon_index:
    input:
        container_runner="workflow/scripts/container_exec.py",
        fasta=REF + '/prepared/gentrome.fa', decoys=REF + '/prepared/decoys.txt', provenance=REF + '/prepared/reference.json'
    output: directory(REF + '/salmon_index')
    params:
        image='tools', image_id=image_identity('tools'), k=CFG['reference']['kmer'],
        decoys=('-d ' + shlex.quote(REF + '/prepared/decoys.txt')) if CFG['reference']['decoys'] else ''
    threads: CFG['resources']['threads']
    resources: mem_mb=CFG['resources']['index_mem_mb']
    log: OUT + '/logs/reference/salmon_index.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'salmon index --transcripts {input.fasta:q} --index {output:q} --threads {threads} '
        '-k {params.k} --sparse --keepDuplicates {params.decoys} > {log:q} 2>&1'

rule salmon_quant:
    input:
        container_runner="workflow/scripts/container_exec.py",
        index=REF + '/salmon_index', r1=lambda w: fastq(w,1), r2=lambda w: fastq(w,2),
        qc=expand(OUT + '/qc/raw/{{sample}}_{mate}_fastqc.zip',mate=[1,2])
    output:
        quant=OUT + '/salmon/{sample}/quant.sf', meta=OUT + '/salmon/{sample}/aux_info/meta_info.json',
        libformat=OUT + '/salmon/{sample}/lib_format_counts.json'
    params:
        image='tools', image_id=image_identity('tools'), outdir=lambda w: OUT + '/salmon/' + w.sample, libtype=CFG['quantification']['library_type'],
        bias=(' --seqBias' if CFG['quantification']['seq_bias'] else '') + (' --gcBias' if CFG['quantification']['gc_bias'] else '')
    threads: CFG['resources']['threads']
    resources: mem_mb=CFG['resources']['quant_mem_mb']
    log: OUT + '/logs/salmon/{sample}.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --threads {threads} --memory-mb {resources.mem_mb} -- '
        'salmon quant --index {input.index:q} --libType {params.libtype} --mates1 {input.r1:q} '
        '--mates2 {input.r2:q} --output {params.outdir:q} --threads {threads} --validateMappings '
        '{params.bias} > {log:q} 2>&1'

rule quantification_qc:
    input:
        container_runner="workflow/scripts/container_exec.py",
        quant=expand(OUT + '/salmon/{sample}/quant.sf',sample=IDS),
        meta=expand(OUT + '/salmon/{sample}/aux_info/meta_info.json',sample=IDS), samples=CFG['samples'],
        code=['workflow/scripts/quant_qc.py','workflow/scripts/common.py']
    output: json=OUT + '/qc/quantification.json', table=OUT + '/qc/quantification.json.tsv'
    params:
        image='tools', image_id=image_identity('tools'), archive_gate='' if CFG['quantification']['trim'] else '--require-archive-counts', root=OUT + '/salmon', pct=CFG['quantification']['min_mapping_percent'], depth=CFG['quantification']['min_assigned_fragments']
    resources: mem_mb=1000
    log: OUT + '/logs/quantification_qc.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'python workflow/scripts/quant_qc.py --samples {input.samples:q} --quant-root {params.root:q} '
        '--output {output.json:q} --min-percent {params.pct} --min-fragments {params.depth} {params.archive_gate} > {log:q} 2>&1'

rule multiqc:
    input:
        container_runner="workflow/scripts/container_exec.py",
        qc=expand(OUT + '/qc/raw/{sample}_{mate}_fastqc.zip',sample=IDS,mate=[1,2]),
        quant=expand(OUT + '/salmon/{sample}/aux_info/meta_info.json',sample=IDS),
        trimmed=expand(OUT + '/qc/trimmed/{sample}_{mate}_fastqc.zip',sample=IDS,mate=[1,2]) if CFG['quantification']['trim'] else []
    output: OUT + '/qc/multiqc_report.html'
    params: image='tools', image_id=image_identity('tools'), root=OUT, outdir=OUT + '/qc'
    resources: mem_mb=2000
    log: OUT + '/logs/multiqc.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'multiqc {input.qc:q} {input.quant:q} {input.trimmed:q} --outdir {params.outdir:q} --filename multiqc_report.html --force > {log:q} 2>&1'
