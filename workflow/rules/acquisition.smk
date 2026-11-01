rule validate_metadata:
    input:
        container_runner="workflow/scripts/container_exec.py",
        samples=CFG['samples'], archive=['metadata/source/ena_runs.tsv','metadata/source/sources.json'],
        code=['workflow/scripts/validate_inputs.py','workflow/scripts/common.py']
    output: OUT + '/validation/metadata.json'
    params:
        image='tools', image_id=image_identity('tools'), free=CFG['resources']['min_free_gib'],
        budget=CFG['resources']['download_limit_gib']
    resources: mem_mb=512
    log: OUT + '/logs/metadata.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'python workflow/scripts/validate_inputs.py --samples {input.samples:q} '
        '--output {output:q} --min-free-gib {params.free} --max-download-gib {params.budget} > {log:q} 2>&1'

rule download_fastq:
    input:
        container_runner="workflow/scripts/container_exec.py",
        validated=OUT + '/validation/metadata.json', samples=CFG['samples'], archive=['metadata/source/ena_runs.tsv','metadata/source/sources.json'],
        code=[ancient('workflow/scripts/download.py'),ancient('workflow/scripts/common.py')]
    output:
        fastq=protected(RAW + '/{sample}_{mate}.fastq.gz'),
        provenance=protected(RAW + '/{sample}_{mate}.fastq.gz.source.json')
    wildcard_constraints: mate='[12]', sample='|'.join(IDS)
    params:
        image='tools', image_id=image_identity('tools'), url=lambda w: BY_ID[w.sample]['fastq_' + w.mate],
        md5=lambda w: BY_ID[w.sample]['md5_' + w.mate],
        nbytes=lambda w: BY_ID[w.sample]['bytes_' + w.mate]
    resources: mem_mb=512, download_slots=1
    log: OUT + '/logs/download/{sample}_{mate}.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'python workflow/scripts/download.py --url {params.url:q} --output {output.fastq:q} '
        '--md5 {params.md5} --bytes {params.nbytes} > {log:q} 2>&1'

rule download_reference:
    input:
        container_runner="workflow/scripts/container_exec.py",
        validated=OUT + '/validation/metadata.json',
        code=[ancient('workflow/scripts/download.py'),ancient('workflow/scripts/common.py')]
    output:
        fasta=protected(REF + '/source/{kind}.gz'),
        provenance=protected(REF + '/source/{kind}.gz.source.json')
    wildcard_constraints: kind='cdna|ncrna|gtf|genome'
    params: image='tools', image_id=image_identity('tools'), url=lambda w: CFG['reference'][w.kind + '_url']
    resources: mem_mb=512, download_slots=1
    log: OUT + '/logs/reference/{kind}.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'python workflow/scripts/download.py --url {params.url:q} --output {output.fasta:q} > {log:q} 2>&1'

rule prepare_reference:
    input:
        container_runner="workflow/scripts/container_exec.py",
        cdna=REF + '/source/cdna.gz', ncrna=REF + '/source/ncrna.gz', gtf=REF + '/source/gtf.gz',
        genome=[REF + '/source/genome.gz'] if CFG['reference']['decoys'] else [],
        code=['workflow/scripts/prepare_reference.py','workflow/scripts/common.py']
    output:
        fasta=REF + '/prepared/gentrome.fa', map=REF + '/prepared/tx2gene.tsv',
        decoys=REF + '/prepared/decoys.txt', provenance=REF + '/prepared/reference.json'
    params:
        image='tools', image_id=image_identity('tools'), outdir=REF + '/prepared',
        genome=('--genome ' + shlex.quote(REF + '/source/genome.gz')) if CFG['reference']['decoys'] else ''
    resources: mem_mb=4000
    log: OUT + '/logs/reference/prepare.log'
    shell:
        '{RUN} --image {params.image} --expected-id {params.image_id} --memory-mb {resources.mem_mb} -- '
        'python workflow/scripts/prepare_reference.py --cdna {input.cdna:q} --ncrna {input.ncrna:q} '
        '--gtf {input.gtf:q} --outdir {params.outdir:q} {params.genome} > {log:q} 2>&1'
