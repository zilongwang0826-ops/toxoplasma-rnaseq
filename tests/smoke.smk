RUN='python3 workflow/scripts/container_exec.py'
BASE='results/tests/salmon'
rule all:
    input: 'results/tests/locus-sensitivity/verified.txt','results/tests/go-integration/verified.json','results/tests/python/verified.json',BASE+'/verified.json', 'results/tests/container-analysis/integration-status.tsv', BASE+'/qc/sample_1_fastqc.zip',BASE+'/qc/sample_2_fastqc.zip'

rule locus_sensitivity_tests:
    input: 'tests/test_locus_sensitivity.R','workflow/scripts/locus_sensitivity.R','workflow/scripts/analyze.R','workflow/scripts/enrichment.R'
    output: 'results/tests/locus-sensitivity/verified.txt'
    threads: 2
    shell: 'mkdir -p results/tests/locus-sensitivity && {RUN} --image statistics --threads {threads} --memory-mb 1000 -- Rscript --vanilla tests/test_locus_sensitivity.R > {output:q}'
rule generate:
    input: 'tests/make_smoke.py'
    output: ref=BASE+'/transcripts.fa', r1=BASE+'/sample_1.fastq.gz',r2=BASE+'/sample_2.fastq.gz'
    threads: 2
    params: image='tools'
    shell: '{RUN} --image {params.image} --threads {threads} --memory-mb 1000 -- python {input:q} --outdir {BASE}'
rule qc:
    input: BASE+'/sample_{mate}.fastq.gz'
    output: zip=BASE+'/qc/sample_{mate}_fastqc.zip',html=BASE+'/qc/sample_{mate}_fastqc.html'
    params: image='tools'
    threads: 2
    shell: '{RUN} --image {params.image} --threads {threads} --memory-mb 2000 -- fastqc --threads {threads} --outdir {BASE}/qc {input:q}'
rule index:
    input: BASE+'/transcripts.fa'
    output: directory(BASE+'/index')
    params: image='tools'
    threads: 2
    shell: '{RUN} --image {params.image} --threads {threads} --memory-mb 2000 -- salmon index -t {input:q} -i {output:q} -p {threads} --keepDuplicates -k 31'
rule quantify:
    input: index=BASE+'/index',r1=BASE+'/sample_1.fastq.gz',r2=BASE+'/sample_2.fastq.gz'
    output: BASE+'/quant/quant.sf'
    params: image='tools'
    threads: 2
    shell: '{RUN} --image {params.image} --threads {threads} --memory-mb 2000 -- salmon quant -i {input.index:q} -l IU -1 {input.r1:q} -2 {input.r2:q} -p {threads} --validateMappings -o {BASE}/quant'
rule assert_quant:
    input: quant=BASE+'/quant/quant.sf',code='tests/assert_smoke.py'
    output: BASE+'/verified.json'
    params: image='tools'
    shell: '{RUN} --image {params.image} --memory-mb 1000 -- python {input.code:q}'
rule statistics:
    input: 'tests/test_analysis.R','workflow/scripts/analyze.R','workflow/scripts/enrichment.R','workflow/scripts/compare_published.R'
    output: 'results/tests/container-analysis/integration-status.tsv'
    params: image='statistics'
    threads: 2
    shell: '{RUN} --image {params.image} --threads {threads} --memory-mb 4000 -- Rscript --vanilla tests/test_analysis.R results/tests/container-analysis'

rule input_integrity:
    input: 'tests/run_python_tests.py','tests/test_inputs.py','workflow/scripts/validate_inputs.py','workflow/scripts/common.py','workflow/scripts/download.py','workflow/scripts/prepare_reference.py','workflow/scripts/quant_qc.py'
    output: 'results/tests/python/verified.json'
    params: image='tools'
    threads: 2
    shell: '{RUN} --image {params.image} --threads {threads} --memory-mb 1000 -- python tests/run_python_tests.py'

rule go_integration:
    input: 'tests/test_go_integration.R','workflow/scripts/analyze.R','workflow/scripts/enrichment.R'
    output: 'results/tests/go-integration/verified.json'
    params: image='statistics'
    threads: 2
    shell: '{RUN} --image {params.image} --threads {threads} --memory-mb 4000 -- Rscript --vanilla tests/test_go_integration.R'
