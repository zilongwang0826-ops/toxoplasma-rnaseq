#!/usr/bin/env bash
# Read-only audit. No remote installations or workspace writes.
set -u
id; hostname; uname -srm; cat /etc/os-release; pwd
nproc; free -h; df -hT / /home /tmp
for tool in docker apptainer singularity snakemake nextflow python3 R Rscript mamba conda git java salmon fastqc fastp multiqc curl rsync; do command -v "$tool" || true; done
docker version
docker info --format 'ServerVersion={{.ServerVersion}} Architecture={{.Architecture}} CPUs={{.NCPU}} Memory={{.MemTotal}} Root={{.DockerRootDir}}'
docker image ls --format '{{.Repository}}:{{.Tag}} {{.Size}}'
ls -ld "$HOME" /var/lib/miniforge /var/run/docker.sock
python3 --version; git --version; java -version 2>&1
/var/lib/miniforge/bin/mamba --version
/var/lib/miniforge/bin/mamba env list
/var/lib/miniforge/bin/mamba config list channel_priority envs_dirs pkgs_dirs
if [[ -f /etc/motd ]]; then cat /etc/motd; fi
for url in 'https://www.ncbi.nlm.nih.gov/bioproject/PRJNA483261' 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA483261&result=read_run&fields=run_accession&format=tsv' 'https://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz'; do
 printf '%s\n' "$url"
 curl -fsSI --connect-timeout 10 --max-time 30 "$url" | sed -n '1p;/^[Cc]ontent-[Ll]ength:/p;/^[Ll]ast-[Mm]odified:/p'
done
