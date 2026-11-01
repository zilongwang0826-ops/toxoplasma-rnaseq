"""Build an Ensembl transcript/decoy FASTA and a verified transcript-to-gene map."""
import argparse
import csv
import gzip
import re
from pathlib import Path
from common import sha256, write_json, utc_now


def fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        header, seq = None, []
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq)


def prepare(cdna, ncrna, gtf, outdir, genome=None):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    genes, names = {}, {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            tid, gid = attrs["transcript_id"], attrs["gene_id"]
            if tid in genes and genes[tid] != gid:
                raise ValueError("Conflicting transcript-to-gene mapping: " + tid)
            genes[tid], names[tid] = gid, attrs.get("gene_name", "")
    seen, duplicates, records = {}, 0, []
    target = out / "gentrome.fa"
    with target.open("w") as f:
        for path in (cdna, ncrna):
            for header, seq in fasta(path):
                versioned = header.split()[0]
                tid = versioned.split(".")[0]
                if tid not in genes:
                    raise ValueError("Transcript missing in matching GTF: " + tid)
                header_gene = re.search(r"(?:^|\s)gene:([^\s]+)", header)
                if header_gene and header_gene.group(1).split(".")[0] != genes[tid]:
                    raise ValueError("FASTA/GTF gene identity mismatch: " + tid)
                if not seq or re.search("[^ACGTNacgtn]", seq):
                    raise ValueError("Invalid transcript sequence: " + tid)
                if tid in seen:
                    if seen[tid] != seq:
                        raise ValueError("Conflicting duplicated transcript sequence: " + tid)
                    duplicates += 1
                    continue
                seen[tid] = seq
                f.write(">" + versioned + "\n" + seq + "\n")
                records.append((versioned, genes[tid], names[tid]))
        # Stream the genome to avoid holding genomic sequences in memory.
        decoys = []
        if genome:
            with gzip.open(genome, "rt") as source:
                for line in source:
                    if line.startswith(">"):
                        chrom = line[1:].split()[0]
                        if chrom in seen or chrom in decoys:
                            raise ValueError("Duplicate decoy identifier")
                        decoys.append(chrom)
                        f.write(">" + chrom + "\n")
                    else:
                        f.write(line)
    if not records:
        raise ValueError("Empty transcript reference")
    with (out / "tx2gene.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["transcript_id", "gene_id", "gene_name"])
        w.writerows(records)
    (out / "decoys.txt").write_text("".join(x + "\n" for x in decoys))
    write_json(out / "reference.json", {"created_at": utc_now(), "transcripts": len(records),
        "genes": len({r[1] for r in records}), "duplicate_identical_transcripts_removed": duplicates,
        "decoy_sequences": len(decoys), "input_sha256": {str(x): sha256(x) for x in (cdna,ncrna,gtf,genome) if x},
        "gentrome_sha256": sha256(target), "tx2gene_sha256": sha256(out / "tx2gene.tsv")})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("cdna", "ncrna", "gtf", "outdir"):
        p.add_argument("--" + key, required=True)
    p.add_argument("--genome")
    a = p.parse_args()
    prepare(a.cdna, a.ncrna, a.gtf, a.outdir, a.genome)
