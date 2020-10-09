#!/usr/bin/env python3

import pyfastaq

# These variants are at 0-based positions
variants = [
    (1000, "A", "T"),
    (1003, "T", "C"),
    (1004, "GGAGA", "G"),
    (1012, "G", "C"),
    (1014, "A", "AC"),
    (1018, "C", "G"),
    (1021, "T", "TAG"),
    (1023, "T", "TGC"),
    (1025, "G", "A"),
    (1027, "GC", "G"),
    (1032, "A", "T"),
    (1034, "A", "G"),
    (1037, "A", "C"),
    (1040, "G", "A"),
    (1043, "G", "A"),
    (1044, "G", "C"),
    (1045, "A", "T"),
    (1047, "C", "CG"),
]


ref = next(pyfastaq.sequences.file_reader("clustered_snp_indel.ref.fa"))

# Make the truth sequence by applying the variants. Start at the end and work
# backwards so indels don't mess up the coordinates.
truth_seq = list(ref.seq)
vcf_lines = []
for (position, ref_allele, alt_allele) in reversed(variants):
    vcf_lines.append("\t".join([ref.id, str(position + 1), ".", ref_allele, alt_allele, ".", "PASS", ".", "GT", "1/1"]))
    assert "".join(truth_seq[position:position + len(ref_allele)]) == ref_allele
    truth_seq[position:position + len(ref_allele)] = alt_allele

# varifier expects the variants to be sorted by position
vcf_lines.reverse()

with open("clustered_snp_indel.in.vcf", "w") as f:
    print("##fileformat=VCFv4.2", file=f)
    print('##FILTER=<ID=PASS,Description="All filters passed">', file=f)
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=f)
    print(f"##contig=<ID=ref,length={len(ref)}>", file=f)
    print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample", sep="\t", file=f)
    print(*vcf_lines, sep="\n", file=f)


# Add in insertion and deletion before and after the cluster, so that the
# probe mapping will have I and/or D in it (this case had to be removed
# from test_annotate_vcf_with_probe_mapping(), so put it in here)
truth_seq[990] = "TG"
truth_seq[1055] = ""

truth_seq = pyfastaq.sequences.Fasta("truth", "".join(truth_seq))
with open("clustered_snp_indel.truth.fa", "w") as f:
    print(truth_seq, file=f)

truth_seq.revcomp()
with open("clustered_snp_indel.truth.revcomp.fa", "w") as f:
    print(truth_seq, file=f)

