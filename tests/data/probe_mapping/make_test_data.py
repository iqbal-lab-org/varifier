#!/usr/bin/env python3

# This script makes a VCF and corresponding reference FASTA, and a truth
# FASTA. Both FASTAs have a single sequence of length 350bp.
# Differences are, using 1-based coords:
# - The first 20bp are different, so probes won't map there
# - SNP at position 40 (so probes won't map with total length) A->G, where
#   bases before/after are not A or G
# - SNP at position 60 with A->G, but 59-61 is AAA
# - SNP at position 80 with T->C, but 79-81 is CCC
# - Three SNPs in a row at 109, 110, 111, GCA -> CTG
# - Deletion GT -> G at position 140 (relative to ref)
# - Insertion A -> AT at position 160
# - complex variant ATGCGTATCG -> AGGTACGAGGCAT at position 300
# - Last 20bp are different, so probes won't map

import pyfastaq
import random

random.seed(42)

acgt = {"A", "C", "G", "T"}


def random_seq(length):
    return [random.choice(["A", "C", "G", "T"]) for _ in range(length)]


def add_fp_snps(positions, vcf_lines, ref_seq, truth_seq):
    for i in positions:
        truth_base = truth_seq[i]
        for base in acgt.difference(truth_base):
            if ref_seq[i] != base:
                vcf_lines.append([ "ref", i + 1, ".", ref_seq[i], base, ".", "PASS", ".", "GT:EXPECT", "1/1:FP"])


vcf_lines = []
ref_seq = random_seq(17) + ["A", "G", "T"]
truth_seq = random_seq(17) + ["T", "C", "A"]

ref_seq += random_seq(17) + ["G", "T"]
truth_seq += ref_seq[-19:]
assert len(ref_seq) == len(truth_seq)


# SNP at 0-based position 39:
ref_seq.append("A")
truth_seq.append("G")
assert len(ref_seq) == len(truth_seq) == 40
vcf_lines.append(["ref", 40, ".", ref_seq[39], truth_seq[39], ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
fp_positions = list(range(36, 43))

# SNP at 0-based position 59
ref_seq += ["C", "A"] + random_seq(16)
truth_seq += ref_seq[-18:]
ref_seq += ["A", "A", "A"]
truth_seq += ["A", "G", "A"]
assert len(ref_seq) == len(truth_seq) == 61
vcf_lines.append(["ref", 60, ".", ref_seq[59], truth_seq[59], ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
fp_positions.extend(range(56, 63))

# SNP at 0-based position 79, with T->C, and 78-80 is CTC
ref_seq += random_seq(17)
truth_seq += ref_seq[-17:]
ref_seq += ["C", "T", "C"]
truth_seq += ["C", "C", "C"]
assert len(ref_seq) == len(truth_seq) == 81
vcf_lines.append(["ref", 80, ".", ref_seq[79], truth_seq[79], ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
fp_positions.extend(range(76, 83))

# Three SNPs in a row at 108, 109, 110, GCA -> CTG
ref_seq += random_seq(27) + ["G", "C", "A"] + random_seq(28)
truth_seq += ref_seq[-58:]
truth_seq[108] = "C"
truth_seq[109] = "T"
truth_seq[110] = "G"
assert len(truth_seq) == len(ref_seq) == 139
vcf_lines.append(["ref", 109, ".", ref_seq[108], truth_seq[108], ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
vcf_lines.append(["ref", 110, ".", ref_seq[109], truth_seq[109], ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
vcf_lines.append(["ref", 111, ".", ref_seq[110], truth_seq[110], ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
fp_positions.extend(range(105, 114))

# - Deletion GT -> G at position 139
ref_seq += ["G", "T"] + random_seq(18)
truth_seq += ["G"] + ref_seq[-18:]
vcf_lines.append(["ref", 140, ".", "".join(ref_seq[139:141]), truth_seq[139], ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
vcf_lines.append(["ref", 140, ".", "".join(ref_seq[139:141]), "AC", ".", "PASS", ".", "GT:EXPECT", "1/1:FP"])
fp_positions.extend(range(136, 140))
assert len(ref_seq) == 159
assert len(truth_seq) == 158

# - Insertion A -> AT at position 159
ref_seq += ["A"] + random_seq(139)
truth_seq += ["A", "T"] + ref_seq[-139:]
vcf_lines.append(["ref", 160, ".", ref_seq[159], "".join(truth_seq[158:160]), ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
vcf_lines.append(["ref", 160, ".", ref_seq[159], "AG", ".", "PASS", ".", "GT:EXPECT", "1/1:Partial_TP"])

# These two make sure we hit lines in the code where the cigar string has I or D
# while the probe position that is tracked is inside the allele of interest.
# The other variants in this test file don't fo that.
vcf_lines.append(["ref", 161, ".", "".join(ref_seq[160:163]), "ACGTACGGGTGGTGTGTTTGAAAGATAG", ".", "PASS", ".", "GT:EXPECT", "1/1:Partial_TP"])
vcf_lines.append(["ref", 161, ".", "".join(ref_seq[160:175]), "ACG", ".", "PASS", ".", "GT:EXPECT", "1/1:Partial_TP"])

fp_positions.extend(range(159, 164))
assert len(ref_seq) == len(truth_seq) == 299

# - complex variant ATGCGTATCG -> AGGTACGAGA at position 299
complex_ref = "ATGCGTATCG"
complex_alt = "AGGTACGAGA"
complex_alt_half_right = "AGGTAACCCC"
ref_seq += list(complex_ref) + random_seq(21)
truth_seq += list(complex_alt) + ref_seq[-21:]
ref_seq += random_seq(100)
truth_seq += random_seq(100)
vcf_lines.append(["ref", 300, ".", complex_ref, complex_alt, ".", "PASS", ".", "GT:EXPECT", "1/1:TP"])
vcf_lines.append(["ref", 300, ".", complex_ref, complex_alt_half_right, ".", "PASS", ".", "GT:EXPECT", "1/1:Partial_TP"])
assert len(ref_seq) == len(truth_seq) == 430
fp_positions.extend(range(297, 304))

# Check that FP deletion does not get called as TP. Here's a toy example
# to explain why:
#   ref   AGTGCAGT
#   call: 3 TG -> T
#   If call is wrong:
#   truth:      AGTGCAGT
#   alt-probe:  AGT-CAGT
#   ref-probe:  AGTGCAGT
# Note that if the call is wrong, the T of the alt probe still maps, so it
# looks correct.
ref_seq[229] = truth_seq[229] = "T"
ref_seq[230] = truth_seq[230] = "G"
vcf_lines.append(["ref", 230, ".", "TG", "T", ".", "PASS", ".", "GT:EXPECT", "1/1:FP"])


# SNP so near the end that probe won't map
vcf_lines.append(["ref", 400, ".", ref_seq[399], "AAAAA", ".", "PASS", ".", "GT:EXPECT", "1/1:FP_PROBE_UNMAPPED"])

# Add in records that can't be evaluated for various reasons
ref_seq[249] = "G"
truth_seq[249] = "C"
vcf_lines.append(["ref", 250, ".", "G", "C", ".", "FILTER_X", ".", "GT:EXPECT", "1/1:TP"])
vcf_lines.append(["ref", 250, ".", "G", "C", ".", "MISMAPPED_UNPLACEABLE", ".", "GT:EXPECT", "1/1:TP"])
vcf_lines.append(["ref", 250, ".", "G", "C", ".", "PASS", ".", "GT:EXPECT", "0/1:VFR_FILTER_CANNOT_USE_GT"])
vcf_lines.append(["ref", 250, ".", "G", "C", ".", "PASS", ".", "GT:EXPECT", "./.:VFR_FILTER_CANNOT_USE_GT"])
vcf_lines.append(["ref", 250, ".", "C", "G", ".", "PASS", ".", "GT:EXPECT", "1/1:VFR_FILTER_REF_STRING_MISMATCH"])
vcf_lines.append(["ref", 250, ".", "G", "C", ".", "PASS", ".", "EXPECT", "VFR_FILTER_NO_GT"])
vcf_lines.append(["ref", 250, ".", "G", ".", ".", "PASS", ".", "GT:EXPECT", "1/1:VFR_FILTER_NO_ALTS"])
vcf_lines.append(["ref", 250, ".", ".", "G", ".", "PASS", ".", "GT:EXPECT", "1/1:VFR_FILTER_NO_REF_SEQ"])


add_fp_snps(fp_positions, vcf_lines, ref_seq, truth_seq)


with open("annotate_vcf_with_probe_mapping.in.vcf", "w") as f:
    print(f"##contig=<ID=ref,length={len(ref_seq)}>", file=f)
    print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "sample", sep="\t", file=f)
    vcf_lines.sort()
    for i, l in enumerate(vcf_lines):
        l[2] = i + 1
        print(*l, sep="\t", file=f)

with open("annotate_vcf_with_probe_mapping.ref.fa", "w") as f:
    fa = pyfastaq.sequences.Fasta("ref", "".join(ref_seq))
    print(fa, file=f)

with open("annotate_vcf_with_probe_mapping.truth.fa", "w") as f:
    fa = pyfastaq.sequences.Fasta("truth", "".join(truth_seq))
    print(fa, file=f)

fa.revcomp()

with open("annotate_vcf_with_probe_mapping.truth.revcomp.fa", "w") as f:
    print(fa, file=f)
