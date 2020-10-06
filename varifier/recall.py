import logging
import operator
import os

import pyfastaq
from cluster_vcf_records import vcf_file_read

from varifier import probe_mapping, truth_variant_finding, utils


def _vcf_file_to_dict(vcf_file):
    """Loads VCF file. Returns a dictionary of sequence name -> sorted list
    by position of variants"""
    records = {}

    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(vcf_file)
    for record in vcf_records:
        if record.CHROM not in records:
            records[record.CHROM] = []
        records[record.CHROM].append(record)

    for l in records.values():
        l.sort(key=operator.attrgetter("POS"))

    return records


def apply_variants_to_genome(ref_fasta, vcf_file, out_fasta):
    """Takes the variants in vcf_file, and applies them to the associated
    reference genome in ref_fasta. Writes a new file out_fasta that has those
    variants applied"""
    ref_sequences = utils.file_to_dict_of_seqs(ref_fasta)
    vcf_dict = _vcf_file_to_dict(vcf_file)
    with open(out_fasta, "w") as f:
        for ref_name, vcf_records in sorted(vcf_dict.items()):
            old_seq = ref_sequences[ref_name]
            new_seq = list(old_seq.seq)
            previous_ref_start = None
            # Applying indels messes up the coords of any subsequent variant,
            # so start at the end and work backwards
            for vcf_record in reversed(vcf_records):
                genotype = set(vcf_record.FORMAT["GT"].split("/"))
                assert len(genotype) == 1
                allele_index = int(genotype.pop())
                if allele_index == 0:
                    continue

                # Some tools report two (or more) variants that overlap.
                # No clear "right" option here.
                # If the current record overlaps the previous one, ignore it.
                # We could try to be cleverer about this (take best records
                # based on likelihoods or whatever else), but every tool is
                # different so no sane consistent way of doing this across tools
                if (
                    previous_ref_start is not None
                    and vcf_record.ref_end_pos() >= previous_ref_start
                ):
                    logging.warn(
                        f"Skipping this record when calculating recall because it overlaps another record: {vcf_record}"
                    )
                    continue

                previous_ref_start = vcf_record.POS
                allele = vcf_record.ALT[allele_index - 1]
                start, end = vcf_record.POS, vcf_record.ref_end_pos() + 1
                assert old_seq[start:end] == "".join(new_seq[start:end])
                new_seq[start:end] = [allele]
            new_seq = pyfastaq.sequences.Fasta(f"{ref_name}.mutated", "".join(new_seq))
            print(new_seq, file=f)


def get_recall(
    ref_fasta,
    vcf_to_test,
    outdir,
    flank_length,
    truth_fasta=None,
    truth_vcf=None,
    debug=False,
    truth_mask=None,
    max_ref_len=None,
):
    os.mkdir(outdir)

    if truth_vcf is None:
        assert truth_fasta is not None
        truth_outdir = os.path.join(outdir, "truth_vcf")
        truth_vcf = truth_variant_finding.make_truth_vcf(
            ref_fasta,
            truth_fasta,
            truth_outdir,
            flank_length,
            debug=debug,
            truth_mask=truth_mask,
            max_ref_len=max_ref_len,
        )
    else:
        assert truth_fasta is None

    mutated_ref_fasta = os.path.join(outdir, "ref_with_mutations_added.fa")
    apply_variants_to_genome(ref_fasta, vcf_to_test, mutated_ref_fasta)

    vcf_out = os.path.join(outdir, "recall.vcf")
    map_outfile = os.path.join(outdir, "probe_map_debug.txt") if debug else None
    probe_mapping.annotate_vcf_with_probe_mapping(
        truth_vcf,
        ref_fasta,
        mutated_ref_fasta,
        flank_length,
        vcf_out,
        map_outfile=map_outfile,
    )
    return vcf_out
