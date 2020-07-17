import operator
import os

import pyfastaq
from cluster_vcf_records import vcf_file_read

from varifier import probe_mapping, truth_variant_finding, utils


def _vcf_file_to_dict(vcf_file, pass_only=True):
    """Loads VCF file. Returns a dictionary of sequence name -> sorted list
    by position of variants"""
    records = {}
    wanted_format = {"PASS"}
    if not pass_only:
        wanted_format.add("FAIL_BUT_TEST")

    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(vcf_file)
    for record in vcf_records:
        if record.FORMAT["VFR_FILTER"] not in wanted_format:
            continue

        if record.CHROM not in records:
            records[record.CHROM] = []
        records[record.CHROM].append(record)

    for l in records.values():
        l.sort(key=operator.attrgetter("POS"))

    return records


def apply_variants_to_genome(ref_fasta, vcf_file, out_fasta, pass_only=True):
    """Takes the variants in vcf_file, and applies them to the associated
    reference genome in ref_fasta. Writes a new file out_fasta that has those
    variants applied"""
    ref_sequences = utils.file_to_dict_of_seqs(ref_fasta)
    vcf_dict = _vcf_file_to_dict(vcf_file, pass_only=pass_only)
    with open(out_fasta, "w") as f:
        for ref_name, vcf_records in sorted(vcf_dict.items()):
            old_seq = ref_sequences[ref_name]
            new_seq = list(old_seq.seq)
            # Applying indels messes up the coords of any subsequent variant,
            # so start at the end and work backwards
            for vcf_record in reversed(vcf_records):
                genotype = set(vcf_record.FORMAT["GT"].split("/"))
                assert len(genotype) == 1
                allele_index = int(genotype.pop())
                if allele_index == 0:
                    continue
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
        # Make truth VCF. This only depends on ref_fasta and truth_fasta, not
        # on VCF to test. In particular, is independent of whether or not
        # were using all records in vcf_to_test, or PASS records only. This means
        # only need to make one truth VCF, which can be used for both cases.
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

    vcfs_out = {}
    for all_or_filt in "ALL", "FILT":
        run_outdir = os.path.join(outdir, all_or_filt)
        os.mkdir(run_outdir)
        mutated_ref_fasta = os.path.join(run_outdir, "00.ref_with_mutations_added.fa")
        apply_variants_to_genome(
            ref_fasta, vcf_to_test, mutated_ref_fasta, pass_only=all_or_filt == "FILT"
        )

        # For each record in the truth VCF, make a probe and map to the mutated genome
        vcfs_out[all_or_filt] = os.path.join(
            run_outdir, "02.truth.probe_mapped_to_mutated_genome.vcf"
        )
        map_outfile = (
            os.path.join(run_outdir, "02.probe_map_debug.txt") if debug else None
        )
        probe_mapping.annotate_vcf_with_probe_mapping(
            truth_vcf,
            ref_fasta,
            mutated_ref_fasta,
            flank_length,
            vcfs_out[all_or_filt],
            map_outfile=map_outfile,
        )
    return vcfs_out["ALL"], vcfs_out["FILT"]
