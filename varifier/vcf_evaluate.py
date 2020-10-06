import json
import logging
import os
import subprocess

from cluster_vcf_records import vcf_file_read, vcf_record

from varifier import probe_mapping, recall, utils, vcf_stats


def _add_overall_precision_and_recall_to_summary_stats(summary_stats):
    # The recall came from the recall VCF. It contains all variants that
    # should be found. TP means it was found, FN means it wasn't found
    # (because these tags were added by probe mapping, which uses TP,FP).
    # So we either have FP or FN for the "wrong" variants from probe mapping.
    for prec_or_recall in "Precision", "Recall":
        fp_key = "FP" if prec_or_recall == "Precision" else "FN"
        d = summary_stats[prec_or_recall]
        tp = d["TP"]["Count"]
        tp_frac = d["TP"]["SUM_ALLELE_MATCH_FRAC"] + d[fp_key]["SUM_ALLELE_MATCH_FRAC"]
        fp = d[fp_key]["Count"]
        total_calls = tp + fp
        if total_calls > 0:
            d[prec_or_recall] = round(tp / total_calls, 8)
            d[f"{prec_or_recall}_frac"] = round(tp_frac / total_calls, 8)
        else:
            d[prec_or_recall] = 0
            d[f"{prec_or_recall}_frac"] = 0

        if d["EDIT_DIST_COUNTS"]["denominator"] > 0:
            d[f"{prec_or_recall}_edit_dist"] = round(
                d["EDIT_DIST_COUNTS"]["numerator"]
                / d["EDIT_DIST_COUNTS"]["denominator"],
                8,
            )
        else:
            d[f"{prec_or_recall}_edit_dist"] = 0


def _filter_vcf(
    infile,
    outfile_keep,
    outfile_exclude,
    ref_seqs,
    filter_pass=None,
    keep_ref_calls=False,
):
    counts = {
        "filter_fail": 0,
        "heterozygous": 0,
        "no_genotype": 0,
        "ref_call": 0,
        "other": 0,
    }
    with vcf_file_read.open_vcf_file_for_reading(infile) as f_in, open(
        outfile_keep, "w"
    ) as f_out_keep, open(outfile_exclude, "w") as f_out_exclude:
        for line in f_in:
            if line.startswith("#"):
                print(line, end="", file=f_out_keep)
                print(line, end="", file=f_out_exclude)
                continue

            record = vcf_record.VcfRecord(line)
            exclude_reason = None

            filter_is_dot_and_fails = (
                filter_pass is not None
                and len(record.FILTER) == 0
                and "." not in filter_pass
            )
            filter_not_dot_and_fails = (
                filter_pass is not None
                and len(record.FILTER) > 0
                and record.FILTER.isdisjoint(filter_pass)
            )
            filter_fails = filter_is_dot_and_fails or filter_not_dot_and_fails

            if "MISMAPPED_UNPLACEABLE" in record.FILTER or filter_fails:
                exclude_reason = "filter_fail"
            elif len(record.ALT) == 0 or record.ALT == ["."]:
                exclude_reason = "other"
            elif record.FORMAT is None:
                exclude_reason = "no_genotype"
            elif record.REF in [".", ""]:
                exclude_reason = "other"
            if (
                ref_seqs[record.CHROM][record.POS : record.POS + len(record.REF)]
                != record.REF
            ):
                exclude_reason = "other"

            if exclude_reason is None:
                gt = set(record.FORMAT.get("GT", ".").split("/"))
                if len(gt) > 1:
                    exclude_reason = "heterozygous"
                elif "." in gt:
                    exclude_reason = "no_genotype"
                elif not keep_ref_calls and "0" in gt:
                    exclude_reason = "ref_call"
                else:
                    print(record, file=f_out_keep)

            if exclude_reason is not None:
                record.set_format_key_value("VFR_EXCLUDE_REASON", exclude_reason)
                print(record, file=f_out_exclude)
                counts[exclude_reason] += 1

    return counts


def evaluate_vcf(
    vcf_to_eval,
    vcf_ref_fasta,
    truth_ref_fasta,
    flank_length,
    outdir,
    truth_vcf=None,
    debug=False,
    force=False,
    ref_mask_bed_file=None,
    truth_mask_bed_file=None,
    discard_ref_calls=True,
    max_recall_ref_len=None,
    filter_pass=None,
):
    if force:
        subprocess.check_output(f"rm -rf {outdir}", shell=True)
    os.mkdir(outdir)

    # Mask if needed
    if ref_mask_bed_file is None:
        vcf_to_filter = vcf_to_eval
    else:
        logging.info("Masking VCF...")
        masked_vcf = os.path.join(outdir, "variants_to_eval.masked.vcf")
        utils.mask_vcf_file(vcf_to_eval, ref_mask_bed_file, masked_vcf)
        vcf_to_filter = masked_vcf
        logging.info("Masked VCF")

    vcf_ref_seqs = utils.file_to_dict_of_seqs(vcf_ref_fasta)
    filtered_vcf = os.path.join(outdir, "variants_to_eval.filtered.vcf")
    excluded_vcf = os.path.join(outdir, "variants_to_eval.excluded.vcf")
    logging.info("Filtering VCF...")
    filtered_counts = _filter_vcf(
        vcf_to_filter,
        filtered_vcf,
        excluded_vcf,
        vcf_ref_seqs,
        filter_pass=filter_pass,
        keep_ref_calls=not discard_ref_calls,
    )
    logging.info("Filtering VCF done")

    vcf_for_precision = os.path.join(outdir, "precision.vcf")
    map_outfile = f"{vcf_for_precision}.debug.map" if debug else None
    if truth_mask_bed_file is None:
        truth_mask = None
    else:
        truth_mask = utils.load_mask_bed_file(truth_mask_bed_file)

    logging.info("Annotating VCF with TP/FP for precision...")
    probe_mapping.annotate_vcf_with_probe_mapping(
        filtered_vcf,
        vcf_ref_fasta,
        truth_ref_fasta,
        flank_length,
        vcf_for_precision,
        map_outfile=map_outfile,
        use_ref_calls=not discard_ref_calls,
        truth_mask=truth_mask,
    )
    logging.info("Annotatiing VCF with with TP/FP for precision done")

    logging.info("Calculating recall...")
    recall_dir = os.path.join(outdir, "recall")
    vcf_for_recall = recall.get_recall(
        vcf_ref_fasta,
        filtered_vcf,
        recall_dir,
        flank_length,
        debug=debug,
        truth_fasta=truth_ref_fasta if truth_vcf is None else None,
        truth_vcf=truth_vcf,
        truth_mask=truth_mask,
        max_ref_len=max_recall_ref_len,
    )
    if ref_mask_bed_file is not None:
        logging.info("Masking recall VCF...")
        utils.mask_vcf_file(
            vcf_for_recall, ref_mask_bed_file, f"{vcf_for_recall}.masked.vcf"
        )
        vcf_for_recall = f"{vcf_for_recall}.masked.vcf"
        logging.info("Masking recall VCF done")
    logging.info("Recall calculation done")

    # Gather stats and make plots
    logging.info("Gathering stats...")
    per_record_recall = vcf_stats.per_record_stats_from_vcf_file(vcf_for_recall)
    recall_stats = vcf_stats.summary_stats_from_per_record_stats(
        per_record_recall, for_recall=True
    )

    per_record_precision = vcf_stats.per_record_stats_from_vcf_file(vcf_for_precision)
    precision_stats = vcf_stats.summary_stats_from_per_record_stats(
        per_record_precision
    )

    summary_stats = {"Recall": recall_stats, "Precision": precision_stats}
    _add_overall_precision_and_recall_to_summary_stats(summary_stats)
    summary_stats["Excluded_record_counts"] = filtered_counts

    summary_stats_json = os.path.join(outdir, "summary_stats.json")
    with open(summary_stats_json, "w") as f:
        json.dump(summary_stats, f, indent=2, sort_keys=True)

    logging.info(f"Done. Results written to {summary_stats_json}")
