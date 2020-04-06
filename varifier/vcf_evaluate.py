import json
import os
import subprocess

from varifier import probe_mapping, recall, utils, vcf_stats


def _add_overall_precision_and_recall_to_summary_stats(summary_stats):
    # The recall came from the recall VCF. It contains all variants that
    # should be found. TP means it was found, FN means it wasn't found
    # (because these tags were added by probe mapping, which uses TP,FP).
    # So we either have FP or FN for the "wrong" variants from probe mapping.
    for prec_or_recall in "Precision", "Recall":
        fp_key = "FP" if prec_or_recall == "Precision" else "FN"
        for all_or_filt in "ALL", "FILT":
            d = summary_stats[prec_or_recall][all_or_filt]
            tp = d["TP"]["Count"]
            fp = d[fp_key]["Count"]
            total_calls = tp + fp
            if total_calls > 0:
                d[prec_or_recall] = round(tp / total_calls, 8)
            else:
                d[prec_or_recall] = 0


def evaluate_vcf(
    vcf_to_eval,
    vcf_ref_fasta,
    truth_ref_fasta,
    flank_length,
    outdir,
    debug=False,
    force=False,
    mask_bed_file=None,
):
    if force:
        subprocess.check_output(f"rm -rf {outdir}", shell=True)
    os.mkdir(outdir)

    # Mask if needed
    if mask_bed_file is not None:
        masked_vcf = os.path.join(outdir, "variants_to_eval.masked.vcf")
        utils.mask_vcf_file(vcf_to_eval, mask_bed_file, masked_vcf)
        vcf_to_eval = masked_vcf

    # Make VCF annotated with TP/FP for precision
    vcf_for_precision = os.path.join(outdir, "precision.vcf")
    if debug:
        map_outfile = f"{vcf_for_precision}.debug.map"
    else:
        map_outfile = None

    probe_mapping.annotate_vcf_with_probe_mapping(
        vcf_to_eval,
        vcf_ref_fasta,
        truth_ref_fasta,
        flank_length,
        vcf_for_precision,
        map_outfile=map_outfile,
    )

    # Make truth VCF annotated with TP/FP for recall
    recall_dir = os.path.join(outdir, "recall")
    vcf_for_recall_all, vcf_for_recall_filtered = recall.get_recall(
        vcf_ref_fasta, vcf_to_eval, truth_ref_fasta, recall_dir, debug=debug
    )
    if mask_bed_file is not None:
        utils.mask_vcf_file(vcf_for_recall_all, mask_bed_file, f"{vcf_for_recall_all}.masked.vcf")
        vcf_for_recall_all = f"{vcf_for_recall_all}.masked.vcf"
        utils.mask_vcf_file(vcf_for_recall_filtered, mask_bed_file, f"{vcf_for_recall_filtered}.masked.vcf")
        vcf_for_recall_filtered = f"{vcf_for_recall_filtered}.masked.vcf"

    # Gather stats and make plots
    per_record_recall_all = vcf_stats.per_record_stats_from_vcf_file(vcf_for_recall_all)
    per_record_recall_filtered = vcf_stats.per_record_stats_from_vcf_file(
        vcf_for_recall_filtered
    )
    recall_stats_all = vcf_stats.summary_stats_from_per_record_stats(
        per_record_recall_all, for_recall=True
    )
    recall_stats_filtered = vcf_stats.summary_stats_from_per_record_stats(
        per_record_recall_filtered, for_recall=True
    )
    recall_stats = {
        "ALL": recall_stats_all["ALL"],
        "FILT": recall_stats_filtered["ALL"],
    }

    per_record_precision = vcf_stats.per_record_stats_from_vcf_file(vcf_for_precision)
    precision_stats = vcf_stats.summary_stats_from_per_record_stats(
        per_record_precision
    )

    summary_stats = {"Recall": recall_stats, "Precision": precision_stats}
    _add_overall_precision_and_recall_to_summary_stats(summary_stats)

    summary_stats_json = os.path.join(outdir, "summary_stats.json")
    with open(summary_stats_json, "w") as f:
        json.dump(summary_stats, f, indent=2, sort_keys=True)
