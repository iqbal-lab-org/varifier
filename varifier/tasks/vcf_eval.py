from varifier import vcf_evaluate


def run(options):
    vcf_evaluate.evaluate_vcf(
        options.vcf_in,
        options.vcf_fasta,
        options.truth_fasta,
        options.flank_length,
        options.outdir,
        truth_vcf=options.truth_vcf,
        debug=options.debug,
        force=options.force,
        ref_mask_bed_file=options.ref_mask,
        truth_mask_bed_file=options.truth_mask,
        discard_ref_calls=not options.use_ref_calls,
        max_recall_ref_len=options.max_recall_ref_len,
    )
