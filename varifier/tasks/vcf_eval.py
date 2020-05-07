from varifier import vcf_evaluate


def run(options):
    vcf_evaluate.evaluate_vcf(
        options.vcf_in,
        options.vcf_fasta,
        options.truth_fasta,
        options.flank_length,
        options.outdir,
        debug=options.debug,
        force=options.force,
        mask_bed_file=options.mask,
        discard_ref_calls=not options.use_ref_calls,
    )
