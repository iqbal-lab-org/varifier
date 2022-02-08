from varifier import vcf_evaluate


def run(options):
    filter_pass = (
        None if options.filter_pass is None else set(options.filter_pass.split(","))
    )
    vcf_evaluate.evaluate_vcf(
        options.vcf_in,
        options.vcf_fasta,
        options.truth_fasta,
        options.flank_length,
        options.outdir,
        truth_vcf=options.truth_vcf,
        debug=options.debug,
        force=options.force,
        filter_pass=filter_pass,
        ref_mask_bed_file=options.ref_mask,
        truth_mask_bed_file=options.truth_mask,
        discard_ref_calls=not options.use_ref_calls,
        max_recall_ref_len=options.max_recall_ref_len,
        split_ref=options.split_ref,
        threads=options.cpus,
        maxmatch=not options.no_maxmatch,
        use_global_align=options.global_align,
        global_align_min_coord=options.global_align_min_coord - 1,
        global_align_max_coord=options.global_align_max_coord - 1,
        ignore_non_acgt=not options.use_non_acgt,
    )
