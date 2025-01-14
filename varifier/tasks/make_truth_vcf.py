from varifier import truth_variant_finding, utils


def run(options):
    if options.truth_mask is None:
        mask = None
    else:
        mask = utils.load_mask_bed_file(options.truth_mask)

    truth_variant_finding.make_truth_vcf(
        options.ref_fasta,
        options.truth_fasta,
        options.outdir,
        options.flank_length,
        debug=options.debug,
        truth_mask=mask,
        max_ref_len=options.max_recall_ref_len,
        split_ref=options.split_ref,
        threads=options.cpus,
        maxmatch=not options.no_maxmatch,
        use_global_align=options.global_align,
        fix_truth_gap_lengths=options.sanitise_truth_gaps,
        hp_min_fix_length=options.hp_min_fix_length,
        fix_indel_max_length=options.indel_max_fix_length,
        global_align_min_coord=options.global_align_min_coord - 1,
        global_align_max_coord=options.global_align_max_coord - 1,
        ignore_non_acgt=not options.use_non_acgt,
        use_mafft=options.use_mafft,
    )
