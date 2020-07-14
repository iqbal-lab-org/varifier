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
    )
