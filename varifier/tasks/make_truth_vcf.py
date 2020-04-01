from varifier import truth_variant_finding


def run(options):
    truth_variant_finding.make_truth_vcf(
        options.ref_fasta, options.truth_fasta, options.outdir, debug=options.debug
    )
