#!/usr/bin/env python3

import argparse
import logging
import varifier


def main(args=None):
    parser = argparse.ArgumentParser(
        prog="varifier",
        usage="varifier <command> <options>",
        description="varifier: variant call adjudication",
    )

    parser.add_argument("--version", action="version", version=varifier.__version__)
    parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )

    subparsers = parser.add_subparsers(title="Available commands", help="", metavar="")

    # ----------- general options common to all tasks ------------------------
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument(
        "--debug",
        help="More verbose logging, and less file cleaning",
        action="store_true",
    )
    common_parser.add_argument(
        "--flank_length",
        help="Length of sequence to add either side of variant when making probe sequences [%(default)s]",
        type=int,
        default=100,
        metavar="INT",
    )
    common_parser.add_argument(
        "--max_recall_ref_len",
        help="Do not include variants where REF length is more than this number. Default is no limit",
        type=int,
        metavar="INT",
    )
    common_parser.add_argument(
        "--split_ref",
        help="When using MUMmer, split the ref genome into one file per sequence, and run MUMmer on each split. Experimental - should improve run time for big genomes",
        action="store_true",
    )
    common_parser.add_argument(
        "--no_maxmatch",
        help="When using nucmer to get expected calls for recall, do not use the --maxmatch option. May reduce sensitivity to finding all variants",
        action="store_true",
    )
    common_parser.add_argument(
        "--global_align",
        help="Only use this with small genomes (ie virus) that have one sequence each, in the same orientation. Instead of using minimap2/nucmer to find variants, do a global alignment for greater accuracy",
        action="store_true",
    )
    common_parser.add_argument(
        "--global_align_min_coord",
        help="Only used if also using --global_align. Do not output variants where the REF allele starts before the given (1-based) coordinate. When running vcf_eval, only applies to recall, not precision [%(default)s]",
        type=int,
        metavar="INT",
        default=1,
    )
    common_parser.add_argument(
        "--global_align_max_coord",
        help="Only used if also using --global_align. Do not output variants where the REF allele ends after the given (1-based) coordinate. When running vcf_eval, only applies to recall, not precision [infinity]",
        type=int,
        metavar="INT",
        default=float("inf"),
    )

    # ---------------------- make_truth_vcf ------------------------------------
    subparser_make_truth_vcf = subparsers.add_parser(
        "make_truth_vcf",
        parents=[common_parser],
        help="Make truth VCF file",
        usage="varifier make_truth_vcf [options] <truth_fasta> <ref_fasta> <outdir>",
        description="Make truth VCF file",
    )

    subparser_make_truth_vcf.add_argument(
        "truth_fasta", help="FASTA file of truth genome"
    )

    subparser_make_truth_vcf.add_argument(
        "ref_fasta", help="FASTA file of reference genome"
    )
    subparser_make_truth_vcf.add_argument(
        "--truth_mask",
        help="BED file of truth genome regions to mask. Any variants in the VCF matching to the mask are flagged and will not count towards precision or recall if the output VCF is used with vcf_eval",
        metavar="FILENAME",
    )
    subparser_make_truth_vcf.add_argument(
        "--cpus",
        help="Number of CPUs to use when running nucmer and minimap2 [%(default)s]",
        type=int,
        default=1,
        metavar="INT",
    )
    subparser_make_truth_vcf.add_argument(
        "--sanitise_truth_gaps",
        help="Only used if also using --global_align. Use the global alignment to change gap lengths in the truth fasta so that gaps are same length as missing sequence from the reference",
        action="store_true",
    )

    subparser_make_truth_vcf.add_argument("outdir", help="Name of output directory")
    subparser_make_truth_vcf.set_defaults(func=varifier.tasks.make_truth_vcf.run)

    # ------------------------ vcf_eval ----------------------------------------
    subparser_vcf_eval = subparsers.add_parser(
        "vcf_eval",
        parents=[common_parser],
        help="Evaluate VCF file",
        usage="varifier vcf_eval [options] <truth_fasta> <vcf_fasta> <vcf_file> <outdir>",
        description="Evaluate VCF file",
    )

    subparser_vcf_eval.add_argument(
        "--force", help="Replace outdir if it already exists", action="store_true"
    )
    subparser_vcf_eval.add_argument(
        "--filter_pass",
        help="Defines how to handle FILTER column of input VCF file. Comma-separated list of filter names. A VCF line is kept if any of its FILTER entries are in the provided list. Put '.' in the list to keep records where the filter column is '.'. Default behaviour is to ignore the filter column and use all records",
        metavar="FILTER1[,FILTER2[,...]]",
    )
    subparser_vcf_eval.add_argument(
        "--ref_mask",
        help="BED file of ref regions to mask. Any variants in the VCF overlapping the mask are removed at the start of the pipeline",
        metavar="FILENAME",
    )
    subparser_vcf_eval.add_argument(
        "--truth_mask",
        help="BED file of truth genome regions to mask. Any variants in the VCF matching to the mask are flagged and do not count towards precision or recall",
        metavar="FILENAME",
    )
    subparser_vcf_eval.add_argument(
        "--truth_vcf",
        help="VCF file of variant calls between vcf_fasta and truth_fasta, where reference of this VCF file is truth_fasta. If provided, used to calculate recall",
        metavar="FILENAME",
    )
    subparser_vcf_eval.add_argument(
        "--use_ref_calls",
        help="Include 0/0 genotype calls when calculating TPs and precision. By default they are ignored",
        action="store_true",
    )
    subparser_vcf_eval.add_argument(
        "--cpus",
        help="Number of CPUs to use when running nucmer and minimap2 to get recall, eveything else is single-core/thread. If you have a big genome, more efficient to run make_truth_vcf with >1 CPU, then use its output with --truth_vcf when running vcf_eval. [%(default)s]",
        type=int,
        default=1,
        metavar="INT",
    )
    subparser_vcf_eval.add_argument("truth_fasta", help="FASTA file of truth genome")
    subparser_vcf_eval.add_argument(
        "vcf_fasta", help="FASTA file corresponding to vcf_file"
    )
    subparser_vcf_eval.add_argument("vcf_in", help="VCF file to evaluate")
    subparser_vcf_eval.add_argument("outdir", help="Name of output directory")
    subparser_vcf_eval.set_defaults(func=varifier.tasks.vcf_eval.run)

    args = parser.parse_args()

    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s]: %(message)s", level=log_level
    )

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
