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

    # ---------------------- make_truth_vcf ------------------------------------
    subparser_make_truth_vcf = subparsers.add_parser(
        "make_truth_vcf",
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
        "--max_recall_ref_len",
        help="Do not include variants where REF length is more than this number. Default is no limit",
        type=int,
        metavar="INT",
    )
    subparser_make_truth_vcf.add_argument(
        "--flank_length",
        help="Length of sequence to add either side of variant when making probe sequences [%(default)s]",
        type=int,
        default=100,
        metavar="INT",
    )
    subparser_make_truth_vcf.add_argument(
        "--truth_mask",
        help="BED file of truth genome regions to mask. Any variants in the VCF matching to the mask are flagged and will not count towards precision or recall if the output VCF is used with vcf_eval",
        metavar="FILENAME",
    )

    subparser_make_truth_vcf.add_argument("outdir", help="Name of output directory")
    subparser_make_truth_vcf.set_defaults(func=varifier.tasks.make_truth_vcf.run)

    # ------------------------ vcf_eval ----------------------------------------
    subparser_vcf_eval = subparsers.add_parser(
        "vcf_eval",
        help="Evaluate VCF file",
        usage="varifier vcf_eval [options] <truth_fasta> <vcf_fasta> <vcf_file> <outdir>",
        description="Evaluate VCF file",
    )

    subparser_vcf_eval.add_argument(
        "--flank_length",
        help="Length of sequence to add either side of variant when making probe sequences [%(default)s]",
        type=int,
        default=100,
        metavar="INT",
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
        "--max_recall_ref_len",
        help="For recall, do not look for expected variants where REF length is more than this number. Default is no limit. This option will not work if you use --truth_vcf",
        type=int,
        metavar="INT",
    )
    subparser_vcf_eval.add_argument(
        "--use_ref_calls",
        help="Include 0/0 genotype calls when calculating TPs and precision. By default they are ignored",
        action="store_true",
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
