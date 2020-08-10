import copy
import logging
import os
import shutil
import subprocess

import pysam
import pysam.bcftools
import re

from cluster_vcf_records import vcf_file_read, vcf_merge

from varifier import dnadiff, probe_mapping, utils


def _check_dependencies_in_path():
    programs = ["minimap2", "paftools.js", "k8"]
    found_all = True
    for program in programs:
        if shutil.which(program) is None:
            found_all = False
            logging.error(f"Not found in PATH: {program}")

    if not found_all:
        raise RuntimeError("Some required programs not found in $PATH. Cannot continue")



qname_matcher = re.compile(".*QNAME=(.+?);")
def fix_minimap2_vcf(input_vcf_file, output_vcf_file, ref_fasta, qry_fasta, snps_only):
    ref_seqs = utils.file_to_dict_of_seqs(ref_fasta)
    qry_seqs = utils.file_to_dict_of_seqs(qry_fasta)

    discarded_variants = []
    with open(input_vcf_file) as input_vcf_filehandler,\
         open(output_vcf_file, "w") as output_vcf_filehandler:
        for line in input_vcf_filehandler:
            if line.startswith("#"):
                output_vcf_filehandler.write(line)
            else:
                line_split = line.strip().split()

                # change QUAL to "." so that it is equal to dnadiff
                line_split[5] = "."

                # get REF_LENGTH and QRY_LENGTH
                ref_chrom = line_split[0]
                ref_length = len(ref_seqs[ref_chrom])
                info = line_split[7]
                qry_chrom = qname_matcher.match(info).group(1)
                qry_length = len(qry_seqs[qry_chrom])

                # change info
                info = f"LENGTH_QRY={qry_length};LENGTH_REF={ref_length};" + info
                line_split[7] = info

                # remake line
                line = "\t".join(line_split)

                ref = line_split[3]
                alts = line_split[4]
                is_snp = ref in ["A", "C", "G", "T"] and alts in ["A", "C", "G", "T"]

                if not snps_only or (snps_only and is_snp):
                    print(line, file=output_vcf_filehandler)
                elif snps_only and not is_snp:
                    discarded_variants.append(line)

    if snps_only:
        with open(f"{output_vcf_file}.discarded_not_snps.vcf", "w") as discarded_variants_fh:
            for discarded_variant in discarded_variants:
                print(discarded_variant, file=discarded_variants_fh)




def _truth_using_minimap2_paftools(ref_fasta, truth_fasta, vcf_file, snps_only):
    _check_dependencies_in_path()
    minimap2_cmd = f"minimap2 -c --cs {ref_fasta} {truth_fasta} | sort -k6,6 -k8,8n"

    script_dir = utils.get_script_dir()
    paftools_cmd = f"k8 {script_dir}/paftools_fixed.js call -l50 -L50 -f {ref_fasta} -"
    cmd = f"{minimap2_cmd} | {paftools_cmd} > {vcf_file}.paftools_raw_output"
    logging.info(f"Running minimap2/paftools with command: {cmd}")
    subprocess.check_output(cmd, shell=True)
    logging.info(f"minimap2/paftools finished ({cmd})")
    fix_minimap2_vcf(f"{vcf_file}.paftools_raw_output", vcf_file, ref_fasta, truth_fasta, snps_only=snps_only)
    os.unlink(f"{vcf_file}.paftools_raw_output")


def _deduplicate_vcf_files(to_merge, disagreement_file):
    ref_chrom_ref_pos_to_vcf_line = {}

    with open(disagreement_file, "w") as disagreement_filehandler:
        for file in to_merge:
            with open(file) as file_handler:
                for line in file_handler:
                    line = line.strip()
                    is_header = line.startswith("#")
                    is_empty = len(line) == 0
                    if not is_header and not is_empty:
                        line_split = line.split("\t")

                        # change FILTER
                        line_split[6] = "PASS"
                        # change FORMAT
                        line_split[8] = "GT:VFR_FILTER"
                        # change sample data
                        line_split[9] = "1/1:PASS"

                        # get some other info
                        ref_chrom, ref_pos = line_split[0], int(line_split[1])

                        # update line
                        line = "\t".join(line_split)

                        this_is_a_variant_present_in_both_VCFs = (ref_chrom, ref_pos) in ref_chrom_ref_pos_to_vcf_line
                        if this_is_a_variant_present_in_both_VCFs:
                            both_VCF_lines_are_identical = ref_chrom_ref_pos_to_vcf_line[(ref_chrom, ref_pos)] == line

                            if not both_VCF_lines_are_identical:
                                # here, dnadiff and minimap2 disagrees on a record, let's remove it for sanity reasons
                                # log it first
                                print(ref_chrom_ref_pos_to_vcf_line[(ref_chrom, ref_pos)], file=disagreement_filehandler)
                                print(line, file=disagreement_filehandler)

                                # and remove it
                                del ref_chrom_ref_pos_to_vcf_line[(ref_chrom, ref_pos)]
                        else:
                            ref_chrom_ref_pos_to_vcf_line[(ref_chrom, ref_pos)] = line

    ordered_ref_chrom_ref_pos = sorted(list(ref_chrom_ref_pos_to_vcf_line.keys()))
    ordered_vcf_lines = [ref_chrom_ref_pos_to_vcf_line[ref_chrom_ref_pos] for ref_chrom_ref_pos in ordered_ref_chrom_ref_pos]
    return ordered_vcf_lines


def _identify_vcf_lines(vcf_lines):
    identified_vcf_lines = []
    for index, vcf_line in enumerate(vcf_lines):
        vcf_line_split = vcf_line.strip().split("\t")
        vcf_line_split[2] = str(index)
        identified_vcf_lines.append("\t".join(vcf_line_split))
    return identified_vcf_lines

def _deduplicate_vcf_files_for_probe_mapping(to_merge, ref_fasta, vcf_out):
    ref_seqs = utils.file_to_dict_of_seqs(ref_fasta)
    vcf_lines = _deduplicate_vcf_files(to_merge, f"{vcf_out}.disagreements_between_dnadiff_and_minimap2")
    vcf_lines = _identify_vcf_lines(vcf_lines)

    with open(vcf_out, "w") as f:
        print("##fileformat=VCFv4.2", file=f)
        for seq in ref_seqs.values():
            print(f"##contig=<ID={seq.id},length={len(seq)}>", file=f)
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=f)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample", file=f)
        print("\n".join(vcf_lines), file=f)

def _merge_vcf_files_for_probe_mapping(list_of_vcf_files, ref_fasta, vcf_out):
    ref_seqs = utils.file_to_dict_of_seqs(ref_fasta)
    # This makes a merged file, where two different ALTs at the same place
    # result in one record with a list of ALTs. For probe mapping, we want
    # a separate record for each allele. Also need genotype to be "1/1"
    vcf_merge.merge_vcf_files(list_of_vcf_files, ref_seqs, vcf_out)
    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(vcf_out)
    with open(vcf_out, "w") as f:
        print("##fileformat=VCFv4.2", file=f)
        for seq in ref_seqs.values():
            print(f"##contig=<ID={seq.id},length={len(seq)}>", file=f)
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', file=f)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample", file=f)

        for i, record in enumerate(vcf_records):
            for alt in record.ALT:
                new_record = copy.copy(record)
                new_record.ID = str(i)
                new_record.ALT = [alt]
                new_record.INFO = {}
                new_record.FILTER = set(["PASS"])
                new_record.FORMAT = {"GT": "1/1", "VFR_FILTER": "PASS"}
                print(new_record, file=f)


def _filter_fps_and_long_vars_from_probe_mapped_vcf(vcf_in, vcf_out, max_ref_len, detailed_VCF=False):
    """vcf_in should be file made by _merge_vcf_files_for_probe_mapping, and
    then annotated using probe_mapping.annotate_vcf_with_probe_mapping().
    Outputs a new VCF file that only contains the TPs, based on probe mapping"""
    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(vcf_in)
    with open(vcf_out, "w") as f:
        for line in header_lines:
            if (
                line.startswith("#CHROM")
                or line.startswith("##fileformat")
                or line.startswith("##contig")
                or line.startswith("##FORMAT=<ID=GT")
            ):
                print(line, file=f)

        i = 0
        while i < len(vcf_records):
            records = [vcf_records[i]]
            i += 1
            while i < len(vcf_records) and vcf_records[i].ID == records[0].ID:
                records.append(vcf_records[i])
                i += 1

            records = [
                x
                for x in records
                if x.FORMAT.get("VFR_RESULT", "FP") == "TP"
                and x.FORMAT.get("VFR_IN_MASK", "0") != "1"
            ]
            if max_ref_len is not None:
                records = [x for x in records if len(x.REF) <= max_ref_len]

            if len(records) > 1:
                logging.warning(
                    "Skipping the following VCF lines. They conflict, but probe mapping thinks they are all true-positives:"
                )
                for record in records:
                    logging.warning(f"  {record}")
            elif len(records) == 1:
                if not detailed_VCF:
                    records[0].FORMAT = {"GT": "1/1"}
                print(records[0], file=f)


def _bcftools_norm(ref_fasta, vcf_in, vcf_out):
    """Runs bcftools norm, to normalise variants and remove
    duplicates. Experimenting has shown that you have to run it twice,
    so that duplicates get removed properly."""
    # The "-o" option doesn't get passed to pysam's bcftools wrapper.
    # Instead it returns a string that is the new vcf file
    options = ["-c", "x", "-d", "any", "-f", ref_fasta]
    with open(vcf_out, "w") as f:
        print(pysam.bcftools.norm(*options, vcf_in), end="", file=f)
    vcf_string = pysam.bcftools.norm(*options, vcf_out)
    with open(vcf_out, "w") as f:
        print(vcf_string, end="", file=f)


def make_truth_vcf(
    ref_fasta,
    truth_fasta,
    outdir,
    flank_length,
    debug=False,
    truth_mask=None,
    max_ref_len=None,
    snps_only=False,
    output_probes=False,
    detailed_VCF=False
):
    _check_dependencies_in_path()
    os.mkdir(outdir)
    minimap2_vcf = os.path.join(outdir, "00.minimap2.vcf")
    dnadiff_vcf = os.path.join(outdir, "00.dnadiff.vcf")
    merged_vcf = os.path.join(outdir, "01.merged.vcf")
    if debug:
        map_debug_file = os.path.join(outdir, "02.debug.map")
    else:
        map_debug_file = None
    probe_mapped_vcf = os.path.join(outdir, "02.merged_and_probe_mapped.vcf")
    probe_filtered_vcf = os.path.join(outdir, "03.probe_filtered.vcf")
    truth_vcf = os.path.join(outdir, "04.truth.vcf")

    dnadiff.make_truth_vcf(ref_fasta, truth_fasta, dnadiff_vcf, snps_only=snps_only, debug=debug)
    _truth_using_minimap2_paftools(ref_fasta, truth_fasta, minimap2_vcf, snps_only=snps_only)
    to_merge = [dnadiff_vcf, minimap2_vcf]

    if snps_only:
        _deduplicate_vcf_files_for_probe_mapping(to_merge, ref_fasta, merged_vcf)
    else:
        _merge_vcf_files_for_probe_mapping(to_merge, ref_fasta, merged_vcf)
    logging.info(f"Made merged VCF file {merged_vcf}")
    logging.info(f"Probe mapping to remove incorrect calls")
    probe_mapping.annotate_vcf_with_probe_mapping(
        merged_vcf,
        ref_fasta,
        truth_fasta,
        flank_length,
        probe_mapped_vcf,
        map_outfile=map_debug_file,
        truth_mask=truth_mask,
        output_probes=output_probes
    )
    _filter_fps_and_long_vars_from_probe_mapped_vcf(
        probe_mapped_vcf, probe_filtered_vcf, max_ref_len, detailed_VCF=detailed_VCF
    )
    logging.info(f"Made filtered VCF file {probe_filtered_vcf}")
    logging.info(f"Using bcftools to normalise and remove duplicates")
    _bcftools_norm(ref_fasta, probe_filtered_vcf, truth_vcf)
    logging.info(f"Finished making truth VCF file {truth_vcf}")
    return truth_vcf
