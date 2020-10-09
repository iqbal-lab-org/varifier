import copy
import logging
import os
import shutil
import subprocess

import pysam
import pysam.bcftools

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


def _truth_using_minimap2_paftools(ref_fasta, truth_fasta, vcf_file):
    _check_dependencies_in_path()
    minimap2_cmd = f"minimap2 -c --cs {ref_fasta} {truth_fasta} | sort -k6,6 -k8,8n"
    paftools_cmd = f"paftools.js call -l50 -L50 -f {ref_fasta} -"
    cmd = f"{minimap2_cmd} | {paftools_cmd} > {vcf_file}"
    logging.info(f"Running minimap2/paftools with command: {cmd}")
    subprocess.check_output(cmd, shell=True)
    logging.info(f"minimap2/paftools finished ({cmd})")


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
                new_record.FORMAT = {"GT": "1/1"}
                print(new_record, file=f)


def _filter_fps_and_long_vars_from_probe_mapped_vcf(vcf_in, vcf_out, max_ref_len):
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

    dnadiff.make_truth_vcf(ref_fasta, truth_fasta, dnadiff_vcf, debug=debug)
    _truth_using_minimap2_paftools(ref_fasta, truth_fasta, minimap2_vcf)
    to_merge = [dnadiff_vcf, minimap2_vcf]
    _merge_vcf_files_for_probe_mapping(to_merge, ref_fasta, merged_vcf)
    logging.info(f"Made merged VCF file {merged_vcf}")
    logging.info("Probe mapping to remove incorrect calls")
    probe_mapping.annotate_vcf_with_probe_mapping(
        merged_vcf,
        ref_fasta,
        truth_fasta,
        flank_length,
        probe_mapped_vcf,
        map_outfile=map_debug_file,
        truth_mask=truth_mask,
    )
    _filter_fps_and_long_vars_from_probe_mapped_vcf(
        probe_mapped_vcf, probe_filtered_vcf, max_ref_len
    )
    logging.info(f"Made filtered VCF file {probe_filtered_vcf}")
    logging.info("Using bcftools to normalise and remove duplicates")
    _bcftools_norm(ref_fasta, probe_filtered_vcf, truth_vcf)
    logging.info(f"Finished making truth VCF file {truth_vcf}")
    return truth_vcf
