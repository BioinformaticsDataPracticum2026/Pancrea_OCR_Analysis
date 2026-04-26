from pathlib import Path
import textwrap

from utils import load_config, submit_job


def make_species_pe_job(
    peak_file,
    tss_file,
    output_dir,
    temp_dir,
    job_name,
    species_label,
    promoter_threshold,
):
    """
    Classify native species OCRs into promoter-like and enhancer-like OCRs.

    Final outputs:
        {species_label}_promoter_ocrs.bed
        {species_label}_enhancer_ocrs.bed
        summary.txt
    """

    peak_file = Path(peak_file)
    tss_file = Path(tss_file)
    output_dir = Path(output_dir)
    temp_dir = Path(temp_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    log_dir = temp_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    species_display = species_label.capitalize()

    job_script = temp_dir / f"{job_name}.job"
    out_log = log_dir / f"{job_name}.out.txt"
    err_log = log_dir / f"{job_name}.err.txt"

    work_dir = temp_dir / f"{job_name}_work"

    ocr_bed = work_dir / f"{species_label}_ocrs.bed"
    ocr_bed_sorted = work_dir / f"{species_label}_ocrs.sorted.bed"
    tss_sorted = work_dir / f"{species_label}_tss.sorted.bed"
    closest_file = work_dir / f"{species_label}_with_nearest_tss.bed"

    promoter_file = output_dir / f"{species_label}_promoter_ocrs.bed"
    enhancer_file = output_dir / f"{species_label}_enhancer_ocrs.bed"
    summary_txt = output_dir / "summary.txt"

    script_text = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -p RM-shared
    #SBATCH -t 04:00:00
    #SBATCH --ntasks-per-node=4
    #SBATCH --job-name={job_name}
    #SBATCH --output={out_log}
    #SBATCH --error={err_log}

    set -euo pipefail

    module load bedtools/2.30.0

    rm -rf {work_dir}
    mkdir -p {work_dir}
    mkdir -p {output_dir}

    echo "Running full species promoter/enhancer classification"
    echo "Species: {species_display}"
    echo "Peak file: {peak_file}"
    echo "TSS file: {tss_file}"
    echo "Promoter threshold: {promoter_threshold} bp"
    echo ""

    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"{species_label}_peak_"NR}}' {peak_file} > {ocr_bed}

    sort -k1,1 -k2,2n {ocr_bed} > {ocr_bed_sorted}
    sort -k1,1 -k2,2n {tss_file} > {tss_sorted}

    bedtools closest \\
      -a {ocr_bed_sorted} \\
      -b {tss_sorted} \\
      -d \\
      > {closest_file}

    awk 'BEGIN{{OFS="\\t"}} $NF <= {promoter_threshold} {{print $1,$2,$3,$4}}' {closest_file} > {promoter_file}
    awk 'BEGIN{{OFS="\\t"}} $NF >  {promoter_threshold} {{print $1,$2,$3,$4}}' {closest_file} > {enhancer_file}

    total_ocrs=$(wc -l < {ocr_bed_sorted})
    promoter_n=$(wc -l < {promoter_file})
    enhancer_n=$(wc -l < {enhancer_file})

    {{
        echo "Species: {species_display}"
        echo "Promoter threshold bp: {promoter_threshold}"
        echo "Total {species_label} OCRs: $total_ocrs"
        echo "{species_display} promoter OCRs: $promoter_n"
        echo "{species_display} enhancer OCRs: $enhancer_n"
    }} > {summary_txt}

    cat {summary_txt}

    rm -rf {work_dir}

    echo "Full species promoter/enhancer classification finished"
    """)

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script


def make_mapped_pe_job(
    mapped_file,
    target_peak_file,
    target_tss_file,
    output_dir,
    temp_dir,
    job_name,
    source_label,
    target_label,
    promoter_threshold,
):
    """
    Classify mapped source OCRs and native target OCRs into promoter/enhancer sets,
    then identify shared and source-specific promoter/enhancer OCRs.

    Final outputs:
        shared_promoter_ocrs.bed
        shared_enhancer_ocrs.bed
        {source_label}_specific_promoter_ocrs.bed
        {source_label}_specific_enhancer_ocrs.bed
        summary.txt
    """

    mapped_file = Path(mapped_file)
    target_peak_file = Path(target_peak_file)
    target_tss_file = Path(target_tss_file)
    output_dir = Path(output_dir)
    temp_dir = Path(temp_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    log_dir = temp_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    source_display = source_label.capitalize()
    target_display = target_label.capitalize()

    job_script = temp_dir / f"{job_name}.job"
    out_log = log_dir / f"{job_name}.out.txt"
    err_log = log_dir / f"{job_name}.err.txt"

    work_dir = temp_dir / f"{job_name}_work"

    mapped_bed = work_dir / "mapped_halper.bed"
    target_bed = work_dir / f"{target_label}_atac.bed"
    target_tss_sorted = work_dir / f"{target_label}_tss.sorted.bed"

    mapped_bed_sorted = work_dir / "mapped_halper.sorted.bed"
    target_bed_sorted = work_dir / f"{target_label}_atac.sorted.bed"

    mapped_closest = work_dir / "mapped_with_nearest_tss.bed"
    target_closest = work_dir / f"{target_label}_with_nearest_tss.bed"

    mapped_promoter = work_dir / f"mapped_{source_label}_promoter_ocrs.bed"
    mapped_enhancer = work_dir / f"mapped_{source_label}_enhancer_ocrs.bed"

    target_promoter = work_dir / f"{target_label}_promoter_ocrs.bed"
    target_enhancer = work_dir / f"{target_label}_enhancer_ocrs.bed"

    shared_promoter = output_dir / "shared_promoter_ocrs.bed"
    shared_enhancer = output_dir / "shared_enhancer_ocrs.bed"

    source_specific_promoter = output_dir / f"{source_label}_specific_promoter_ocrs.bed"
    source_specific_enhancer = output_dir / f"{source_label}_specific_enhancer_ocrs.bed"

    summary_txt = output_dir / "summary.txt"

    script_text = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -p RM-shared
    #SBATCH -t 04:00:00
    #SBATCH --ntasks-per-node=4
    #SBATCH --job-name={job_name}
    #SBATCH --output={out_log}
    #SBATCH --error={err_log}

    set -euo pipefail

    module load bedtools/2.30.0

    rm -rf {work_dir}
    mkdir -p {work_dir}
    mkdir -p {output_dir}

    echo "Running mapped promoter/enhancer BEDTools analysis"
    echo "Direction: {source_display} to {target_display}"
    echo "Mapped HALPER file: {mapped_file}"
    echo "Native target ATAC peak file: {target_peak_file}"
    echo "Target TSS file: {target_tss_file}"
    echo "Promoter threshold: {promoter_threshold} bp"
    echo ""

    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$5}}' {mapped_file} > {mapped_bed}
    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"{target_label}_peak_"NR}}' {target_peak_file} > {target_bed}

    sort -k1,1 -k2,2n {mapped_bed} > {mapped_bed_sorted}
    sort -k1,1 -k2,2n {target_bed} > {target_bed_sorted}
    sort -k1,1 -k2,2n {target_tss_file} > {target_tss_sorted}

    bedtools closest \\
      -a {mapped_bed_sorted} \\
      -b {target_tss_sorted} \\
      -d \\
      > {mapped_closest}

    bedtools closest \\
      -a {target_bed_sorted} \\
      -b {target_tss_sorted} \\
      -d \\
      > {target_closest}

    awk 'BEGIN{{OFS="\\t"}} $NF <= {promoter_threshold} {{print $1,$2,$3,$4}}' {mapped_closest} > {mapped_promoter}
    awk 'BEGIN{{OFS="\\t"}} $NF >  {promoter_threshold} {{print $1,$2,$3,$4}}' {mapped_closest} > {mapped_enhancer}

    awk 'BEGIN{{OFS="\\t"}} $NF <= {promoter_threshold} {{print $1,$2,$3,$4}}' {target_closest} > {target_promoter}
    awk 'BEGIN{{OFS="\\t"}} $NF >  {promoter_threshold} {{print $1,$2,$3,$4}}' {target_closest} > {target_enhancer}

    bedtools intersect \\
      -a {mapped_promoter} \\
      -b {target_promoter} \\
      -u \\
      > {shared_promoter}

    bedtools intersect \\
      -a {mapped_enhancer} \\
      -b {target_enhancer} \\
      -u \\
      > {shared_enhancer}

    bedtools intersect \\
      -a {mapped_promoter} \\
      -b {target_promoter} \\
      -v \\
      > {source_specific_promoter}

    bedtools intersect \\
      -a {mapped_enhancer} \\
      -b {target_enhancer} \\
      -v \\
      > {source_specific_enhancer}

    mapped_promoter_n=$(wc -l < {mapped_promoter})
    mapped_enhancer_n=$(wc -l < {mapped_enhancer})
    target_promoter_n=$(wc -l < {target_promoter})
    target_enhancer_n=$(wc -l < {target_enhancer})
    shared_promoter_n=$(wc -l < {shared_promoter})
    shared_enhancer_n=$(wc -l < {shared_enhancer})
    source_specific_promoter_n=$(wc -l < {source_specific_promoter})
    source_specific_enhancer_n=$(wc -l < {source_specific_enhancer})

    {{
        echo "Direction: {source_display} to {target_display}"
        echo "Promoter threshold bp: {promoter_threshold}"
        echo "Mapped {source_label} promoter OCRs: $mapped_promoter_n"
        echo "Mapped {source_label} enhancer OCRs: $mapped_enhancer_n"
        echo "Native {target_label} promoter OCRs: $target_promoter_n"
        echo "Native {target_label} enhancer OCRs: $target_enhancer_n"
        echo "Shared promoter OCRs: $shared_promoter_n"
        echo "Shared enhancer OCRs: $shared_enhancer_n"
        echo "{source_display}-specific promoter OCRs: $source_specific_promoter_n"
        echo "{source_display}-specific enhancer OCRs: $source_specific_enhancer_n"
    }} > {summary_txt}

    cat {summary_txt}

    rm -rf {work_dir}

    echo "Mapped promoter/enhancer BEDTools analysis finished"
    """)

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script


def main():
    config = load_config(
        path_keys=[
            "mapped_htm_file",
            "mapped_mth_file",
            "mouse_peak_file",
            "human_peak_file",
            "mouse_tss_file",
            "human_tss_file",
            "bed_pe_output_dir_htm",
            "bed_pe_output_dir_mth",
            "bed_pe_output_dir_human",
            "bed_pe_output_dir_mouse",
            "bed_pe_temp_dir",
        ],
        required_paths=[
            "mapped_htm_file",
            "mapped_mth_file",
            "mouse_peak_file",
            "human_peak_file",
            "mouse_tss_file",
            "human_tss_file",
        ],
        mkdir_keys=[
            "bed_pe_output_dir_htm",
            "bed_pe_output_dir_mth",
            "bed_pe_output_dir_human",
            "bed_pe_output_dir_mouse",
            "bed_pe_temp_dir",
        ],
    )

    print("Promoter/enhancer BEDTools config loaded successfully.")

    human_full_job = make_species_pe_job(
        peak_file=config["human_peak_file"],
        tss_file=config["human_tss_file"],
        output_dir=config["bed_pe_output_dir_human"],
        temp_dir=config["bed_pe_temp_dir"],
        job_name="bed_pe_human_full",
        species_label="human",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    mouse_full_job = make_species_pe_job(
        peak_file=config["mouse_peak_file"],
        tss_file=config["mouse_tss_file"],
        output_dir=config["bed_pe_output_dir_mouse"],
        temp_dir=config["bed_pe_temp_dir"],
        job_name="bed_pe_mouse_full",
        species_label="mouse",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    htm_job = make_mapped_pe_job(
        mapped_file=config["mapped_htm_file"],
        target_peak_file=config["mouse_peak_file"],
        target_tss_file=config["mouse_tss_file"],
        output_dir=config["bed_pe_output_dir_htm"],
        temp_dir=config["bed_pe_temp_dir"],
        job_name="bed_pe_human_pancreas_to_mouse",
        source_label="human",
        target_label="mouse",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    mth_job = make_mapped_pe_job(
        mapped_file=config["mapped_mth_file"],
        target_peak_file=config["human_peak_file"],
        target_tss_file=config["human_tss_file"],
        output_dir=config["bed_pe_output_dir_mth"],
        temp_dir=config["bed_pe_temp_dir"],
        job_name="bed_pe_mouse_pancreas_to_human",
        source_label="mouse",
        target_label="human",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    jobs = [
        human_full_job,
        mouse_full_job,
        htm_job,
        mth_job,
    ]

    for job in jobs:
        print(f"Generated: {job}")

    for job in jobs:
        submit_job(job)


if __name__ == "__main__":
    main()