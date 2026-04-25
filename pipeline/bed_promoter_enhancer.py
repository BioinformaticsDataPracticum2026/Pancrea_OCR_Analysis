from pathlib import Path
import textwrap
from utils import load_config, submit_job

def make_job_script(
    mapped_file,
    target_peak_file,
    target_tss_file,
    output_dir,
    temp_dir,
    log_dir,
    job_name,
    source_label,
    target_label,
    promoter_threshold,
):
    mapped_file = Path(mapped_file)
    target_peak_file = Path(target_peak_file)
    target_tss_file = Path(target_tss_file)
    output_dir = Path(output_dir)
    temp_dir = Path(temp_dir)
    log_dir = Path(log_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
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

    mapped_promoter = work_dir / "mapped_promoter_ocrs.bed"
    mapped_enhancer = work_dir / "mapped_enhancer_ocrs.bed"
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

    mkdir -p {work_dir}
    mkdir -p {output_dir}

    echo "Running promoter/enhancer BEDTools analysis"
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

    echo "Promoter/enhancer BEDTools analysis finished"
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
        "bed_pe_temp_dir",
    ],
    )
    print("BEDTools config loaded successfully.")

    htm_job = make_job_script(
        mapped_file=config["mapped_htm_file"],
        target_peak_file=config["mouse_peak_file"],
        target_tss_file=config["mouse_tss_file"],
        output_dir=config["bed_pe_output_dir_htm"],
        temp_dir=config["bed_pe_temp_dir"],
        log_dir=config["log_dir"],
        job_name="bed_pe_human_pancreas_to_mouse",
        source_label="human",
        target_label="mouse",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    mth_job = make_job_script(
        mapped_file=config["mapped_mth_file"],
        target_peak_file=config["human_peak_file"],
        target_tss_file=config["human_tss_file"],
        output_dir=config["bed_pe_output_dir_mth"],
        temp_dir=config["bed_pe_temp_dir"],
        log_dir=config["log_dir"],
        job_name="bed_pe_mouse_pancreas_to_human",
        source_label="mouse",
        target_label="human",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    print(f"Generated: {htm_job}")
    print(f"Generated: {mth_job}")

    submit_job(htm_job)
    submit_job(mth_job)


if __name__ == "__main__":
    main()