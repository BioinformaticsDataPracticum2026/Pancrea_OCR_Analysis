from pathlib import Path
import textwrap
from utils import load_config, submit_job

def make_job_script(
    mapped_file,
    target_peak_file,
    output_dir,
    temp_dir,
    log_dir,
    job_name,
    source_label,
    target_label,
):
    """
    Input:
        mapped_file:
            HALPER output mapped from source species to target species.
            Example: human-to-mouse HALPER narrowPeak.

        target_peak_file:
            Native ATAC peak file from the target species.
            Example: mouse idr.optimal_peak.narrowPeak.

    Final outputs:
        shared_ocrs.bed
            Mapped source OCRs that overlap native target OCRs.

        {source_label}_specific_ocrs.bed
            Mapped source OCRs that do not overlap native target OCRs.

        mapped_vs_{target_label}.coverage.txt
            Coverage statistics comparing mapped source OCRs against native target OCRs.

        summary.txt
            Count summary and definitions.

    Temporary files:
        mapped_halper.bed
        target_atac.bed
        sorted versions

    Temporary files are placed in output/temp/{job_name}_work and removed after completion.
    """
    
    mapped_file = Path(mapped_file)
    target_peak_file = Path(target_peak_file)
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

    mapped_bed_sorted = work_dir / "mapped_halper.sorted.bed"
    target_bed_sorted = work_dir / f"{target_label}_atac.sorted.bed"

    shared_bed = output_dir / "shared_ocrs.bed"
    species_specific_bed = output_dir / f"{source_label}_specific_ocrs.bed"
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

    echo "Running BEDTools comparison"
    echo "Direction: {source_display} to {target_display}"
    echo "Mapped HALPER file: {mapped_file}"
    echo "Native target ATAC peak file: {target_peak_file}"
    echo ""

    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$5}}' {mapped_file} > {mapped_bed}
    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"{target_label}_peak_"NR}}' {target_peak_file} > {target_bed}

    sort -k1,1 -k2,2n {mapped_bed} > {mapped_bed_sorted}
    sort -k1,1 -k2,2n {target_bed} > {target_bed_sorted}

    bedtools intersect \\
      -a {mapped_bed_sorted} \\
      -b {target_bed_sorted} \\
      -u \\
      > {shared_bed}

    bedtools intersect \\
      -a {mapped_bed_sorted} \\
      -b {target_bed_sorted} \\
      -v \\
      > {species_specific_bed}

    total_mapped=$(wc -l < {mapped_bed_sorted})
    total_target=$(wc -l < {target_bed_sorted})
    total_shared=$(wc -l < {shared_bed})
    total_species_specific=$(wc -l < {species_specific_bed})

    {{
        echo "Direction: {source_display} to {target_display}"
        echo "Total mapped {source_label} HALPER regions: $total_mapped"
        echo "Total native {target_label} ATAC peaks: $total_target"
        echo "Shared OCRs: $total_shared"
        echo "{source_display}-specific OCRs: $total_species_specific"
    }} > {summary_txt}

    cat {summary_txt}

    rm -rf {work_dir}

    echo "BEDTools job finished"
    """)

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script


def main():
    config = load_config(
        path_keys=[
            "mapped_htm_file",
            "mouse_peak_file",
            "mapped_mth_file",
            "human_peak_file",
            "bed_output_dir_htm",
            "bed_output_dir_mth",
            "bed_g_temp_dir",
        ],
        required_paths=[
            "mapped_htm_file",
            "mouse_peak_file",
            "mapped_mth_file",
            "human_peak_file",
        ],
        mkdir_keys=[
            "bed_output_dir_htm",
            "bed_output_dir_mth",
            "bed_g_temp_dir",
        ],
    )

    print("BEDTools config loaded successfully.")

    htm_job = make_job_script(
        mapped_file=config["mapped_htm_file"],
        target_peak_file=config["mouse_peak_file"],
        output_dir=config["bed_output_dir_htm"],
        temp_dir=config["temp_dir"],
        log_dir=config["log_dir"],
        job_name="bed_human_pancreas_to_mouse",
        source_label="human",
        target_label="mouse",
    )

    mth_job = make_job_script(
        mapped_file=config["mapped_mth_file"],
        target_peak_file=config["human_peak_file"],
        output_dir=config["bed_output_dir_mth"],
        temp_dir=config["temp_dir"],
        log_dir=config["log_dir"],
        job_name="bed_mouse_pancreas_to_human",
        source_label="mouse",
        target_label="human",
    )

    print(f"Generated: {htm_job}")
    print(f"Generated: {mth_job}")

    submit_job(htm_job)
    submit_job(mth_job)


if __name__ == "__main__":
    main()