from pathlib import Path
import subprocess
import textwrap
import yaml

# run with:
# python3 /ocean/projects/bio230007p/jji5/pipeline/bed.py

def load_config():
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config.yaml"

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    required_paths = [
        ("mapped_htm_file", Path(config["mapped_htm_file"])),
        ("mouse_peak_file", Path(config["mouse_peak_file"])),
        ("mapped_mth_file", Path(config["mapped_mth_file"])),
        ("human_peak_file", Path(config["human_peak_file"])),
    ]

    # Check if files exist
    for label, path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    Path(config["bed_output_dir_htm"]).mkdir(parents=True, exist_ok=True)
    Path(config["bed_output_dir_mth"]).mkdir(parents=True, exist_ok=True)
    Path(config["temp_dir"]).mkdir(parents=True, exist_ok=True)

    return config

# generates job bash job scripts
def make_job_script(mapped_file, target_peak_file, output_dir, temp_dir, job_name, target_label):
    """
    Input: mapped source_to_target.halper.narrowPeak file, target peak file (original)
    Output: lots of intermediate output, important: shared_ocrs.bed 
    bedtools intersect find shared OCR. 
    """
    mapped_file = Path(mapped_file)
    target_peak_file = Path(target_peak_file)
    output_dir = Path(output_dir)
    temp_dir = Path(temp_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    job_script = temp_dir / f"{job_name}.job"
    out_log = output_dir / f"{job_name}.out.txt"
    err_log = output_dir / f"{job_name}.err.txt"

    mapped_bed = output_dir / "mapped_halper.bed"
    target_bed = output_dir / f"{target_label}_atac.bed"
    shared_bed = output_dir / "shared_ocrs.bed"
    mapped_only_bed = output_dir / f"mapped_not_open_in_{target_label}.bed"
    coverage_txt = output_dir / f"mapped_vs_{target_label}.coverage.txt"

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

    echo "Running bedtools comparison"
    echo "Mapped HALPER file: {mapped_file}"
    echo "Target peak file: {target_peak_file}"

    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$5}}' {mapped_file} > {mapped_bed}
    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"{target_label}_peak_"NR}}' {target_peak_file} > {target_bed}

    sort -k1,1 -k2,2n {mapped_bed} > {mapped_bed}.sorted
    sort -k1,1 -k2,2n {target_bed} > {target_bed}.sorted

    bedtools intersect -a {mapped_bed}.sorted -b {target_bed}.sorted -u > {shared_bed}
    bedtools intersect -a {mapped_bed}.sorted -b {target_bed}.sorted -v > {mapped_only_bed}
    bedtools coverage -a {mapped_bed}.sorted -b {target_bed}.sorted > {coverage_txt}

    echo "Mapped HALPER regions:" $(wc -l < {mapped_bed}.sorted)
    echo "Shared OCRs:" $(wc -l < {shared_bed})
    echo "Mapped but not open in {target_label}:" $(wc -l < {mapped_only_bed})
    echo "Bedtools job finished"
    """)

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script


def submit_job(job_script):
    result = subprocess.run(
        ["sbatch", str(job_script)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    print(result.stdout)
    if result.stderr:
        print(result.stderr)


def main():
    config = load_config()
    print("Bedtools config loaded successfully.")
    # Human to Mouse bed.
    htm_job = make_job_script(
        mapped_file=config["mapped_htm_file"],
        target_peak_file=config["mouse_peak_file"],
        output_dir=config["bed_output_dir_htm"],
        temp_dir=config["temp_dir"],
        job_name="bed_human_pancreas_to_mouse",
        target_label="mouse",
    )
    # Mouse to Human bed.
    mth_job = make_job_script(
        mapped_file=config["mapped_mth_file"],
        target_peak_file=config["human_peak_file"],
        output_dir=config["bed_output_dir_mth"],
        temp_dir=config["temp_dir"],
        job_name="bed_mouse_pancreas_to_human",
        target_label="human",
    )

    print(f"Generated: {htm_job}")
    print(f"Generated: {mth_job}")

    submit_job(htm_job)
    submit_job(mth_job)


if __name__ == "__main__":
    main()