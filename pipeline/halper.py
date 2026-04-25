from pathlib import Path
import textwrap
from utils import load_config, submit_job

# generates job bash job scripts
def make_job_script(config, source_species, target_species, peak_file, job_name, output_dir):
    """
    Final useful outputs:
        1. *.HALPER.narrowPeak
        2. *.HALPER.narrowPeak.failed
        3. *.HALPER.narrowPeak.png
        4. *.HALPER.narrowPeak-peak.png, if generated

    Temporary files used internally:
        - converted peak BED
        - summit BED
        - lifted peak BED
        - lifted summit BED
    These are deleted after the run.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    temp_dir = Path(config["halper_temp_dir"])
    temp_dir.mkdir(parents=True, exist_ok=True)

    hal_file = Path(config["hal_file"])
    halper_repo = Path(config["halper_repo"])
    ortholog_script = halper_repo / "orthologFind.py"

    peak_file = Path(peak_file)
    peak_base = peak_file.stem

    job_script = temp_dir / f"{job_name}.job"

    # Keep logs outside the final HALPER output folder
    log_dir = temp_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    out_log = log_dir / f"{job_name}.out.txt"
    err_log = log_dir / f"{job_name}.err.txt"

    # Temporary working directory for files needed only during execution
    work_dir = temp_dir / f"{job_name}_work"

    # Final HALPER output file
    output_file = output_dir / f"{peak_base}.{source_species}To{target_species}.HALPER.narrowPeak"

    # Temporary intermediate files
    peak_bed = work_dir / f"{peak_base}.{source_species}To{target_species}.bed"
    summit_bed = work_dir / f"{peak_base}.{source_species}To{target_species}.summits.bed"
    lifted_peak_bed = work_dir / f"{peak_base}.{source_species}To{target_species}.halLiftover.tFile.bed"
    lifted_summit_bed = work_dir / f"{peak_base}.{source_species}To{target_species}.halLiftover.sFile.bed"

    script_text = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -p RM-shared
    #SBATCH -t 12:00:00
    #SBATCH --ntasks-per-node=4
    #SBATCH --job-name={job_name}
    #SBATCH --output={out_log}
    #SBATCH --error={err_log}

    set -euo pipefail

    module load anaconda3
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate {config["conda_env"]}

    export PATH=/ocean/projects/bio230007p/jji5/tools/hal/bin:$PATH

    mkdir -p {work_dir}
    mkdir -p {output_dir}

    echo "Running HALPER: {source_species} -> {target_species}"
    echo "Peak file: {peak_file}"
    echo "HAL file: {hal_file}"
    echo "Final output file: {output_file}"

    which halLiftover
    which python

    echo "Creating temporary 4-column peak BED"
    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"peak_"NR}}' {peak_file} > {peak_bed}

    echo "Creating temporary summit BED"
    awk 'BEGIN{{OFS="\\t"}} {{summit=$2+$10; print $1,summit,summit+1,"peak_"NR}}' {peak_file} > {summit_bed}

    echo "Running halLiftover on peaks"
    halLiftover {hal_file} {source_species} {peak_bed} {target_species} {lifted_peak_bed}

    echo "Running halLiftover on summits"
    halLiftover {hal_file} {source_species} {summit_bed} {target_species} {lifted_summit_bed}

    echo "Running orthologFind.py"
    python {ortholog_script} \\
      -max_len 1000 \\
      -min_len 50 \\
      -protect_dist 5 \\
      -mult_keepone \\
      -narrowPeak \\
      -qFile {peak_bed} \\
      -tFile {lifted_peak_bed} \\
      -sFile {lifted_summit_bed} \\
      -oFile {output_file}

    if [ ! -s {output_file} ]; then
        echo "ERROR: Output file missing or empty: {output_file}" >&2
        exit 1
    fi

    echo "Cleaning temporary intermediate files"
    rm -rf {work_dir}

    echo "Finished HALPER: {source_species} -> {target_species}"
    echo "Final files should be:"
    echo "{output_file}"
    echo "{output_file}.failed"
    echo "{output_file}.png"
    echo "{str(output_file).replace('.narrowPeak', '.narrowPeak-peak.png')}"
    """)

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script


def main():
    config = load_config(
    path_keys=[
        "halper_repo",
        "hal_file",
        "species_1_peak_file",
        "species_2_peak_file",
        "halper_output_dir_htm",
        "halper_output_dir_mth",
        "halper_temp_dir",
        "conda_env",
    ],
    required_paths=[
        "halper_repo",
        "hal_file",
        "species_1_peak_file",
        "species_2_peak_file",
        "conda_env",
    ],
    mkdir_keys=[
        "halper_output_dir_htm",
        "halper_output_dir_mth",
        "halper_temp_dir",
    ],
    )
    print("Config loaded successfully.")

    human_to_mouse_job = make_job_script(
        config=config,
        source_species=config["species_1"],
        target_species=config["species_2"],
        peak_file=config["species_1_peak_file"],
        job_name="halper_human_pancreas_to_mouse",
        output_dir=config["halper_output_dir_htm"],
    )

    mouse_to_human_job = make_job_script(
        config=config,
        source_species=config["species_2"],
        target_species=config["species_1"],
        peak_file=config["species_2_peak_file"],
        job_name="halper_mouse_pancreas_to_human",
        output_dir=config["halper_output_dir_mth"],
    )

    print(f"Generated: {human_to_mouse_job}")
    print(f"Generated: {mouse_to_human_job}")

    submit_job(human_to_mouse_job)
    submit_job(mouse_to_human_job)


if __name__ == "__main__":
    main()