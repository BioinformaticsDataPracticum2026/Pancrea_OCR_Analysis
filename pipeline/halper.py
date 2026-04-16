from pathlib import Path
import subprocess
import textwrap
import yaml

# run with:
# python3 /ocean/projects/bio230007p/jji5/pipeline/halper.py

def load_config():
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config.yaml"

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    required_paths = [
        ("halper_repo", Path(config["halper_repo"])),
        ("hal_file", Path(config["hal_file"])),
        ("species_1_peak_file", Path(config["species_1_peak_file"])),
        ("species_2_peak_file", Path(config["species_2_peak_file"])),
    ]
    # Check if files exist
    for label, path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")
    # Check if file format is correct (must be .narrowPeak)
    for peak_key in ["species_1_peak_file", "species_2_peak_file"]:
        peak_path = Path(config[peak_key])
        if peak_path.suffix != ".narrowPeak":
            raise ValueError(f"{peak_key} must be a .narrowPeak file: {peak_path}")

    Path(config["halper_output_dir_htm"]).mkdir(parents=True, exist_ok=True)
    Path(config["halper_output_dir_mth"]).mkdir(parents=True, exist_ok=True)
    Path(config["temp_dir"]).mkdir(parents=True, exist_ok=True)

    ortholog_script = Path(config["halper_repo"]) / "orthologFind.py"
    if not ortholog_script.exists():
        raise FileNotFoundError(f"orthologFind.py not found: {ortholog_script}")

    return config

# generates job bash job scripts
def make_job_script(config, source_species, target_species, peak_file, job_name, output_dir):
    """
    Input: config file path, source species name, target species name, peak file of source species, job name (customize), output dir.
    Output: Not that important (peak bed, summit bed, lifted peak and summit from source to target species);
            Important (SourceToTarget.HALPER.narrowPeak)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(config["temp_dir"])
    hal_file = Path(config["hal_file"])
    halper_repo = Path(config["halper_repo"])
    ortholog_script = halper_repo / "orthologFind.py"

    peak_file = Path(peak_file)
    peak_base = peak_file.stem

    job_script = temp_dir / f"{job_name}.job"
    out_log = output_dir / f"{job_name}.out.txt"
    err_log = output_dir / f"{job_name}.err.txt"

    peak_bed = output_dir / f"{peak_base}.{source_species}To{target_species}.bed"
    summit_bed = output_dir / f"{peak_base}.{source_species}To{target_species}.summits.bed"
    lifted_peak_bed = output_dir / f"{peak_base}.{source_species}To{target_species}.halLiftover.tFile.bed"
    lifted_summit_bed = output_dir / f"{peak_base}.{source_species}To{target_species}.halLiftover.sFile.bed"
    output_file = output_dir / f"{peak_base}.{source_species}To{target_species}.HALPER.narrowPeak"

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

    echo "Running HALPER: {source_species} -> {target_species}"
    echo "Peak file: {peak_file}"
    echo "HAL file: {hal_file}"

    python -c "import numpy, matplotlib; print('Python deps ok')"
    which halLiftover || true
    which python || true

    echo "Creating 4-column peak BED"
    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"peak_"NR}}' {peak_file} > {peak_bed}

    echo "Creating summit BED"
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
      -qFile {peak_bed} \\
      -tFile {lifted_peak_bed} \\
      -sFile {lifted_summit_bed} \\
      -oFile {output_file}

    if [ ! -s {output_file} ]; then
        echo "ERROR: Output file missing or empty: {output_file}" >&2
        exit 1
    fi

    echo "Finished HALPER: {source_species} -> {target_species}"
    """)

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script

# Submit created job to cluster (output/temp)
def submit_job(job_script):
    result = subprocess.run(
        ["sbatch", str(job_script)],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    print(result.stdout)
    if result.stderr:
        print(result.stderr)

def main():
    config = load_config()
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