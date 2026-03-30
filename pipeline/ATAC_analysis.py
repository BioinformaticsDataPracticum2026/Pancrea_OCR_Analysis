# run with command: python3 /ocean/projects/bio230007p/jji5/pipeline/ATAC_analysis.py
# sbatch /ocean/projects/bio230007p/jji5/output/temp/halper_human_pancreas_to_mouse.job

from dataclasses import dataclass
from pathlib import Path
import subprocess
import yaml

@dataclass
class HalperOneRunConfig:
    source_species: str
    target_species: str
    organ: str
    halper_script: Path
    hal_file: Path
    peak_file: Path
    output_dir: Path
    temp_dir: Path
    job_name: str = "halper_pancreas"

    # check if halper script and files exist.
    def __post_init__(self):
        if not self.halper_script.exists():
            raise FileNotFoundError(f"HALPER script not found: {self.halper_script}")
        if not self.hal_file.exists():
            raise FileNotFoundError(f"HAL file not found: {self.hal_file}")
        if not self.peak_file.exists():
            raise FileNotFoundError(f"Peak file not found: {self.peak_file}")
        if self.peak_file.suffix != ".narrowPeak":
            raise ValueError(f"Peak file must be .narrowPeak: {self.peak_file}")

        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)

# load config file -> species names, halper_script path, data path. 
def load_config(config_path: Path) -> HalperOneRunConfig:
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    return HalperOneRunConfig(
        source_species=config["species_1"],
        target_species=config["species_2"],
        organ=config["organ"],
        halper_script=Path(config["halper_script"]),
        hal_file=Path(config["hal_file"]),
        peak_file=Path(config["species_1_peak_file"]),
        output_dir=Path(config["halper_output_dir"]),
        temp_dir=Path(config["temp_dir"]),
        job_name=config.get("job_name", "halper_pancreas")
    )

# Generate sbatch script based on names and path given in config. 
def generate_slurm_script(config: HalperOneRunConfig) -> Path:
    job_script = config.temp_dir / f"{config.job_name}.job"
    out_log = config.output_dir / f"{config.job_name}.out.txt"
    err_log = config.output_dir / f"{config.job_name}.err.txt"

    # see output/temp/halper_human_pancreas_mouse.job for the script generated.
    script_text = f"""#!/bin/bash
        #SBATCH -p RM-shared
        #SBATCH -t 12:00:00
        #SBATCH --ntasks-per-node=4
        #SBATCH --job-name={config.job_name}
        #SBATCH --output={out_log}
        #SBATCH --error={err_log}

        set -euo pipefail

        export PATH=/ocean/projects/bio230007p/jji5/tools/hal/bin:$PATH
        export PATH=$HOME/bin:$PATH

        echo "Running HALPER"
        echo "Source species: {config.source_species}"
        echo "Target species: {config.target_species}"
        echo "Organ: {config.organ}"
        echo "Peak file: {config.peak_file}"
        echo "HAL file: {config.hal_file}"
        echo "Output dir: {config.output_dir}"

        which halLiftover || true
        which python || true

        bash {config.halper_script} \\
        -b {config.peak_file} \\
        -o {config.output_dir} \\
        -s {config.source_species} \\
        -t {config.target_species} \\
        -c {config.hal_file}

        echo "HALPER job finished"
        """

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script

# submit job to cluster. 
def submit_job(job_script: Path) -> None:
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
    script_dir = Path(__file__).resolve().parent
    config_path = script_dir / "config.yaml"

    config = load_config(config_path)
    print("Config loaded successfully.")

    job_script = generate_slurm_script(config)
    print(f"Generated SLURM script: {job_script}")

    submit_job(job_script)


if __name__ == "__main__":
    main()