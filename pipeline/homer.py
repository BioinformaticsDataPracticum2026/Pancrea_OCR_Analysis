from pathlib import Path
import subprocess
import textwrap
import yaml


# run with:
# python3 /ocean/projects/bio230007p/jji5/pipeline/homer.py


def load_config():
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config.yaml"

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    required_paths = [
        ("mouse_enhancer_file", Path(config["mouse_enhancer_file"])),
        ("human_enhancer_file", Path(config["human_enhancer_file"])),
        ("mouse_promoter_file", Path(config["mouse_promoter_file"])),
        ("human_promoter_file", Path(config["human_promoter_file"])),
        ("shared_enhancer_htm_file", Path(config["shared_enhancer_htm_file"])),
        ("shared_enhancer_mth_file", Path(config["shared_enhancer_mth_file"])),
        ("specific_enhancer_htm_file", Path(config["specific_enhancer_htm_file"])),
        ("specific_enhancer_mth_file", Path(config["specific_enhancer_mth_file"])),
        ("mouse_genome_fasta", Path(config["mouse_genome_fasta"])),
        ("human_genome_fasta", Path(config["human_genome_fasta"])),
    ]

    for label, path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    Path(config["homer_output_dir"]).mkdir(parents=True, exist_ok=True)
    Path(config["temp_dir"]).mkdir(parents=True, exist_ok=True)

    return config

def make_motif_job(foreground, genome_fasta, output_dir, temp_dir, job_name, background=None):
    foreground = Path(foreground)
    genome_fasta = Path(genome_fasta)
    output_dir = Path(output_dir)
    temp_dir = Path(temp_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    job_script = temp_dir / f"{job_name}.job"
    out_log = output_dir / f"{job_name}.out.txt"
    err_log = output_dir / f"{job_name}.err.txt"

    if background:
        motif_cmd = f"findMotifsGenome.pl {foreground} {genome_fasta} {output_dir} -size 200 -p 4 -bg {background}"
        background_echo = f'echo "Background: {background}"'
    else:
        motif_cmd = f"findMotifsGenome.pl {foreground} {genome_fasta} {output_dir} -size 200 -p 4"
        background_echo = 'echo "Background: default HOMER background"'

    script_text = textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH -p RM-shared
    #SBATCH -t 08:00:00
    #SBATCH --ntasks-per-node=4
    #SBATCH --job-name={job_name}
    #SBATCH --output={out_log}
    #SBATCH --error={err_log}

    set -euo pipefail

    module load homer

    echo "Running HOMER motif analysis"
    echo "Foreground: {foreground}"
    echo "Genome FASTA: {genome_fasta}"
    echo "Output dir: {output_dir}"
    {background_echo}

    {motif_cmd}

    echo "HOMER motif analysis finished"
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
    print("HOMER config loaded successfully.")

    jobs = [
        # enhancer for each species
        make_motif_job(
            foreground=config["mouse_enhancer_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "mouse_enhancer",
            temp_dir=config["temp_dir"],
            job_name="homer_mouse_enhancer",
        ),
        make_motif_job(
            foreground=config["human_enhancer_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "human_enhancer",
            temp_dir=config["temp_dir"],
            job_name="homer_human_enhancer",
        ),

        # promoter for each species
        make_motif_job(
            foreground=config["mouse_promoter_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "mouse_promoter",
            temp_dir=config["temp_dir"],
            job_name="homer_mouse_promoter",
        ),
        make_motif_job(
            foreground=config["human_promoter_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "human_promoter",
            temp_dir=config["temp_dir"],
            job_name="homer_human_promoter",
        ),

        # shared enhancers vs specific enhancers
        make_motif_job(
            foreground=config["shared_enhancer_htm_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "shared_enhancer_htm",
            temp_dir=config["temp_dir"],
            job_name="homer_shared_enhancer_htm",
            background=config["specific_enhancer_htm_file"],
        ),
        make_motif_job(
            foreground=config["shared_enhancer_mth_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "shared_enhancer_mth",
            temp_dir=config["temp_dir"],
            job_name="homer_shared_enhancer_mth",
            background=config["specific_enhancer_mth_file"],
        ),

        # specific enhancers vs shared enhancers
        make_motif_job(
            foreground=config["specific_enhancer_htm_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "specific_enhancer_htm",
            temp_dir=config["temp_dir"],
            job_name="homer_specific_enhancer_htm",
            background=config["shared_enhancer_htm_file"],
        ),
        make_motif_job(
            foreground=config["specific_enhancer_mth_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=Path(config["homer_output_dir"]) / "specific_enhancer_mth",
            temp_dir=config["temp_dir"],
            job_name="homer_specific_enhancer_mth",
            background=config["shared_enhancer_mth_file"],
        ),
    ]

    for job in jobs:
        print(f"Generated: {job}")

    for job in jobs:
        submit_job(job)


if __name__ == "__main__":
    main()