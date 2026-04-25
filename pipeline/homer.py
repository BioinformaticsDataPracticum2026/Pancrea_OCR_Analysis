from pathlib import Path
import textwrap
from utils import load_config, submit_job

def make_motif_job(
    foreground,
    genome_fasta,
    output_dir,
    temp_dir,
    job_name,
    background=None,
):
    foreground = Path(foreground)
    genome_fasta = Path(genome_fasta)
    output_dir = Path(output_dir)
    temp_dir = Path(temp_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    job_script = temp_dir / f"{job_name}.job"

    log_dir = temp_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    out_log = log_dir / f"{job_name}.out.txt"
    err_log = log_dir / f"{job_name}.err.txt"

    if background:
        background = Path(background)
        motif_cmd = (
            f"findMotifsGenome.pl {foreground} {genome_fasta} {output_dir} "
            f"-size 200 -p 4 -bg {background}"
        )
        background_echo = f'echo "Background: {background}"'
    else:
        motif_cmd = (
            f"findMotifsGenome.pl {foreground} {genome_fasta} {output_dir} "
            f"-size 200 -p 4"
        )
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
    echo ""

    {motif_cmd}

    echo "HOMER motif analysis finished"
    """)

    with open(job_script, "w") as f:
        f.write(script_text)

    job_script.chmod(0o755)
    return job_script


def main():
    config = load_config(
        path_keys=[
            "mouse_enhancer_file",
            "human_enhancer_file",
            "mouse_promoter_file",
            "human_promoter_file",
            "shared_enhancer_htm_file",
            "shared_enhancer_mth_file",
            "human_specific_enhancer_htm_file",
            "mouse_specific_enhancer_mth_file",
            "mouse_genome_fasta",
            "human_genome_fasta",
            "homer_output_dir",
            "homer_temp_dir",
        ],
        required_paths=[
            "mouse_enhancer_file",
            "human_enhancer_file",
            "mouse_promoter_file",
            "human_promoter_file",
            "shared_enhancer_htm_file",
            "shared_enhancer_mth_file",
            "human_specific_enhancer_htm_file",
            "mouse_specific_enhancer_mth_file",
            "mouse_genome_fasta",
            "human_genome_fasta",
        ],
        mkdir_keys=[
            "homer_output_dir",
            "homer_temp_dir",
        ],
    )

    print("HOMER config loaded successfully.")

    homer_output_dir = Path(config["homer_output_dir"])

    jobs = [
        # Native enhancer sets
        make_motif_job(
            foreground=config["mouse_enhancer_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=homer_output_dir / "mouse_enhancer",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_mouse_enhancer",
        ),
        make_motif_job(
            foreground=config["human_enhancer_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=homer_output_dir / "human_enhancer",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_human_enhancer",
        ),

        # Native promoter sets
        make_motif_job(
            foreground=config["mouse_promoter_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=homer_output_dir / "mouse_promoter",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_mouse_promoter",
        ),
        make_motif_job(
            foreground=config["human_promoter_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=homer_output_dir / "human_promoter",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_human_promoter",
        ),

        # HtM: shared enhancer and human-specific enhancer
        # These are in mouse coordinates, so use mouse genome.
        make_motif_job(
            foreground=config["shared_enhancer_htm_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=homer_output_dir / "shared_enhancer_htm",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_shared_enhancer_htm",
            background=config["human_specific_enhancer_htm_file"],
        ),
        make_motif_job(
            foreground=config["human_specific_enhancer_htm_file"],
            genome_fasta=config["mouse_genome_fasta"],
            output_dir=homer_output_dir / "human_specific_enhancer_htm",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_human_specific_enhancer_htm",
            background=config["shared_enhancer_htm_file"],
        ),

        # MtH: shared enhancer and mouse-specific enhancer
        # These are in human coordinates, so use human genome.
        make_motif_job(
            foreground=config["shared_enhancer_mth_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=homer_output_dir / "shared_enhancer_mth",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_shared_enhancer_mth",
            background=config["mouse_specific_enhancer_mth_file"],
        ),
        make_motif_job(
            foreground=config["mouse_specific_enhancer_mth_file"],
            genome_fasta=config["human_genome_fasta"],
            output_dir=homer_output_dir / "mouse_specific_enhancer_mth",
            temp_dir=config["homer_temp_dir"],
            job_name="homer_mouse_specific_enhancer_mth",
            background=config["shared_enhancer_mth_file"],
        ),
    ]

    for job in jobs:
        print(f"Generated: {job}")

    for job in jobs:
        submit_job(job)

if __name__ == "__main__":
    main()