"""
Run part or all of the pipeline:

    python3 main.py --run halper
    python3 main.py --run bed_g homer
    python3 main.py --run rgreat all
    python3 main.py --run rgreat HtM_shared
    python3 main.py --run rgreat MtH_mouse_specific
    python3 main.py --run rgreat all plot_rgreat
    python3 main.py --run all
"""

import argparse
from pathlib import Path
import subprocess
import sys

project_root = Path(__file__).resolve().parent
sys.path.append(str(project_root / "pipeline"))

import pipeline.halper as halper
import pipeline.bed_genome as bed_genome
import pipeline.bed_promoter_enhancer as bed_promoter_enhancer
import pipeline.homer as homer


PIPELINE_STEPS = {
    "halper",
    "bed_g",
    "bed_pe",
    "homer",
    "rgreat",
    "plot_rgreat",
    "all",
}

RGREAT_TASKS = [
    "Human_full",
    "Mouse_full",
    "HtM_shared",
    "HtM_human_specific",
    "MtH_shared",
    "MtH_mouse_specific",
]

ALL_PIPELINE_STEPS = [
    "halper",
    "bed_g",
    "bed_pe",
    "homer",
    "rgreat",
    "plot_rgreat",
]


def run_r_script(script_name: str, script_args=None) -> None:
    if script_args is None:
        script_args = []

    script_path = project_root / "pipeline" / "r_scripts" / script_name

    if not script_path.exists():
        raise FileNotFoundError(f"R script not found: {script_path}")

    cmd = ["Rscript", str(script_path)] + script_args

    print(f"Running: {' '.join(cmd)}")

    subprocess.run(
        cmd,
        check=True,
        cwd=project_root,
    )


def run_rgreat(task: str) -> None:
    if task == "all":
        for rgreat_task in RGREAT_TASKS:
            run_r_script("rGREAT.R", [rgreat_task])
        return

    if task not in RGREAT_TASKS:
        valid = ", ".join(["all"] + RGREAT_TASKS)
        raise ValueError(
            f"Unknown rGREAT task: {task}\n"
            f"Valid rGREAT tasks are: {valid}"
        )

    run_r_script("rGREAT.R", [task])


def parse_run_steps(run_tokens):
    """
    Parse flexible --run input.

    Examples:
        --run all
        --run halper bed_g
        --run rgreat all
        --run rgreat HtM_shared
        --run rgreat all plot_rgreat
    """
    parsed_steps = []
    i = 0

    while i < len(run_tokens):
        token = run_tokens[i]

        if token == "all":
            parsed_steps.extend(ALL_PIPELINE_STEPS)
            i += 1
            continue

        if token == "rgreat":
            # Default behavior:
            #   --run rgreat
            # runs all rGREAT tasks.
            rgreat_task = "all"

            # Special behavior:
            #   --run rgreat HtM_shared
            #   --run rgreat all
            # consumes the next token as the rGREAT task.
            if i + 1 < len(run_tokens):
                next_token = run_tokens[i + 1]

                if next_token == "all" or next_token in RGREAT_TASKS:
                    rgreat_task = next_token
                    i += 1

            parsed_steps.append(("rgreat", rgreat_task))
            i += 1
            continue

        if token in PIPELINE_STEPS:
            parsed_steps.append((token, None))
            i += 1
            continue

        valid_pipeline = ", ".join(sorted(PIPELINE_STEPS))
        valid_rgreat = ", ".join(["all"] + RGREAT_TASKS)

        raise ValueError(
            f"Unknown --run token: {token}\n"
            f"Valid pipeline steps are: {valid_pipeline}\n"
            f"Valid rGREAT tasks after 'rgreat' are: {valid_rgreat}"
        )

    return parsed_steps


def main():
    parser = argparse.ArgumentParser(
        description="Run selected parts of the ATAC-seq comparative pipeline."
    )

    parser.add_argument(
        "--run",
        nargs="+",
        required=True,
        help=(
            "Pipeline steps to run. Examples: "
            "--run halper, "
            "--run bed_g homer, "
            "--run rgreat all, "
            "--run rgreat HtM_shared, "
            "--run all"
        ),
    )

    args = parser.parse_args()
    steps = parse_run_steps(args.run)

    for step, option in steps:
        if step == "halper":
            halper.main()

        elif step == "bed_g":
            bed_genome.main()

        elif step == "bed_pe":
            bed_promoter_enhancer.main()

        elif step == "homer":
            homer.main()

        elif step == "rgreat":
            run_rgreat(option)

        elif step == "plot_rgreat":
            run_r_script("plot_rGREAT.R")

        else:
            raise ValueError(f"Unhandled pipeline step: {step}")


if __name__ == "__main__":
    main()