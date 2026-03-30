# run with command: python3 /ocean/projects/bio230007p/jji5/pipeline/ATAC_analysis.py

from pathlib import Path
import os
import yaml
import subprocess


# First load config and check file path. 
def load_and_validate_config():
    script_dir = Path(__file__).resolve().parent
    config_path = script_dir / "config.yaml"

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    print("Config loaded successfully.")

    paths_to_check = [
        ("halper_script", config["halper_script"]),
        ("hal_file", config["hal_file"]),
        ("species_1_peak_file", config["species_1_peak_file"]),
        ("species_2_peak_file", config["species_2_peak_file"]),
    ]

    for label, p in paths_to_check:
        path = Path(p)
        print(f"{label}: exists={path.exists()} | {path}")
        if not path.exists():
            raise FileNotFoundError(f"Missing required path for {label}: {path}")

    output_dir = Path(config["halper_output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output dir ready: {output_dir}")

    return config

if __name__ == "__main__":
    config = load_and_validate_config()


def build_halper_command(config):
    cmd = [
        "bash",
        config["halper_script"],
        "-b", config["species_1_peak_file"],
        "-o", config["halper_output_dir"],
        "-s", config["species_1"],
        "-t", config["species_2"],
        "-c", config["hal_file"],
    ]
    return cmd


def run_halper(config):
    env = os.environ.copy()

    hal_bin = "/ocean/projects/bio230007p/jji5/tools/hal/bin"
    user_bin = os.path.expanduser("~/bin")

    env["PATH"] = f"{user_bin}:{hal_bin}:{env['PATH']}"

    cmd = build_halper_command(config)

    print("\nRunning HALPER:")
    print(" ".join(cmd))
    print("PATH =", env["PATH"])

    subprocess.run(cmd, check=True, env=env)

if __name__ == "__main__":
    config = load_and_validate_config()
    run_halper(config)