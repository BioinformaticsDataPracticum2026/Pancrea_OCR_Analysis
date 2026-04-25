from pathlib import Path
import subprocess
import yaml


def resolve_path(project_root, path_str):
    path = Path(path_str)
    if path.is_absolute():
        return path
    return project_root / path


def load_config(path_keys=None, required_paths=None, mkdir_keys=None):
    script_path = Path(__file__).resolve()
    project_root = script_path.parents[1]
    config_path = project_root / "config.yaml"

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    if "project_root" in config:
        project_root = Path(config["project_root"]).resolve()

    path_keys = path_keys or []
    required_paths = required_paths or []
    mkdir_keys = mkdir_keys or []

    for key in path_keys:
        if key not in config:
            raise KeyError(f"Missing required config key: {key}")
        config[key] = str(resolve_path(project_root, config[key]))

    for key in required_paths:
        path = Path(config[key])
        if not path.exists():
            raise FileNotFoundError(f"{key} not found: {path}")

    for key in mkdir_keys:
        Path(config[key]).mkdir(parents=True, exist_ok=True)

    config["project_root"] = str(project_root)

    return config


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