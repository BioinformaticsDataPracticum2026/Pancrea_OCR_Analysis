from pathlib import Path
import subprocess
import textwrap
import yaml


# run with:
# python3 /ocean/projects/bio230007p/jji5/pipeline/bed_promoter_enhancer.py


def load_config():
    project_root = Path(__file__).resolve().parents[1]
    config_path = project_root / "config.yaml"

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    required_paths = [
        ("mapped_htm_file", Path(config["mapped_htm_file"])),
        ("mapped_mth_file", Path(config["mapped_mth_file"])),
        ("mouse_peak_file", Path(config["mouse_peak_file"])),
        ("human_peak_file", Path(config["human_peak_file"])),
        ("mouse_tss_file", Path(config["mouse_tss_file"])),
        ("human_tss_file", Path(config["human_tss_file"])),
    ]

    for label, path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    Path(config["bed_pe_output_dir_htm"]).mkdir(parents=True, exist_ok=True)
    Path(config["bed_pe_output_dir_mth"]).mkdir(parents=True, exist_ok=True)
    Path(config["temp_dir"]).mkdir(parents=True, exist_ok=True)

    return config


def make_job_script(
    mapped_file,
    target_peak_file,
    target_tss_file,
    output_dir,
    temp_dir,
    job_name,
    target_label,
    promoter_threshold,
):
    mapped_file = Path(mapped_file)
    target_peak_file = Path(target_peak_file)
    target_tss_file = Path(target_tss_file)
    output_dir = Path(output_dir)
    temp_dir = Path(temp_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    job_script = temp_dir / f"{job_name}.job"
    out_log = output_dir / f"{job_name}.out.txt"
    err_log = output_dir / f"{job_name}.err.txt"

    mapped_bed = output_dir / "mapped_halper.bed"
    target_bed = output_dir / f"{target_label}_atac.bed"

    mapped_closest = output_dir / "mapped_with_nearest_tss.bed"
    target_closest = output_dir / f"{target_label}_with_nearest_tss.bed"

    mapped_promoter = output_dir / "mapped_promoter_ocrs.bed"
    mapped_enhancer = output_dir / "mapped_enhancer_ocrs.bed"
    target_promoter = output_dir / f"{target_label}_promoter_ocrs.bed"
    target_enhancer = output_dir / f"{target_label}_enhancer_ocrs.bed"

    shared_promoter = output_dir / "shared_promoter_ocrs.bed"
    shared_enhancer = output_dir / "shared_enhancer_ocrs.bed"

    specific_promoter = output_dir / "specific_promoter_ocrs.bed"
    specific_enhancer = output_dir / "specific_enhancer_ocrs.bed"

    summary_txt = output_dir / "promoter_vs_enhancer_summary.txt"

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

    echo "Running promoter vs enhancer conservation analysis"
    echo "Mapped file: {mapped_file}"
    echo "Target peak file: {target_peak_file}"
    echo "Target TSS file: {target_tss_file}"
    echo "Target label: {target_label}"
    echo "Promoter threshold: {promoter_threshold} bp"

    # 1) Convert to BED4
    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$5}}' {mapped_file} > {mapped_bed}
    awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,"{target_label}_peak_"NR}}' {target_peak_file} > {target_bed}

    # 2) Sort inputs
    sort -k1,1 -k2,2n {mapped_bed} > {mapped_bed}.sorted
    sort -k1,1 -k2,2n {target_bed} > {target_bed}.sorted
    sort -k1,1 -k2,2n {target_tss_file} > {target_tss_file}.sorted

    # 3) Find nearest TSS and distance
    bedtools closest -a {mapped_bed}.sorted -b {target_tss_file}.sorted -d > {mapped_closest}
    bedtools closest -a {target_bed}.sorted -b {target_tss_file}.sorted -d > {target_closest}

    # 4) Classify by nearest-TSS distance
    awk 'BEGIN{{OFS="\\t"}} $NF <= {promoter_threshold} {{print $1,$2,$3,$4}}' {mapped_closest} > {mapped_promoter}
    awk 'BEGIN{{OFS="\\t"}} $NF > {promoter_threshold}  {{print $1,$2,$3,$4}}' {mapped_closest} > {mapped_enhancer}

    awk 'BEGIN{{OFS="\\t"}} $NF <= {promoter_threshold} {{print $1,$2,$3,$4}}' {target_closest} > {target_promoter}
    awk 'BEGIN{{OFS="\\t"}} $NF > {promoter_threshold}  {{print $1,$2,$3,$4}}' {target_closest} > {target_enhancer}

    # 5) Shared OCRs within class
    bedtools intersect -a {mapped_promoter} -b {target_promoter} -u > {shared_promoter}
    bedtools intersect -a {mapped_enhancer} -b {target_enhancer} -u > {shared_enhancer}

    # 6) Specific OCRs within class = mapped minus shared
    bedtools intersect -a {mapped_promoter} -b {shared_promoter} -v > {specific_promoter}
    bedtools intersect -a {mapped_enhancer} -b {shared_enhancer} -v > {specific_enhancer}

    # 7) Summarize
    mapped_promoter_n=$(wc -l < {mapped_promoter})
    mapped_enhancer_n=$(wc -l < {mapped_enhancer})
    shared_promoter_n=$(wc -l < {shared_promoter})
    shared_enhancer_n=$(wc -l < {shared_enhancer})
    specific_promoter_n=$(wc -l < {specific_promoter})
    specific_enhancer_n=$(wc -l < {specific_enhancer})

    promoter_pct=$(awk -v a="$shared_promoter_n" -v b="$mapped_promoter_n" 'BEGIN{{if(b==0) print "NA"; else printf "%.4f", 100*a/b}}')
    enhancer_pct=$(awk -v a="$shared_enhancer_n" -v b="$mapped_enhancer_n" 'BEGIN{{if(b==0) print "NA"; else printf "%.4f", 100*a/b}}')

    {{
        echo "target_label\\t{target_label}"
        echo "promoter_threshold_bp\\t{promoter_threshold}"
        echo "mapped_promoter_n\\t$mapped_promoter_n"
        echo "shared_promoter_n\\t$shared_promoter_n"
        echo "specific_promoter_n\\t$specific_promoter_n"
        echo "shared_promoter_pct\\t$promoter_pct"
        echo "mapped_enhancer_n\\t$mapped_enhancer_n"
        echo "shared_enhancer_n\\t$shared_enhancer_n"
        echo "specific_enhancer_n\\t$specific_enhancer_n"
        echo "shared_enhancer_pct\\t$enhancer_pct"
    }} > {summary_txt}

    echo "Mapped promoter OCRs: $mapped_promoter_n"
    echo "Shared promoter OCRs: $shared_promoter_n"
    echo "Specific promoter OCRs: $specific_promoter_n"
    echo "Shared promoter %: $promoter_pct"
    echo "Mapped enhancer OCRs: $mapped_enhancer_n"
    echo "Shared enhancer OCRs: $shared_enhancer_n"
    echo "Specific enhancer OCRs: $specific_enhancer_n"
    echo "Shared enhancer %: $enhancer_pct"
    echo "Promoter vs enhancer conservation analysis finished"
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
    print("Promoter/enhancer BEDTools config loaded successfully.")

    htm_job = make_job_script(
        mapped_file=config["mapped_htm_file"],
        target_peak_file=config["mouse_peak_file"],
        target_tss_file=config["mouse_tss_file"],
        output_dir=config["bed_pe_output_dir_htm"],
        temp_dir=config["temp_dir"],
        job_name="bed_pe_human_pancreas_to_mouse",
        target_label="mouse",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    mth_job = make_job_script(
        mapped_file=config["mapped_mth_file"],
        target_peak_file=config["human_peak_file"],
        target_tss_file=config["human_tss_file"],
        output_dir=config["bed_pe_output_dir_mth"],
        temp_dir=config["temp_dir"],
        job_name="bed_pe_mouse_pancreas_to_human",
        target_label="human",
        promoter_threshold=config.get("promoter_threshold", 2000),
    )

    print(f"Generated: {htm_job}")
    print(f"Generated: {mth_job}")

    submit_job(htm_job)
    submit_job(mth_job)


if __name__ == "__main__":
    main()