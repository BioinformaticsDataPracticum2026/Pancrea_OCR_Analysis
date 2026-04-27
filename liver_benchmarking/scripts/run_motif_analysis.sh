#!/usr/bin/bash
#SBATCH -J motif_analysis
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH -o /jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/task5_motif_%j.out
#SBATCH -e /jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/logs/task5_motif_%j.err

set -euo pipefail

module purge
module load bedtools/2.30.0
module load MEME-suite/5.4.1

ROOT="/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR"
PE_OUT="$ROOT/PE_classification/output"
ROWCOUNT="$PE_OUT/rowcount"
UNIQUE="$PE_OUT/unique"

OUT="$ROOT/Motif_analysis/output"
BEDS="$OUT/beds"
FASTAS="$OUT/fastas"

MM10_GENOME="/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR/data/genome/mm10.fa"
JASPAR_DB="$ROOT/Motif_analysis/motif_dbs/JASPAR2024_vertebrates.meme"
WINDOW=100

mkdir -p "$BEDS" "$FASTAS" "$OUT/logs"

echo "Performing Motif Analysis"
echo "Root:       $ROOT"
echo "Genome:     $MM10_GENOME"
echo "JASPAR:     $JASPAR_DB"
echo ""

if [[ ! -f "$MM10_GENOME" ]]; then
  echo "ERROR: mm10 genome not found at $MM10_GENOME"; exit 1
fi
if [[ ! -f "$JASPAR_DB" ]]; then
  echo "ERROR: JASPAR DB not found at $JASPAR_DB"; exit 1
fi

echo "[Step 1] Generating mouse-native BED files..."

bedtools intersect -v \
  -a "$ROWCOUNT/mouse_native.sorted.bed" \
  -b "$ROWCOUNT/mouse_promoter_pm2kb.merged.bed" \
  | cut -f1-3 | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/mouse_native_enhancer.bed"

bedtools intersect -u \
  -a "$ROWCOUNT/mouse_native.sorted.bed" \
  -b "$ROWCOUNT/mouse_promoter_pm2kb.merged.bed" \
  | cut -f1-3 | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/mouse_native_promoter.bed"

cat "$UNIQUE/shared_enhancer.bed" "$UNIQUE/human_specific_enhancer.bed" \
  | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/human_all_enhancer.bed"

cat "$UNIQUE/shared_promoter.bed" "$UNIQUE/human_specific_promoter.bed" \
  | sort -k1,1 -k2,2n | uniq \
  > "$BEDS/human_all_promoter.bed"

cp "$UNIQUE/shared_enhancer.bed"          "$BEDS/shared_enhancer.bed"
cp "$UNIQUE/human_specific_enhancer.bed"  "$BEDS/human_specific_enhancer.bed"
cp "$UNIQUE/mouse_specific_enhancer.bed"  "$BEDS/mouse_specific_enhancer.bed"

echo "[Step 1] Region counts:"
for f in "$BEDS"/*.bed; do
  printf "  %-35s %d regions\n" "$(basename "$f")" "$(wc -l < "$f")"
done
echo ""

echo "[Step 2] Resizing peaks to $(( WINDOW * 2 ))bp windows..."

resize_bed() {
  local in_bed="$1"
  local out_bed="$2"
  awk -v w="$WINDOW" 'BEGIN{OFS="\t"} {
    mid = int(($2 + $3) / 2)
    s = mid - w
    if (s < 0) s = 0
    print $1, s, mid + w
  }' "$in_bed" | sort -k1,1 -k2,2n | uniq > "$out_bed"
}

declare -A PEAK_SETS
PEAK_SETS["human_enhancer"]="$BEDS/human_all_enhancer.bed"
PEAK_SETS["mouse_enhancer"]="$BEDS/mouse_native_enhancer.bed"
PEAK_SETS["human_promoter"]="$BEDS/human_all_promoter.bed"
PEAK_SETS["mouse_promoter"]="$BEDS/mouse_native_promoter.bed"
PEAK_SETS["shared_enhancer"]="$BEDS/shared_enhancer.bed"
PEAK_SETS["human_specific_enhancer"]="$BEDS/human_specific_enhancer.bed"
PEAK_SETS["mouse_specific_enhancer"]="$BEDS/mouse_specific_enhancer.bed"

for name in "${!PEAK_SETS[@]}"; do
  resize_bed "${PEAK_SETS[$name]}" "$BEDS/${name}.resized.bed"
  echo "  $name: $(wc -l < "$BEDS/${name}.resized.bed") windows"
done
echo ""

echo "[Step 3] Extracting FASTA sequences from mm10..."

for name in "${!PEAK_SETS[@]}"; do
  bedtools getfasta \
    -fi "$MM10_GENOME" \
    -bed "$BEDS/${name}.resized.bed" \
    -fo "$FASTAS/${name}.fa"
  echo "  $name: $(grep -c '^>' "$FASTAS/${name}.fa") sequences"
done
echo ""

echo "[Step 4] Running MEME-ChIP (this will take a while)..."

run_meme_chip() {
  local name="$1"
  local fasta="$FASTAS/${name}.fa"
  local outdir="$OUT/meme_chip_${name}"
  local logfile="$OUT/logs/${name}.log"

  echo "  [$name] Starting..."
  export OMPI_MCA_rmaps_base_oversubscribe=1
  meme-chip \
    -oc "$outdir" \
    -maxw 20 \
    -meme-nmotifs 10 \
    -meme-p 4 \
    -db "$JASPAR_DB" \
    "$fasta" \
    > "$logfile" 2>&1
  echo "  [$name] Done. Results: $outdir"
}

for name in "${!PEAK_SETS[@]}"; do
  run_meme_chip "$name"
done

echo ""
echo "===== All done ====="
echo "Results in: $OUT"
