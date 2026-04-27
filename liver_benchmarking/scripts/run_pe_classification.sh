#!/usr/bin/env bash
#SBATCH --job-name=pe_classify
#SBATCH --output=slurm-pe_classify-%j.out
#SBATCH --error=slurm-pe_classify-%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000M
# SBATCH --account=YOUR_BRIDGES2_ALLOCATION

set -euo pipefail

ROOT="/jet/home/zhuang21/TRACE_test/liver-ATAC-OCR"
IN="$ROOT/PE_classification/input"
OUT_ROOT="$ROOT/PE_classification/output"
ROW_OUT="$OUT_ROOT/rowcount"
UNIQ_OUT="$OUT_ROOT/unique"

HALPER="$IN/human_liver.HumanToMouse.HALPER.narrowPeak.gz"
MOUSE_PEAK="$IN/idr.conservative_peak.narrowPeak.gz"
TSS="$IN/gencode.vM15.annotation.protTranscript.geneNames_TSSWithStrand_sorted.bed"

mkdir -p "$ROW_OUT" "$UNIQ_OUT"

if ! command -v module >/dev/null 2>&1; then
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/Modules/init/bash
  fi
fi

if command -v module >/dev/null 2>&1; then
  module load bedtools/2.31.1 2>/dev/null || module load bedtools 2>/dev/null || true
fi

if ! command -v bedtools >/dev/null 2>&1; then
  echo "ERROR: bedtools not found. On Bridges-2, run with module load bedtools." >&2
  exit 1
fi

for f in "$HALPER" "$MOUSE_PEAK" "$TSS"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing input file: $f" >&2
    exit 1
  fi
done

dedup_by_coord() {
  local in_bed="$1"
  local out_bed="$2"
  cut -f1-3 "$in_bed" | sort -k1,1 -k2,2n -k3,3n | uniq > "$out_bed"
}

echo "Step 1: sort HALPER output and mouse native peaks ..."
gzip -dc "$HALPER" | sort -k1,1 -k2,2n > "$ROW_OUT/halper_mouse.sorted.bed"
gzip -dc "$MOUSE_PEAK" | sort -k1,1 -k2,2n > "$ROW_OUT/mouse_native.sorted.bed"

echo "Step 2: create promoter regions from TSS (TSS ± 2 kb) ..."
awk 'BEGIN{OFS="\t"} {s=$2-2000; if(s<0)s=0; e=$3+2000; print $1,s,e}' "$TSS" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - \
  > "$ROW_OUT/mouse_promoter_pm2kb.merged.bed"

echo "Step 3: identify shared-open OCRs and species-specific OCRs ..."
bedtools intersect -u \
  -a "$ROW_OUT/halper_mouse.sorted.bed" \
  -b "$ROW_OUT/mouse_native.sorted.bed" \
  > "$ROW_OUT/shared_open.bed"

bedtools intersect -v \
  -a "$ROW_OUT/halper_mouse.sorted.bed" \
  -b "$ROW_OUT/mouse_native.sorted.bed" \
  > "$ROW_OUT/human_open_mouse_closed.bed"

bedtools intersect -v \
  -a "$ROW_OUT/mouse_native.sorted.bed" \
  -b "$ROW_OUT/halper_mouse.sorted.bed" \
  > "$ROW_OUT/mouse_open_human_closed.bed"

echo "Step 4: classify shared OCRs into promoter-like vs enhancer-like ..."
bedtools intersect -u \
  -a "$ROW_OUT/shared_open.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/shared_promoter.bed"

bedtools intersect -v \
  -a "$ROW_OUT/shared_open.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/shared_enhancer.bed"

echo "Step 5: classify human-specific OCRs into promoter-like vs enhancer-like ..."
bedtools intersect -u \
  -a "$ROW_OUT/human_open_mouse_closed.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/human_specific_promoter.bed"

bedtools intersect -v \
  -a "$ROW_OUT/human_open_mouse_closed.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/human_specific_enhancer.bed"

echo "Step 6: classify mouse-specific OCRs into promoter-like vs enhancer-like ..."
bedtools intersect -u \
  -a "$ROW_OUT/mouse_open_human_closed.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/mouse_specific_promoter.bed"

bedtools intersect -v \
  -a "$ROW_OUT/mouse_open_human_closed.bed" \
  -b "$ROW_OUT/mouse_promoter_pm2kb.merged.bed" \
  > "$ROW_OUT/mouse_specific_enhancer.bed"

echo "Step 7: summarize row-based counts (reference) ..."
shared_total=$(wc -l < "$ROW_OUT/shared_open.bed")
shared_prom=$(wc -l < "$ROW_OUT/shared_promoter.bed")
shared_enh=$(wc -l < "$ROW_OUT/shared_enhancer.bed")

human_spec_total=$(wc -l < "$ROW_OUT/human_open_mouse_closed.bed")
human_spec_prom=$(wc -l < "$ROW_OUT/human_specific_promoter.bed")
human_spec_enh=$(wc -l < "$ROW_OUT/human_specific_enhancer.bed")

mouse_spec_total=$(wc -l < "$ROW_OUT/mouse_open_human_closed.bed")
mouse_spec_prom=$(wc -l < "$ROW_OUT/mouse_specific_promoter.bed")
mouse_spec_enh=$(wc -l < "$ROW_OUT/mouse_specific_enhancer.bed")

shared_prom_frac=$(awk -v p="$shared_prom" -v t="$shared_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
shared_enh_frac=$(awk -v e="$shared_enh" -v t="$shared_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

human_prom_frac=$(awk -v p="$human_spec_prom" -v t="$human_spec_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
human_enh_frac=$(awk -v e="$human_spec_enh" -v t="$human_spec_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

mouse_prom_frac=$(awk -v p="$mouse_spec_prom" -v t="$mouse_spec_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
mouse_enh_frac=$(awk -v e="$mouse_spec_enh" -v t="$mouse_spec_total" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

{
  echo -e "category\ttotal\tpromoter\tenhancer\tpromoter_fraction\tenhancer_fraction"
  echo -e "shared_open\t$shared_total\t$shared_prom\t$shared_enh\t$shared_prom_frac\t$shared_enh_frac"
  echo -e "human_open_mouse_closed\t$human_spec_total\t$human_spec_prom\t$human_spec_enh\t$human_prom_frac\t$human_enh_frac"
  echo -e "mouse_open_human_closed\t$mouse_spec_total\t$mouse_spec_prom\t$mouse_spec_enh\t$mouse_prom_frac\t$mouse_enh_frac"
} > "$ROW_OUT/summary.tsv"

echo "Step 8: deduplicate BED files by genomic coordinates and summarize unique OCR counts ..."
dedup_by_coord "$ROW_OUT/shared_open.bed" "$UNIQ_OUT/shared_open.bed"
dedup_by_coord "$ROW_OUT/human_open_mouse_closed.bed" "$UNIQ_OUT/human_open_mouse_closed.bed"
dedup_by_coord "$ROW_OUT/mouse_open_human_closed.bed" "$UNIQ_OUT/mouse_open_human_closed.bed"
dedup_by_coord "$ROW_OUT/shared_promoter.bed" "$UNIQ_OUT/shared_promoter.bed"
dedup_by_coord "$ROW_OUT/shared_enhancer.bed" "$UNIQ_OUT/shared_enhancer.bed"
dedup_by_coord "$ROW_OUT/human_specific_promoter.bed" "$UNIQ_OUT/human_specific_promoter.bed"
dedup_by_coord "$ROW_OUT/human_specific_enhancer.bed" "$UNIQ_OUT/human_specific_enhancer.bed"
dedup_by_coord "$ROW_OUT/mouse_specific_promoter.bed" "$UNIQ_OUT/mouse_specific_promoter.bed"
dedup_by_coord "$ROW_OUT/mouse_specific_enhancer.bed" "$UNIQ_OUT/mouse_specific_enhancer.bed"

shared_total_u=$(wc -l < "$UNIQ_OUT/shared_open.bed")
shared_prom_u=$(wc -l < "$UNIQ_OUT/shared_promoter.bed")
shared_enh_u=$(wc -l < "$UNIQ_OUT/shared_enhancer.bed")

human_spec_total_u=$(wc -l < "$UNIQ_OUT/human_open_mouse_closed.bed")
human_spec_prom_u=$(wc -l < "$UNIQ_OUT/human_specific_promoter.bed")
human_spec_enh_u=$(wc -l < "$UNIQ_OUT/human_specific_enhancer.bed")

mouse_spec_total_u=$(wc -l < "$UNIQ_OUT/mouse_open_human_closed.bed")
mouse_spec_prom_u=$(wc -l < "$UNIQ_OUT/mouse_specific_promoter.bed")
mouse_spec_enh_u=$(wc -l < "$UNIQ_OUT/mouse_specific_enhancer.bed")

shared_prom_frac_u=$(awk -v p="$shared_prom_u" -v t="$shared_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
shared_enh_frac_u=$(awk -v e="$shared_enh_u" -v t="$shared_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

human_prom_frac_u=$(awk -v p="$human_spec_prom_u" -v t="$human_spec_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
human_enh_frac_u=$(awk -v e="$human_spec_enh_u" -v t="$human_spec_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

mouse_prom_frac_u=$(awk -v p="$mouse_spec_prom_u" -v t="$mouse_spec_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", p/t}')
mouse_enh_frac_u=$(awk -v e="$mouse_spec_enh_u" -v t="$mouse_spec_total_u" 'BEGIN{if(t==0) print "NA"; else printf "%.4f", e/t}')

{
  echo -e "category\ttotal\tpromoter\tenhancer\tpromoter_fraction\tenhancer_fraction"
  echo -e "shared_open\t$shared_total_u\t$shared_prom_u\t$shared_enh_u\t$shared_prom_frac_u\t$shared_enh_frac_u"
  echo -e "human_open_mouse_closed\t$human_spec_total_u\t$human_spec_prom_u\t$human_spec_enh_u\t$human_prom_frac_u\t$human_enh_frac_u"
  echo -e "mouse_open_human_closed\t$mouse_spec_total_u\t$mouse_spec_prom_u\t$mouse_spec_enh_u\t$mouse_prom_frac_u\t$mouse_enh_frac_u"
} > "$UNIQ_OUT/summary.tsv"

echo "Done."
echo "Row-based output folder: $ROW_OUT"
echo "Unique-based output folder: $UNIQ_OUT"
