#!/bin/bash
## Merge the fragment files of one simulated mixture into a .pat file, keep a
## random half of its reads, and convert it to a .beta file with wgbs_tools
## (https://github.com/nloyfer/wgbs_tools), which must be on PATH.
##
##   generate_beta.sh <fragment_dir> <param> <replicate>
##
## Reads fragments_<param>_<replicate>_{1..14}.txt, one file per cell type with
## a non-zero proportion, as written by generate_pseudo_plasma.R.
set -eo pipefail
if [ $# -ne 3 ]; then
  echo "Usage: generate_beta.sh <fragment_dir> <param> <replicate>" >&2
  exit 1
fi
cd "$1"
c="$2"
d="$3"

files=(fragments_${c}_${d}_{1..14}.txt)
for f in "${files[@]}"; do
  [[ -f $f ]] || { echo "Missing $f" >&2; exit 2; }
done

out="merged_${c}_${d}_sampled_new.pat"
awk '
BEGIN {OFS="\t"}
{
  chr=$1; sub(/^chr/,"",chr)
  if      (chr=="X") key=23
  else if (chr=="Y") key=24
  else if (chr=="M") key=25
  else               key=chr+0
  print key,$0
}' "${files[@]}" \
| sort -t$'\t' -k1,1n -k3,3n \
| cut -f2- \
> "$out"
awk 'BEGIN{srand();} rand() < 0.5' "$out" > "${out}.tmp" && mv "${out}.tmp" "$out"
if command -v pigz &>/dev/null; then
  pigz -f "$out"
else
  gzip -f "$out"
fi

wgbstools pat2beta "${out}.gz"
