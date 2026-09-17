#!/usr/bin/env bash
set -euo pipefail

MERGED_DIR="/tahoma/emsl61054/AMG_eco-evo_dir/master_annotation_for_all_v3_dir/"
PY_FUNC="function_filter.py"
PY_TAXA="taxonomy_filter_v7.py"

KO="NA"
PF="NA"
TARGETS="Viruses,Bacteria,Archaea,Fungi,Viridiplantae,Metazoa,SAR"
TAXA_MODE="lenient"
CHAIN="yes"
OUTDIR="."
TAXA_GREEDY="no"
SKIP_TAXA="no"

SAMPLE_LIST=""
LOG_PREFIX="[batch]"

ts() { date +"%Y-%m-%d_%H:%M:%S"; }
usage() { echo "Usage: $0 <sample_list_file> [options]"; exit 1; }

[[ $# -lt 1 ]] && usage
SAMPLE_LIST="$1"; shift

while [[ $# -gt 0 ]]; do
	case "$1" in
		--merged-dir) MERGED_DIR="${2:?}"; shift 2;;
		--ko) KO="${2:?}"; shift 2;;
		--pf) PF="${2:?}"; shift 2;;
		--targets) TARGETS="${2:?}"; shift 2;;
		--taxa-mode) TAXA_MODE="${2:?}"; shift 2;;
		--chain) CHAIN="${2:?}"; shift 2;;
		--outdir) OUTDIR="${2:?}"; shift 2;;
		--taxa-greedy) TAXA_GREEDY="${2:?}"; shift 2;;
		--skip-taxa) SKIP_TAXA="${2:?}"; shift 2;;
		--py-func) PY_FUNC="${2:?}"; shift 2;;
		--py-taxa) PY_TAXA="${2:?}"; shift 2;;
		-h|--help) usage;;
		*) echo "Unknown arg: $1"; usage;;
	esac
done

if [[ ! -f "$SAMPLE_LIST" ]]; then echo "$(ts) ${LOG_PREFIX} sample list not found: $SAMPLE_LIST" >&2; exit 2; fi
if [[ ! -d "$MERGED_DIR" ]]; then echo "$(ts) ${LOG_PREFIX} merged dir not found: $MERGED_DIR" >&2; exit 2; fi

FUNC_DIR="${OUTDIR%/}/func_filter"
CHAIN_DIR="${OUTDIR%/}/func_then_taxa"
mkdir -p "$FUNC_DIR" "$CHAIN_DIR"

echo "$(ts) ${LOG_PREFIX} start"
echo "$(ts)       list=$SAMPLE_LIST"
echo "$(ts)       merged_dir=${MERGED_DIR%/}/"
echo "$(ts)       KO=$KO PF=$PF"
echo "$(ts)       targets=$TARGETS"
echo "$(ts)       taxa_mode=$TAXA_MODE greedy=$TAXA_GREEDY skip_taxa=$SKIP_TAXA"
echo "$(ts)       chain=$CHAIN"
echo "$(ts)       outdir=$OUTDIR"
echo "$(ts)       py_func=$PY_FUNC"
echo "$(ts)       py_taxa=$PY_TAXA"

ok=0
skipped=0
missing_merged=0
func_hits=0
func_zero=0

first_col() { awk -F',' '{print $1}' "$1" | awk '{print $1}'; }

find_merged_file() {
	local sample="$1"
	local path="${MERGED_DIR%/}/${sample}-merged_kofam_pfam_union_gtdb_nr_viral.txt"
	if [[ -f "$path" ]]; then echo "$path"; return 0; fi
	local any
	any="$(ls ${MERGED_DIR%/}/${sample}-merged_kofam_pfam_* 2>/dev/null | head -n1 || true)"
	if [[ -n "$any" && -f "$any" ]]; then echo "$any"; return 0; fi
	return 1
}

strict_flag=()
if [[ "$TAXA_MODE" == "strict" ]]; then strict_flag+=(--strict-taxid-only); fi

while IFS= read -r sample; do
	sample="${sample//[$'\r\n\t ']/}"
	[[ -z "$sample" ]] && continue

	echo "$(ts) [RUN ] sample=${sample}"

	merged_file="$(find_merged_file "$sample" || true)"
	if [[ -z "$merged_file" ]]; then
		echo "$(ts)       merged file missing"
		missing_merged=$((missing_merged+1))
		skipped=$((skipped+1))
		continue
	fi
	echo "$(ts)       merged_file=${merged_file}"

	# optional, quick PF sanity count
	pf_count=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="pfamsearch:target_accession") p=i; next} p && $p ~ /'"${PF}"'/{c++} END{print c+0}' "$merged_file")
	echo "$(ts)       [INFO] PF ${PF} rows in merged: ${pf_count}"

	# function filter
	func_tag="KO=${KO}, PF=${PF}"
	expect_new="${FUNC_DIR%/}/${sample}_func_KO-${KO}_PF-${PF}.txt"
	expect_legacy="${FUNC_DIR%/}/${sample}-filtered_KO-${KO}_PF-${PF}.txt"

	echo "$(ts)       function filter ${func_tag} -> ${expect_new}"
	func_cmd=(python3 "$PY_FUNC" --merged "$merged_file" --ko "$KO" --pf "$PF" --outdir "$FUNC_DIR")
	if ! "${func_cmd[@]}" >/dev/null 2>&1; then
		echo "$(ts)       function filter error, skipping sample"
		skipped=$((skipped+1))
		continue
	fi

	func_file=""
	if [[ -s "$expect_new" ]]; then
		func_file="$expect_new"
	elif [[ -s "$expect_legacy" ]]; then
		func_file="$expect_legacy"
	else
		found="$(ls "${FUNC_DIR%/}/${sample}"*.txt 2>/dev/null | head -n1 || true)"
		[[ -n "$found" ]] && func_file="$found"
	fi

	if [[ "$CHAIN" == "yes" ]]; then
		if [[ -z "$func_file" ]]; then
			echo "$(ts)       [INFO] function step wrote no file, zero matches, skipping taxonomy in chain mode"
			func_zero=$((func_zero+1))
			skipped=$((skipped+1))
			continue
		fi
		if [[ ! -s "$func_file" ]]; then
			echo "$(ts)       [INFO] function file exists but is empty, skipping taxonomy in chain mode"
			func_zero=$((func_zero+1))
			skipped=$((skipped+1))
			continue
		fi
		echo "$(ts)       [OK] function step found hits -> ${func_file}"
		func_hits=$((func_hits+1))
	fi

	# choose input for taxonomy
	if [[ "$CHAIN" == "yes" ]]; then
		input_for_taxa="$func_file"
	else
		input_for_taxa="$merged_file"
	fi

	target_tag="${TARGETS//,/+}"
	taxa_out="${CHAIN_DIR%/}/${sample}-filtered_taxa-${target_tag}.txt"

	if [[ "$CHAIN" == "yes" ]]; then
		echo "$(ts)       taxa filter chained after function"
	else
		echo "$(ts)       taxa filter on merged"
	fi

	taxa_cmd=(python3 "$PY_TAXA"
		--merged "$input_for_taxa"
		--sample "$sample"
		--targets "$TARGETS"
		"${strict_flag[@]}"
		--use_greedy_filter "$TAXA_GREEDY"
		--skip_taxa_filter "$SKIP_TAXA"
		--outdir "$CHAIN_DIR"
	)
	if ! "${taxa_cmd[@]}" >/dev/null 2>&1; then
		echo "$(ts)       taxa filter error, skipping sample"
		skipped=$((skipped+1))
		continue
	fi

	if [[ -s "$taxa_out" ]]; then
		ok=$((ok+1))
	elif compgen -G "${CHAIN_DIR%/}/${sample}-filtered_taxa-*.txt" > /dev/null; then
		ok=$((ok+1))
	else
		echo "$(ts)       [WARN] taxa produced no output file"
		skipped=$((skipped+1))
	fi

done < <(first_col "$SAMPLE_LIST")

echo "$(ts) [DONE] success=${ok} skipped=${skipped} missing_merged=${missing_merged} func_hits=${func_hits} func_zero=${func_zero}"