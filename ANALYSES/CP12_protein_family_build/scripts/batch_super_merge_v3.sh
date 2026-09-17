#!/usr/bin/env bash
set -euo pipefail

# -------- config, edit if needed --------
KOF_DIR="/tahoma/emsl61054/AMG_eco-evo_dir/renamed_sample_pfam_kofam_dir/KO_out_dir"
PFAM_DIR="/tahoma/emsl61054/AMG_eco-evo_dir/renamed_sample_pfam_kofam_dir/pfam_out_dir"
GTDB_DIR="/tahoma/emsl61054/AMG_eco-evo_dir/renamed_1kb_contig_CAT_dir/gtdb_lineage_out_dir"
NR_DIR="/tahoma/emsl61054/AMG_eco-evo_dir/renamed_1kb_contig_CAT_dir/nr_lineage_withgtdb_format_out_dir"
VIRUS_DIR="/tahoma/emsl61054/AMG_eco-evo_dir/renamed_virus_summary_dir"

PY_SCRIPT="super_merge_v10.py"
# ---------------------------------------

# usage: run_all_clean_merge_cat_virus.sh samples.txt [outer|inner] [outdir]
if [[ $# -lt 1 || $# -gt 3 ]]; then
	echo "Usage: $0 <samples.txt> [outer|inner] [outdir]" >&2
	exit 1
fi

SAMPLE_LIST="$1"
JOIN_TYPE="${2:-outer}"
OUTDIR="${3:-}"

[[ "$JOIN_TYPE" == "outer" || "$JOIN_TYPE" == "inner" ]] || { echo "Join must be outer or inner" >&2; exit 1; }
[[ -f "$SAMPLE_LIST" ]] || { echo "Sample list not found, $SAMPLE_LIST" >&2; exit 1; }
[[ -f "$PY_SCRIPT"   ]] || { echo "Python script not found, $PY_SCRIPT" >&2; exit 1; }
[[ -d "$KOF_DIR"     ]] || { echo "KOF dir not found, $KOF_DIR" >&2; exit 1; }
[[ -d "$PFAM_DIR"    ]] || { echo "PFAM dir not found, $PFAM_DIR" >&2; exit 1; }
[[ -z "$GTDB_DIR" || -d "$GTDB_DIR" ]] || { echo "GTDB dir not found, $GTDB_DIR" >&2; exit 1; }
[[ -z "$NR_DIR"   || -d "$NR_DIR"   ]] || { echo "NR dir not found, $NR_DIR" >&2; exit 1; }
[[ -z "$VIRUS_DIR" || -d "$VIRUS_DIR" ]] || { echo "VIRUS dir not found, $VIRUS_DIR" >&2; exit 1; }

# output dir and log
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
base_outdir="${OUTDIR:-$script_dir}"
mkdir -p "$base_outdir"
ts="$(date +%Y%m%d-%H%M%S)"
LOG="${base_outdir}/batch_clean_merge_${ts}.log"
: > "$LOG"

log(){ printf "%s %s\n" "$(date +%F_%T)" "$*" | tee -a "$LOG" >&2; }

# escape for ERE
escape_ere(){ printf '%s' "$1" | sed -e 's/[][^$.*/+?(){}/|]/\\&/g'; }

# hyphen aware boundary match, boundaries are non alnum, so '-' is a separator
# optional require_substr to narrow file types
choose_best_file(){
	dir="$1"; token_raw="$2"; require_substr="${3:-}"
	[[ -z "$token_raw" || -z "$dir" || ! -d "$dir" ]] && { printf ""; return; }

	token_ere="$(escape_ere "$token_raw")"
	pat="(^|[^[:alnum:]])${token_ere}([^[:alnum:]]|$)"

	# case sensitive pass, newest first
	while IFS= read -r f; do
		[[ -z "$f" ]] && continue
		base="${f##*/}"
		[[ -n "$require_substr" && "$base" != *"$require_substr"* ]] && continue
		if printf '%s' "$base" | grep -Eq "$pat"; then
			printf '%s\n' "$f"
			return
		fi
	done < <(ls -1t "$dir" 2>/dev/null | awk -v d="$dir" '{print d "/" $0}')

	# case insensitive fallback
	lpat="$(printf '%s' "$pat" | tr '[:upper:]' '[:lower:]')"
	while IFS= read -r f; do
		[[ -z "$f" ]] && continue
		base="${f##*/}"
		[[ -n "$require_substr" && "$base" != *"$require_substr"* ]] && continue
		if printf '%s' "$base" | tr '[:upper:]' '[:lower:]' | grep -Eq "$lpat"; then
			printf '%s\n' "$f"
			return
		fi
	done < <(ls -1t "$dir" 2>/dev/null | awk -v d="$dir" '{print d "/" $0}')

	printf ""
}

suggest_nearby(){
	dir="$1"; token_raw="$2"; require_substr="${3:-}"
	core="$token_raw"
	if [[ "$token_raw" =~ ^([0-9]{8}-PNNL-[A-Za-z0-9]+) ]]; then
		core="${BASH_REMATCH[1]}"
	elif [[ "$token_raw" =~ ^(PNNL-[A-Za-z0-9]+) ]]; then
		core="${BASH_REMATCH[1]}"
	fi
	if [[ -n "$require_substr" ]]; then
		ls -1t "$dir" 2>/dev/null | grep -F "$require_substr" | grep -E "$core" | head -3 | paste -sd' ' - || true
	else
		ls -1t "$dir" 2>/dev/null | grep -E "$core" | head -3 | paste -sd' ' - || true
	fi
}

ok=0; skipped=0; miss_kof=0; miss_pfa=0

log "[START] join=${JOIN_TYPE} outdir=${base_outdir}"
log "[PATHS] kof_dir=${KOF_DIR} pfam_dir=${PFAM_DIR} gtdb_dir=${GTDB_DIR:-NA} nr_dir=${NR_DIR:-NA} virus_dir=${VIRUS_DIR:-NA}"
log "[NOTE]  sample list supports 'sample' or 'sample,kofam_token,pfam_token,gtdb_token,nr_token,viral_token'"

# read list with CR stripping
while IFS= read -r line || [[ -n "$line" ]]; do
	line="${line%$'\r'}"
	line="${line#"${line%%[![:space:]]*}"}"
	line="${line%"${line##*[![:space:]]}"}"
	[[ -z "$line" || "${line:0:1}" == "#" ]] && continue

	# parse up to 6 CSV fields
	IFS=',' read -r sample tok_kof tok_pfa tok_gtdb tok_nr tok_vir <<< "$line"
	sample="${sample:-}"
	tok_kof="${tok_kof:-$sample}"
	tok_pfa="${tok_pfa:-$sample}"
	tok_gtdb="${tok_gtdb:-$sample}"
	tok_nr="${tok_nr:-$sample}"
	tok_vir="${tok_vir:-$sample}"

	kof_file="$(choose_best_file "$KOF_DIR" "$tok_kof" "kofamscan_out")"
	pfa_file="$(choose_best_file "$PFAM_DIR" "$tok_pfa" "pfam_out")"

	# GTDB optional, look for contig2classification and gtdb_ prefix
	gtdb_file=""
	if [[ -n "$GTDB_DIR" ]]; then
		gtdb_file="$(choose_best_file "$GTDB_DIR" "$tok_gtdb" "contig2classification")"
		# prefer gtdb_ over nr_ if both exist
		if [[ -n "$gtdb_file" && "$(basename "$gtdb_file")" != gtdb_* ]]; then
		 cand="$(ls -1t "$GTDB_DIR" 2>/dev/null | grep -E "^gtdb_.*${tok_gtdb}.*contig2classification" | head -1 || true)"
		 [[ -n "$cand" ]] && gtdb_file="${GTDB_DIR}/${cand}"
		fi
	fi

	# NR optional
	nr_file=""
	if [[ -n "$NR_DIR" ]]; then
		nr_file="$(choose_best_file "$NR_DIR" "$tok_nr" "contig2classification")"
		if [[ -n "$nr_file" && "$(basename "$nr_file")" != nr_* ]]; then
		 cand="$(ls -1t "$NR_DIR" 2>/dev/null | grep -E "^nr_.*${tok_nr}.*contig2classification" | head -1 || true)"
		 [[ -n "$cand" ]] && nr_file="${NR_DIR}/${cand}"
		fi
	fi

	# Viral optional, summary tsv
	vir_file=""
	if [[ -n "$VIRUS_DIR" ]]; then
		vir_file="$(choose_best_file "$VIRUS_DIR" "$tok_vir" "virus_summary.tsv")"
	fi

	if [[ -z "$kof_file" || -z "$pfa_file" ]]; then
		if [[ -z "$kof_file" && -z "$pfa_file" ]]; then
			log "[SKIP] ${sample} no kofam or pfam match"
		elif [[ -z "$kof_file" ]]; then
			near="$(suggest_nearby "$KOF_DIR" "$tok_kof" "kofamscan_out")"
			log "[SKIP] ${sample} pfam ok, kofam missing, kofam_token=${tok_kof} nearby=[${near}]"
			miss_kof=$((miss_kof+1))
		else
			near="$(suggest_nearby "$PFAM_DIR" "$tok_pfa" "pfam_out")"
			log "[SKIP] ${sample} kofam ok, pfam missing, pfam_token=${tok_pfa} nearby=[${near}]"
			miss_pfa=$((miss_pfa+1))
		fi
		skipped=$((skipped+1))
		continue
	fi

	log "[RUN]  sample=${sample}"
	log "       kofam_file=${kof_file}"
	log "       pfam_file=${pfa_file}"
	[[ -n "$gtdb_file" ]] && log "       gtdb_file=${gtdb_file}" || log "       gtdb_file=NA"
	[[ -n "$nr_file"   ]] && log "       nr_file=${nr_file}"     || log "       nr_file=NA"
	[[ -n "$vir_file"  ]] && log "       virus_file=${vir_file}" || log "       virus_file=NA"
	log "       join=${JOIN_TYPE}"

	# build args
	args=( --kofam "$kof_file" --pfam "$pfa_file" --join "$JOIN_TYPE" --sample "$sample" --outdir "$base_outdir" )
	[[ -n "$gtdb_file" ]] && args+=( --cat_gtdb "$gtdb_file" )
	[[ -n "$nr_file"   ]] && args+=( --cat_nr   "$nr_file" )
	[[ -n "$vir_file"  ]] && args+=( --virus    "$vir_file" )

	# run python and tee output
	if python3 "$PY_SCRIPT" "${args[@]}" 2>&1 | tee -a "$LOG" >&2; then
		ok=$((ok+1))
	else
		log "[ERR]  sample=${sample} python failed"
		skipped=$((skipped+1))
	fi
done < <(sed -e 's/\r$//' "$SAMPLE_LIST")

log "[DONE] success=${ok} skipped=${skipped} missing_kof=${miss_kof} missing_pfam=${miss_pfa}"