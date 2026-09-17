#!/usr/bin/env bash
set -euo pipefail

# ---------- defaults, edit if needed ----------
PY_EXTRACT="scripts/extract_by_id_split_partial_full_v2.py"
# ---------------------------------------------

usage(){
	cat <<'EOF' 1>&2
Usage:
	batch_extract_fasta_by_new_gene_v6.1.sh <filter_outdir> <fasta_dir> [subdirs]

Arguments
	filter_outdir  root that contains func_filter, taxa_filter, func_then_taxa
	fasta_dir      directory that holds the full protein FASTA files
	subdirs        optional space separated list, default "func_filter taxa_filter func_then_taxa"

Behavior
	1) For each subdir, run the Python extractor per sample table to produce:
				*_full_sequences.fasta
				*_partial_sequences.fasta
				*_full_sequences.txt
				*_partial_sequences.txt
	2) Build MASTER strictly by scanning those four sidecar patterns.
			MASTER includes:
				master_full_sequences.fasta
				master_partial_sequences.fasta
				master_full_annotation.tsv
				master_partial_annotation.tsv
				master_counts_full_greedy.tsv
				master_counts_full_conservative.tsv
				master_counts_partial_greedy.tsv
				master_counts_partial_conservative.tsv
				summary.log
EOF
	exit 1
}

[[ $# -ge 2 && $# -le 3 ]] || usage

OUTROOT="$1"
FASTA_DIR="$2"
SUBDIRS="${3:-func_filter taxa_filter func_then_taxa}"

[[ -d "$OUTROOT"   ]] || { echo "filter_outdir not found, $OUTROOT" >&2; exit 1; }
[[ -d "$FASTA_DIR" ]] || { echo "fasta_dir not found, $FASTA_DIR" >&2; exit 1; }
[[ -f "$PY_EXTRACT" ]] || { echo "extract script not found, $PY_EXTRACT" >&2; exit 1; }

ts(){ date +%F_%T; }
log(){ printf "%s %s\n" "$(ts)" "$*" >&2; }
escape_ere(){ printf '%s' "$1" | sed -e 's/[][^$.*/+?(){}/|]/\\&/g'; }

choose_fasta(){
	dir="$1"; sample="$2"
	for ext in ".fasta" ".fasta.gz"; do
		exact="${dir}/renamed_prodigal_out_${sample}_megahit_contigs_filtered_1kb_prot${ext}"
		[[ -f "$exact" ]] && { printf '%s\n' "$exact"; return; }
	done
	token_ere="$(escape_ere "$sample")"
	while IFS= read -r f; do
		[[ -z "$f" ]] && continue
		base="${f##*/}"
		case "$base" in
			*.fasta|*.fasta.gz)
				if printf '%s' "$base" | grep -Eq "$token_ere"; then
					printf '%s\n' "$f"; return
				fi
				;;
		esac
	done < <(ls -1t "$dir" 2>/dev/null | awk -v d="$dir" '{print d "/" $0}')
	lpat="$(printf '%s' "$token_ere" | tr '[:upper:]' '[:lower:]')"
	while IFS= read -r f; do
		[[ -z "$f" ]] && continue
		base="${f##*/}"
		case "$base" in
			*.fasta|*.fasta.gz)
				if printf '%s' "$base" | tr '[:upper:]' '[:lower:]' | grep -Eq "$lpat"; then
					printf '%s\n' "$f"; return
				fi
				;;
		esac
	done < <(ls -1t "$dir" 2>/dev/null | awk -v d="$dir" '{print d "/" $0}')
	printf ""
}

append_tsv_body(){
	src="$1"; dst="$2"
	[[ ! -s "$src" ]] && return 0
	if [[ ! -s "$dst" ]]; then
		cat "$src" > "$dst"
	else
		awk 'NR>1' "$src" >> "$dst"
	fi
}

count_by_lineage_domain(){
	local src="$1"; local lineage_col="$2"; local dst="$3"
	awk -F'\t' -v lin="$lineage_col" '
		BEGIN{ci=-1; li=-1;}
		NR==1{
			for(i=1;i<=NF;i++){
				if($i=="Category") ci=i;
				if($i==lin) li=i;
			}
			next
		}
		{
			val = (li>0 ? $li : "");
			dom = "Unclassified";
			if(val ~ /d__Archaea/) dom="Archaea";
			else if(val ~ /d__Bacteria/) dom="Bacteria";
			else if(val ~ /(Viridiplantae)/) dom="Viridiplantae";
			else if(val ~ /(Metazoa)/) dom="Metazoa";
			else if(val ~ /(^|;)[^;]*SAR[^;]*/) dom="SAR";
			else if(val ~ /(Viruses|virus|Realm|kingdom__Viruses)/) dom="Viruses";
			else if(ci>0){
				if($ci ~ /Archaea/) dom="Archaea";
				else if($ci ~ /Bacteria/) dom="Bacteria";
				else if($ci ~ /Viridiplantae/) dom="Viridiplantae";
				else if($ci ~ /Metazoa/) dom="Metazoa";
				else if($ci ~ /SAR/) dom="SAR";
				else if($ci ~ /Viruses/) dom="Viruses";
			}
			cnt[dom]++;
		}
		END{
			order[1]="Archaea"; order[2]="Bacteria"; order[3]="SAR"; order[4]="Viridiplantae"; order[5]="Metazoa"; order[6]="Viruses"; order[7]="Unclassified";
			printf "Domain\tCount\n" > "'"$dst"'";
			for(i=1;i<=7;i++){k=order[i]; printf "%s\t%d\n", k, (cnt[k]?cnt[k]:0) >> "'"$dst"'";}
		}
	' "$src"
}

total=0; ok=0; skipped=0; missing_fa=0

for sd in $SUBDIRS; do
	subdir="$OUTROOT/$sd"
	[[ -d "$subdir" ]] || { log "[NOTE] skip, subdir not found [$sd]"; continue; }
	log "[SCAN] $subdir"

	# extractor pass, unchanged behavior
	while IFS= read -r -d '' tbl; do
		total=$((total+1))
		base="$(basename "$tbl")"
		sample="${base%%-filtered_*}"

		fa="$(choose_fasta "$FASTA_DIR" "$sample")"
		if [[ -z "$fa" ]]; then
			log "[SKIP] sample=$sample, fasta missing in $FASTA_DIR"
			missing_fa=$((missing_fa+1))
			continue
		fi

		log "[RUN ] sample=$sample, tbl=$(basename "$tbl"), fa=$(basename "$fa")"
		if python3 "$PY_EXTRACT" --filtered "$tbl" --fasta "$fa" --write-missing; then
			ok=$((ok+1))
		else
			log "[ERR ] extractor failed for $tbl"
			skipped=$((skipped+1))
		fi
	done < <(
		find "$subdir" -type f \( -name "*-filtered_*.txt" -o -name "*-filtered_*.tsv" \) \
			! -name "*.matched_headers.txt" \
			! -name "*_full_sequences.txt" \
			! -name "*_partial_sequences.txt" \
			! -name "*.missing_ids.txt" \
			! -path "$subdir/MASTER/*" \
			-print0 | sort -z
	)

	# MASTER assembly, excludes MASTER from inputs to avoid self-cat
	MASTER="$subdir/MASTER"
	mkdir -p "$MASTER"

	M_FULL_FASTA="$MASTER/master_full_sequences.fasta"
	M_PART_FASTA="$MASTER/master_partial_sequences.fasta"
	M_FULL_ANN="$MASTER/master_full_annotation.tsv"
	M_PART_ANN="$MASTER/master_partial_annotation.tsv"
	C_FULL_G="$MASTER/master_counts_full_greedy.tsv"
	C_FULL_C="$MASTER/master_counts_full_conservative.tsv"
	C_PART_G="$MASTER/master_counts_partial_greedy.tsv"
	C_PART_C="$MASTER/master_counts_partial_conservative.tsv"
	SUMMARY="$MASTER/summary.log"

	: > "$M_FULL_FASTA"
	: > "$M_PART_FASTA"
	rm -f "$M_FULL_ANN" "$M_PART_ANN" "$C_FULL_G" "$C_FULL_C" "$C_PART_G" "$C_PART_C"
	: > "$SUMMARY"

	# concat full FASTAs
	while IFS= read -r -d '' f; do
		cat "$f" >> "$M_FULL_FASTA"
	done < <(find "$subdir" -type f -name "*_full_sequences.fasta" ! -path "$subdir/MASTER/*" -print0 | sort -z)

	# concat partial FASTAs
	while IFS= read -r -d '' f; do
		cat "$f" >> "$M_PART_FASTA"
	done < <(find "$subdir" -type f -name "*_partial_sequences.fasta" ! -path "$subdir/MASTER/*" -print0 | sort -z)

	# concat full annotations with one header
	while IFS= read -r -d '' t; do
		append_tsv_body "$t" "$M_FULL_ANN"
	done < <(find "$subdir" -type f -name "*_full_sequences.txt" ! -path "$subdir/MASTER/*" -print0 | sort -z)

	# concat partial annotations with one header
	while IFS= read -r -d '' t; do
		append_tsv_body "$t" "$M_PART_ANN"
	done < <(find "$subdir" -type f -name "*_partial_sequences.txt" ! -path "$subdir/MASTER/*" -print0 | sort -z)

	# counts
	if [[ -s "$M_FULL_ANN" ]]; then
		count_by_lineage_domain "$M_FULL_ANN" "MatchBasis-greedy" "$C_FULL_G"
		count_by_lineage_domain "$M_FULL_ANN" "MatchBasis-conservative" "$C_FULL_C"
	fi
	if [[ -s "$M_PART_ANN" ]]; then
		count_by_lineage_domain "$M_PART_ANN" "MatchBasis-greedy" "$C_PART_G"
		count_by_lineage_domain "$M_PART_ANN" "MatchBasis-conservative" "$C_PART_C"
	fi

	full_n=$(grep -c '^>' "$M_FULL_FASTA" || true)
	part_n=$(grep -c '^>' "$M_PART_FASTA" || true)
	full_ann_n=$([[ -s "$M_FULL_ANN" ]] && awk 'END{print NR-1}' "$M_FULL_ANN" || echo 0)
	part_ann_n=$([[ -s "$M_PART_ANN" ]] && awk 'END{print NR-1}' "$M_PART_ANN" || echo 0)

	{
		printf "[kofam; pfam]:\n"
		printf "a total of %d samples:\n" "$total"
		printf "after extraction and aggregation:\n"
		printf "\t%4d samples failed\n" "$((skipped+missing_fa))"
		printf "\t%4d samples passed with sidecar fasta and annotation files\n\n" "$ok"

		printf "Full sequences: %d\n" "${full_n:-0}"
		printf "Partial sequences: %d\n" "${part_n:-0}"
		printf "Full annotation rows: %d\n" "$full_ann_n"
		printf "Partial annotation rows: %d\n" "$part_ann_n"
		printf "masters: %s %s %s %s\n" \
			"$(basename "$M_FULL_FASTA")" \
			"$(basename "$M_PART_FASTA")" \
			"$(basename "$M_FULL_ANN")" \
			"$(basename "$M_PART_ANN")"
		printf "counts: %s %s %s %s\n" \
			"$(basename "$C_FULL_G")" \
			"$(basename "$C_FULL_C")" \
			"$(basename "$C_PART_G")" \
			"$(basename "$C_PART_C")"
	} | tee -a "$SUMMARY"

	log "[DONE] $(basename "$subdir") MASTER written"
done

log "[ALL DONE] total=$total ok=$ok skipped=$skipped missing_fasta=$missing_fa"