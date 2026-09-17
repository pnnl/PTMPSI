#!/usr/bin/env bash
set -euo pipefail


OUTDIR="cp12_public_refs"

# Explicit RefSeq CP12 accessions
MANUAL_REFS=("YP_004323596.1" "WP_019478051.1")

# Length window for CP12 RefSeq title search
MIN_LEN=40
MAX_LEN=120

# UniProt query (current REST API syntax)
UNIPROT_QUERY='(xref:pfam-PF02672) AND (reviewed:true)'


mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

TAX_CACHE="taxid_lineage.cache.tsv"
touch "${TAX_CACHE}"

get_lineage () {
	local taxid="$1"

	# taxonomy DB accepts numeric TaxID only
	if [[ -z "${taxid}" || ! "${taxid}" =~ ^[0-9]+$ ]]; then
		echo "NA"
		return
	fi

	# check cache
	local cached
	cached=$(awk -F"\t" -v id="${taxid}" '$1==id{print $2}' "${TAX_CACHE}")
	if [[ -n "${cached}" ]]; then
		echo "${cached}"
		return
	fi

	# fetch lineage
	local lineage
	lineage=$(efetch -db taxonomy -id "${taxid}" -format xml 2>/dev/null \
		| grep -oPm1 "(?<=<Lineage>)[^<]+")

	if [[ -z "${lineage}" ]]; then
		lineage="NA"
	fi

	echo -e "${taxid}\t${lineage}" >> "${TAX_CACHE}"
	echo "${lineage}"
}

echo "[1/5] Fetching explicit RefSeq CP12 accessions"

NCBI_MANUAL_FASTA="ncbi_refseq_manual.fasta"
NCBI_MANUAL_META="ncbi_refseq_manual.metadata.tsv"

echo -e "Accession\tOrganism\tTaxID\tLineage\tProteinTitle\tSourceDB" \
	> "${NCBI_MANUAL_META}"

for acc in "${MANUAL_REFS[@]}"; do
	efetch -db protein -id "${acc}" -format fasta >> "${NCBI_MANUAL_FASTA}"

	efetch -db protein -id "${acc}" -format docsum \
	| xtract -pattern DocumentSummary \
		-element AccessionVersion Organism TaxId Title \
	| while IFS=$'\t' read -r accession organism taxid title; do
		lineage=$(get_lineage "${taxid}")
		printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
			"${accession}" \
			"${organism}" \
			"${taxid:-NA}" \
			"${lineage}" \
			"${title}" \
			"NCBI_refseq_manual"
	done >> "${NCBI_MANUAL_META}"
done

echo "[2/5] Searching NCBI RefSeq for CP12-named proteins"

NCBI_SEARCH_FASTA="ncbi_refseq_cp12_named.fasta"
NCBI_SEARCH_META="ncbi_refseq_cp12_named.metadata.tsv"

echo -e "Accession\tOrganism\tTaxID\tLineage\tProteinTitle\tSourceDB" \
	> "${NCBI_SEARCH_META}"

QUERY="CP12[Title] AND ${MIN_LEN}:${MAX_LEN}[SLEN] AND srcdb_refseq[PROP]"

esearch -db protein -query "${QUERY}" \
| efetch -format fasta \
> "${NCBI_SEARCH_FASTA}"

esearch -db protein -query "${QUERY}" \
| efetch -format docsum \
| xtract -pattern DocumentSummary \
	-element AccessionVersion Organism TaxId Title \
| while IFS=$'\t' read -r accession organism taxid title; do
	lineage=$(get_lineage "${taxid}")
	printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
		"${accession}" \
		"${organism}" \
		"${taxid:-NA}" \
		"${lineage}" \
		"${title}" \
		"NCBI_refseq_CP12"
done >> "${NCBI_SEARCH_META}"

echo "[3/5] Fetching UniProt reviewed CP12 references"

UNIPROT_FASTA="uniprot_cp12_refs.fasta"
UNIPROT_META_RAW="uniprot_cp12_metadata.raw.tsv"
UNIPROT_META="uniprot_cp12_metadata.tsv"

curl -sG https://rest.uniprot.org/uniprotkb/stream \
	--data-urlencode "query=${UNIPROT_QUERY}" \
	--data-urlencode "format=fasta" \
	> "${UNIPROT_FASTA}"

curl -sG https://rest.uniprot.org/uniprotkb/stream \
	--data-urlencode "query=${UNIPROT_QUERY}" \
	--data-urlencode "format=tsv" \
	--data-urlencode "fields=accession,protein_name,organism_name,organism_id,lineage" \
	> "${UNIPROT_META_RAW}"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {
	print "Accession","Organism","TaxID","Lineage","ProteinTitle","SourceDB"; next
}
{
	# normalize delimiter to NCBI style
	gsub(/, /,"; ",$5)

	# remove rank annotations, keep taxon names only
	gsub(/\s*\([^)]*\)/,"",$5)

	# normalize spacing
	gsub(/;\s*/,"; ",$5)
	gsub(/^\s+|\s+$/,"",$5)

	print $1,$3,$4,$5,$2,"UniProt"
}' "${UNIPROT_META_RAW}" > "${UNIPROT_META}"


echo "[4/5] Combining and cleaning FASTA"

ALL_FASTA_RAW="cp12_public_refs.raw.fasta"
ALL_FASTA_CLEAN="cp12_public_refs.clean.fasta"

cat \
"${NCBI_MANUAL_FASTA}" \
"${NCBI_SEARCH_FASTA}" \
"${UNIPROT_FASTA}" \
> "${ALL_FASTA_RAW}"

awk '
/^>/ {
	split($0,a," ")
	print a[1]
	next
}
{
	gsub(/\*/, "", $0)
	print
}
' "${ALL_FASTA_RAW}" > "${ALL_FASTA_CLEAN}"


echo "[5/5] Combining metadata tables"

ALL_META="cp12_public_refs.metadata.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==1 || FNR>1' \
"${NCBI_MANUAL_META}" \
"${NCBI_SEARCH_META}" \
"${UNIPROT_META}" \
> "${ALL_META}"

echo "----------------------------------------"
echo "CP12 public reference extraction complete"
echo "FASTA (clean): ${OUTDIR}/${ALL_FASTA_CLEAN}"
echo "Metadata:      ${OUTDIR}/${ALL_META}"
echo "----------------------------------------"