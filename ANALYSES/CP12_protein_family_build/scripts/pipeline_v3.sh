#!/usr/bin/env bash
set -euo pipefail

FULL_FASTA="master_full_sequences.fasta"
PARTIAL_FASTA="master_partial_sequences.fasta"
FULL_ANN="master_full_annotation.tsv"
PARTIAL_ANN="master_partial_annotation.tsv"

COMBINED_FASTA="combined_sequences.fasta"
COMBINED_ANN="combined_annotation.tsv"
DEREP_FASTA="combined_sequences_derep.fasta"
CLSTR_FILE="${DEREP_FASTA}.clstr"
DEREP_IDS="derep_ids.txt"
DEREP_ANN="combined_annotation_derep.tsv"
CP12_TSV="CP12_motif_and_classification.tsv"
FINAL_ANN="combined_annotation_derep_cp12.tsv"

echo "Starting CP12 pipeline"


echo "[1/6] Combining FASTA files"
cat "${FULL_FASTA}" "${PARTIAL_FASTA}" > "${COMBINED_FASTA}"

echo "[2/6] Merging annotation tables"

python << 'EOF'
import csv

full_ann = "master_full_annotation.tsv"
partial_ann = "master_partial_annotation.tsv"
output_ann = "combined_annotation.tsv"

with open(full_ann) as f:
	reader = csv.DictReader(f, delimiter="\t")
	header = reader.fieldnames
	full_rows = list(reader)

with open(partial_ann) as f:
	reader = csv.DictReader(f, delimiter="\t")
	partial_header = reader.fieldnames
	partial_rows = list(reader)

if header != partial_header:
	raise SystemExit("ERROR: Annotation headers do not match")

with open(output_ann, "w", newline="") as out:
	writer = csv.DictWriter(out, delimiter="\t",
		fieldnames=header + ["SequenceStatus"])
	writer.writeheader()

	for row in full_rows:
		row["SequenceStatus"] = "full"
		writer.writerow(row)

	for row in partial_rows:
		row["SequenceStatus"] = "partial"
		writer.writerow(row)
EOF

echo "[3/6] Running CD-HIT dereplication"

cd-hit \
	-i "${COMBINED_FASTA}" \
	-o "${DEREP_FASTA}" \
	-c 1.0 \
	-n 5 \
	-d 0 \
	-T 8 \
	-M 16000


echo "[4/6] Extracting dereplicated IDs"

grep "^>" "${DEREP_FASTA}" \
	| sed 's/^>//' \
	| awk '{print $1}' \
	> "${DEREP_IDS}"


echo "[5/6] Filtering annotation table to dereplicated proteins"

python << 'EOF'
import csv

annotation_file = "combined_annotation.tsv"
id_file = "derep_ids.txt"
output_file = "combined_annotation_derep.tsv"

with open(id_file) as f:
	keep_ids = set(line.strip() for line in f)

with open(annotation_file) as f:
	reader = csv.DictReader(f, delimiter="\t")
	fieldnames = reader.fieldnames
	rows = list(reader)

if "new_gene_name" not in fieldnames:
	raise SystemExit("ERROR: new_gene_name column not found")

with open(output_file, "w", newline="") as out:
	writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
	writer.writeheader()

	for row in rows:
		if row["new_gene_name"] in keep_ids:
			writer.writerow(row)
EOF


echo "[6/6] Running CP12 motif classification"

python << 'EOF'
from Bio import SeqIO
from itertools import combinations

input_fasta = "combined_sequences_derep.fasta"
output_tsv = "CP12_motif_and_classification.tsv"

def pairwise_cys_spacings(cys_positions):
	return [j - i for i, j in combinations(cys_positions, 2)]

with open(output_tsv, "w") as out:
	out.write(
		"ProteinID\tLength\tNumCys\tCysPositions\tCysSpacings\t"
		"AcidicFraction\tCP12_Class\tCP12_Reason\n"
	)

	for record in SeqIO.parse(input_fasta, "fasta"):
		protein_id = record.id.split()[0]
		seq = str(record.seq)

		length = len(seq)
		cys_positions = [i + 1 for i, aa in enumerate(seq) if aa == "C"]
		num_cys = len(cys_positions)
		acidic_fraction = (seq.count("E") + seq.count("D")) / length
		cys_spacings = pairwise_cys_spacings(cys_positions)

		if (
			length < 50 or length > 250 or
			num_cys < 2 or
			acidic_fraction < 0.15
		):
			cp12_class = "Likely_false_positive"
			reason = "Fails CP12 biochemical expectations"

		elif (
			60 <= length <= 150 and
			num_cys == 2 and
			acidic_fraction >= 0.20 and
			len(cys_spacings) == 1 and
			15 <= cys_spacings[0] <= 80
		):
			cp12_class = "Canonical_CP12"
			reason = "Matches bacterial-like CP12"

		else:
			cp12_class = "Divergent_CP12"
			reason = "CP12-like features"

		out.write(
			f"{protein_id}\t{length}\t{num_cys}\t{cys_positions}\t"
			f"{cys_spacings}\t{round(acidic_fraction,3)}\t"
			f"{cp12_class}\t{reason}\n"
		)
EOF

echo "[FINAL] Appending CP12 annotations"

python << 'EOF'
import csv

ann_file = "combined_annotation_derep.tsv"
cp12_file = "CP12_motif_and_classification.tsv"
output_file = "combined_annotation_derep_cp12.tsv"

cp12 = {}
with open(cp12_file) as f:
	reader = csv.DictReader(f, delimiter="\t")
	for row in reader:
		cp12[row["ProteinID"]] = row

with open(ann_file) as f:
	reader = csv.DictReader(f, delimiter="\t")
	headers = reader.fieldnames
	rows = list(reader)

cp12_map = {
	"CP12_Class": "CP12_Class",
	"CP12_Reason": "CP12_Reason",
	"CP12_Length": "Length",
	"CP12_NumCys": "NumCys",
	"CP12_CysPositions": "CysPositions",
	"CP12_CysSpacings": "CysSpacings",
	"CP12_AcidicFraction": "AcidicFraction"
}

with open(output_file, "w", newline="") as out:
	writer = csv.writer(out, delimiter="\t")
	writer.writerow(headers + list(cp12_map.keys()))

	for row in rows:
		pid = row["new_gene_name"]
		cp12_row = cp12.get(pid, {})
		writer.writerow(
			[row[h] for h in headers] +
			[cp12_row.get(v, "NA") for v in cp12_map.values()]
		)
EOF

echo "Pipeline completed successfully"