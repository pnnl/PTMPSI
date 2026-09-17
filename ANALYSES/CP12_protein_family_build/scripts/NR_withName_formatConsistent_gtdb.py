#!/usr/bin/env python3
import sys
from pathlib import Path

RANK2PREFIX = {
	"superkingdom": "d__",
	"phylum": "p__",
	"class": "c__",
	"order": "o__",
	"family": "f__",
	"genus": "g__",
	"species": "s__",
}

def to_gtdb_style_lineage(full_lineage_fields):
	# full_lineage_fields are cols 6.., each like "Name (rank): score"
	parts = ["root"]
	seen = set()
	for fld in full_lineage_fields:
		fld = fld.strip()
		if not fld or ":" not in fld:
			continue
		left = fld.split(":", 1)[0].strip()  # "Name (rank)"
		# extract name and rank
		if left.endswith(")") and " (" in left:
			name, rank = left.rsplit(" (", 1)
			rank = rank[:-1]  # drop trailing ")"
			name = name.strip()
			rank = rank.strip().lower()
		else:
			# unexpected format, skip
			continue
		# keep root explicitly
		if name.lower() == "root":
			if "root" not in seen:
				parts = ["root"]
				seen.add("root")
			continue
		# map ranks to GTDB-style prefixes
		prefix = RANK2PREFIX.get(rank)
		if prefix:
			key = (rank, name)
			if key not in seen:
				parts.append(f"{prefix}{name}")
				seen.add(key)
		# ignore "no rank", "cellular organisms", and unmapped ranks
	return ";".join(parts)

def convert_one(in_path: Path):
	out_path = Path(str(in_path).replace("_withNames.txt", "_withNames_gtdbStyle.txt"))
	with in_path.open("r", encoding="utf-8") as fin, out_path.open("w", encoding="utf-8") as fout:
		for line in fin:
			line = line.rstrip("\n")
			if not line:
				continue
			if line.startswith("# contig"):
				# Normalize header to match GTDB 5-column header
				fout.write("# contig\tclassification\treason\tlineage\tlineage scores (f: 0.5)\n")
				continue
			fields = line.split("\t")
			if len(fields) < 5:
				continue  # malformed line
			contig, cls, reason, lineage_ids, lineage_scores = fields[:5]
			full_lineage_fields = fields[5:]  # names with ranks and scores
			lineage_gtdb_style = to_gtdb_style_lineage(full_lineage_fields)
			# Write 5 columns to mirror GTDB file
			fout.write("\t".join([contig, cls, reason, lineage_gtdb_style, lineage_scores]) + "\n")
	print(f"[OK] {in_path.name} -> {out_path.name}")

def main():
	if len(sys.argv) < 2:
		print("Usage: python nr_withNames_to_gtdb_style.py <nr_*withNames.txt> [more files...]", file=sys.stderr)
		sys.exit(2)
	for arg in sys.argv[1:]:
		convert_one(Path(arg))

if __name__ == "__main__":
	letters = main()