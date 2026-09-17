#!/usr/bin/env python3
"""
Build cp12_combined.metadata_v3.1.txt and cp12_combined_pf02672.sto automatically.

Inputs expected from Step 11 workflow:
- Metagenome annotation: combined_annotation_derep_cp12.tsv (from pipeline_v3.sh)
- Metagenome FASTA: combined_sequences_derep.fasta (from pipeline_v3.sh)
- Public reference metadata: cp12_public_refs.metadata.tsv (from public_seq_fetching_v6.sh)
- Public reference FASTA: cp12_public_refs.clean.fasta (from public_seq_fetching_v6.sh)
- PF02672 HMM: PF02672.hmm

Outputs:
- cp12_combined.metadata_v3.1.txt
- cp12_combined_combined.faa
- cp12_combined_pf02672.sto

Notes:
- Metadata columns are harmonized to the v3.1 schema used in downstream scripts.
- Use --exclude-accessions to apply one-off curated removals.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


V31_COLUMNS = [
	"Accession",
	"SourceType",
	"SourceDB",
	"Organism",
	"TaxID",
	"Lineage_ncbi",
	"Lineage_consensus",
	"pfamsearch:E-value",
	"pfamsearch:query_name",
	"pfamsearch:score",
	"pfamsearch:target_accession",
	"pfamsearch:target_name",
	"SequenceStatus",
	"Length",
	"NumCys",
	"CysPositions",
	"CysSpacings",
	"AcidicFraction",
	"CP12_Class",
	"CP12_ClassReason",
	"Lineage_consensus_cleaned",
]


@dataclass
class SeqRec:
	seq_id: str
	header: str
	seq: str


def parse_args() -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Build cp12_combined metadata + sto automatically")
	p.add_argument("--metagenome-annotation", required=True, help="combined_annotation_derep_cp12.tsv")
	p.add_argument("--metagenome-fasta", required=True, help="combined_sequences_derep.fasta")
	p.add_argument("--reference-metadata", required=True, help="cp12_public_refs.metadata.tsv")
	p.add_argument("--reference-fasta", required=True, help="cp12_public_refs.clean.fasta")
	p.add_argument("--hmm", required=True, help="PF02672.hmm")
	p.add_argument("--outdir", required=True, help="Output directory")
	p.add_argument("--metadata-name", default="cp12_combined.metadata_v3.1.txt")
	p.add_argument("--combined-fasta-name", default="cp12_combined_combined.faa")
	p.add_argument("--sto-name", default="cp12_combined_pf02672.sto")
	p.add_argument(
		"--exclude-accessions",
		default=None,
		help="Optional text file with one accession per line to drop (curated exclusions)",
	)
	p.add_argument(
		"--strict-pfam-cp12",
		action="store_true",
		help="Keep metagenome rows only when pfamsearch:target_accession starts with PF02672",
	)
	return p.parse_args()


def read_fasta(path: Path) -> Dict[str, SeqRec]:
	out: Dict[str, SeqRec] = {}
	header: Optional[str] = None
	chunks: List[str] = []

	with path.open("r", encoding="utf-8") as fh:
		for raw in fh:
			line = raw.rstrip("\n")
			if not line:
				continue
			if line.startswith(">"):
				if header is not None:
					sid = header.split()[0]
					out[sid] = SeqRec(seq_id=sid, header=header, seq="".join(chunks))
				header = line[1:].strip()
				chunks = []
			else:
				chunks.append(line.strip())

	if header is not None:
		sid = header.split()[0]
		out[sid] = SeqRec(seq_id=sid, header=header, seq="".join(chunks))

	return out


def calc_cys_positions(seq: str) -> List[int]:
	return [i + 1 for i, aa in enumerate(seq) if aa == "C"]


def pairwise_spacings(pos: Sequence[int]) -> List[int]:
	out: List[int] = []
	for i in range(len(pos)):
		for j in range(i + 1, len(pos)):
			out.append(pos[j] - pos[i])
	return out


def classify_cp12(seq: str) -> Tuple[str, str]:
	length = len(seq)
	cys = calc_cys_positions(seq)
	n_cys = len(cys)
	acidic = ((seq.count("E") + seq.count("D")) / length) if length else 0.0
	spacings = pairwise_spacings(cys)

	if length < 50 or length > 250 or n_cys < 2 or acidic < 0.15:
		return "Likely_false_positive", "Fails CP12 biochemical expectations"

	if 60 <= length <= 150 and n_cys == 2 and acidic >= 0.20 and len(spacings) == 1 and 15 <= spacings[0] <= 80:
		return "Canonical_CP12", "Matches bacterial-like CP12"

	return "Divergent_CP12", "CP12-like features"


def normalize_taxid(v: str) -> str:
	s = (v or "").strip()
	return s if s else "NA"


def pick_first_nonempty(row: Dict[str, str], keys: Sequence[str], fallback: str = "NA") -> str:
	for k in keys:
		v = (row.get(k, "") or "").strip()
		if v:
			return v
	return fallback


def parse_domain_phylum(lineage: str) -> Tuple[str, str]:
	s = (lineage or "").strip()
	if not s or s == "NA":
		return ("Unclassified", "Unclassified")

	toks = [t.strip() for t in s.split(";") if t.strip()]

	domain = "Unclassified"
	phylum = "Unclassified"

	for t in toks:
		if t.startswith("d__"):
			domain = t.replace("d__", "", 1) or "Unclassified"
			break

	for t in toks:
		if t.startswith("p__"):
			phylum = t.replace("p__", "", 1) or "Unclassified"
			break

	if domain == "Unclassified":
		low = s.lower()
		if "viruses" in low or "virus" in low:
			domain = "Viruses"
		elif "archaea" in low:
			domain = "Archaea"
		elif "bacteria" in low:
			domain = "Bacteria"
		elif "fungi" in low:
			domain = "Fungi"
		elif "viridiplantae" in low or "plantae" in low:
			domain = "Viridiplantae"
		elif "metazoa" in low or "animalia" in low:
			domain = "Metazoa"
		elif re.search(r"(^|[;\s])sar([;\s]|$)", low):
			domain = "SAR"

	if phylum == "Unclassified":
		# For NCBI-style lineages, choose first token after domain that looks like a phylum-level name.
		for i, t in enumerate(toks):
			plain = t.replace("d__", "").strip()
			if plain.lower() == domain.lower() and i + 1 < len(toks):
				cands = toks[i + 1 :]
				for c in cands:
					c2 = c.replace("p__", "").strip()
					if re.search(r"(ota|phyta|mycota|zoa)$", c2, flags=re.IGNORECASE):
						phylum = c2
						break
				if phylum == "Unclassified" and cands:
					phylum = cands[min(1, len(cands) - 1)].replace("p__", "").strip()
				break

	return (domain or "Unclassified", phylum or "Unclassified")


def lineage_cleaned(lineage: str) -> str:
	d, p = parse_domain_phylum(lineage)
	return f"{d}:{p}"


def load_exclude_ids(path: Optional[Path]) -> Set[str]:
	if path is None:
		return set()
	out: Set[str] = set()
	with path.open("r", encoding="utf-8") as fh:
		for line in fh:
			v = line.strip()
			if v and not v.startswith("#"):
				out.add(v)
	return out


def format_float(x: float, nd: int = 3) -> str:
	if math.isnan(x):
		return "NA"
	return f"{x:.{nd}f}".rstrip("0").rstrip(".")


def build_metagenome_rows(
	ann_path: Path,
	metagenome_fasta: Dict[str, SeqRec],
	exclude_ids: Set[str],
	strict_pfam_cp12: bool,
) -> List[Dict[str, str]]:
	rows: List[Dict[str, str]] = []

	with ann_path.open("r", encoding="utf-8", newline="") as fh:
		rdr = csv.DictReader(fh, delimiter="\t")

		for r in rdr:
			acc = pick_first_nonempty(r, ["new_gene_name", "Accession"]) 
			if acc == "NA" or acc in exclude_ids:
				continue
			if acc not in metagenome_fasta:
				continue

			if strict_pfam_cp12:
				targ = pick_first_nonempty(r, ["pfamsearch:target_accession"], fallback="")
				if not targ.startswith("PF02672"):
					continue

			seq = metagenome_fasta[acc].seq
			cys_pos = calc_cys_positions(seq)
			spacings = pairwise_spacings(cys_pos)
			acidic = ((seq.count("E") + seq.count("D")) / len(seq)) if seq else float("nan")

			cp12_class = pick_first_nonempty(r, ["CP12_Class"], fallback="NA")
			cp12_reason = pick_first_nonempty(r, ["CP12_Reason", "CP12_ClassReason"], fallback="NA")
			if cp12_class == "NA":
				cp12_class, cp12_reason = classify_cp12(seq)

			lineage_ncbi = pick_first_nonempty(r, ["cat_nr:lineage", "Lineage_ncbi"], fallback="root")
			lineage_cons = pick_first_nonempty(
				r,
				["MatchBasis-greedy", "MatchBasis-conservative", "cat_gtdb:lineage", "Lineage_consensus"],
				fallback=lineage_ncbi,
			)

			out = {
				"Accession": acc,
				"SourceType": "Metagenome",
				"SourceDB": "Metagenome",
				"Organism": "NA",
				"TaxID": "NA",
				"Lineage_ncbi": lineage_ncbi,
				"Lineage_consensus": lineage_cons,
				"pfamsearch:E-value": pick_first_nonempty(r, ["pfamsearch:E-value"]),
				"pfamsearch:query_name": pick_first_nonempty(r, ["pfamsearch:query_name"]),
				"pfamsearch:score": pick_first_nonempty(r, ["pfamsearch:score"]),
				"pfamsearch:target_accession": pick_first_nonempty(r, ["pfamsearch:target_accession"]),
				"pfamsearch:target_name": pick_first_nonempty(r, ["pfamsearch:target_name"]),
				"SequenceStatus": pick_first_nonempty(r, ["SequenceStatus"], fallback="NA"),
				"Length": str(len(seq)),
				"NumCys": str(len(cys_pos)),
				"CysPositions": str(cys_pos),
				"CysSpacings": str(spacings),
				"AcidicFraction": format_float(acidic),
				"CP12_Class": cp12_class,
				"CP12_ClassReason": cp12_reason,
				"Lineage_consensus_cleaned": lineage_cleaned(lineage_cons),
			}
			rows.append(out)

	return rows


def build_reference_rows(
	ref_meta_path: Path,
	ref_fasta: Dict[str, SeqRec],
	exclude_ids: Set[str],
) -> List[Dict[str, str]]:
	rows: List[Dict[str, str]] = []

	with ref_meta_path.open("r", encoding="utf-8", newline="") as fh:
		rdr = csv.DictReader(fh, delimiter="\t")

		for r in rdr:
			acc = (r.get("Accession", "") or "").strip()
			if not acc or acc in exclude_ids:
				continue
			if acc not in ref_fasta:
				continue

			seq = ref_fasta[acc].seq
			cys_pos = calc_cys_positions(seq)
			spacings = pairwise_spacings(cys_pos)
			acidic = ((seq.count("E") + seq.count("D")) / len(seq)) if seq else float("nan")
			cp12_class, cp12_reason = classify_cp12(seq)

			lineage = (r.get("Lineage", "") or "").strip() or "NA"

			out = {
				"Accession": acc,
				"SourceType": "Reference",
				"SourceDB": (r.get("SourceDB", "") or "").strip() or "Reference",
				"Organism": (r.get("Organism", "") or "").strip() or "NA",
				"TaxID": normalize_taxid((r.get("TaxID", "") or "").strip()),
				"Lineage_ncbi": lineage,
				"Lineage_consensus": lineage,
				"pfamsearch:E-value": "NA",
				"pfamsearch:query_name": "NA",
				"pfamsearch:score": "NA",
				"pfamsearch:target_accession": "NA",
				"pfamsearch:target_name": "NA",
				"SequenceStatus": "NA",
				"Length": str(len(seq)),
				"NumCys": str(len(cys_pos)),
				"CysPositions": str(cys_pos),
				"CysSpacings": str(spacings),
				"AcidicFraction": format_float(acidic),
				"CP12_Class": cp12_class,
				"CP12_ClassReason": cp12_reason,
				"Lineage_consensus_cleaned": lineage_cleaned(lineage),
			}
			rows.append(out)

	return rows


def dedupe_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
	# Keep first occurrence by accession for deterministic behavior.
	seen: Set[str] = set()
	out: List[Dict[str, str]] = []
	for r in rows:
		acc = r["Accession"]
		if acc in seen:
			continue
		seen.add(acc)
		out.append(r)
	return out


def write_metadata(rows: List[Dict[str, str]], out_path: Path) -> None:
	with out_path.open("w", encoding="utf-8", newline="") as fh:
		w = csv.DictWriter(fh, delimiter="\t", fieldnames=V31_COLUMNS, extrasaction="ignore")
		w.writeheader()
		for r in rows:
			w.writerow(r)


def write_combined_fasta(rows: List[Dict[str, str]], seq_map: Dict[str, SeqRec], out_path: Path) -> None:
	with out_path.open("w", encoding="utf-8") as fh:
		for r in rows:
			acc = r["Accession"]
			rec = seq_map.get(acc)
			if rec is None:
				continue
			fh.write(f">{rec.header}\n")
			s = rec.seq
			for i in range(0, len(s), 60):
				fh.write(s[i : i + 60] + "\n")


def run_hmmalign(hmm: Path, fasta: Path, out_sto: Path) -> None:
	cmd = ["hmmalign", str(hmm), str(fasta)]
	with out_sto.open("w", encoding="utf-8") as out_fh:
		proc = subprocess.run(cmd, stdout=out_fh, stderr=subprocess.PIPE, text=True)
	if proc.returncode != 0:
		raise RuntimeError(f"hmmalign failed: {proc.stderr.strip()}")


def main() -> int:
	args = parse_args()

	outdir = Path(args.outdir).resolve()
	outdir.mkdir(parents=True, exist_ok=True)

	metagenome_ann = Path(args.metagenome_annotation).resolve()
	metagenome_fa = Path(args.metagenome_fasta).resolve()
	ref_meta = Path(args.reference_metadata).resolve()
	ref_fa = Path(args.reference_fasta).resolve()
	hmm = Path(args.hmm).resolve()

	for p in [metagenome_ann, metagenome_fa, ref_meta, ref_fa, hmm]:
		if not p.exists():
			print(f"[ERR] Missing required input: {p}", file=sys.stderr)
			return 2

	exclude_ids = load_exclude_ids(Path(args.exclude_accessions).resolve() if args.exclude_accessions else None)

	metagenome_seqs = read_fasta(metagenome_fa)
	reference_seqs = read_fasta(ref_fa)
	all_seq_map = {**reference_seqs, **metagenome_seqs}

	mg_rows = build_metagenome_rows(
		ann_path=metagenome_ann,
		metagenome_fasta=metagenome_seqs,
		exclude_ids=exclude_ids,
		strict_pfam_cp12=args.strict_pfam_cp12,
	)
	ref_rows = build_reference_rows(ref_meta, reference_seqs, exclude_ids)

	rows = dedupe_rows(mg_rows + ref_rows)

	meta_out = outdir / args.metadata_name
	fasta_out = outdir / args.combined_fasta_name
	sto_out = outdir / args.sto_name

	write_metadata(rows, meta_out)
	write_combined_fasta(rows, all_seq_map, fasta_out)
	run_hmmalign(hmm, fasta_out, sto_out)

	print(f"[OK] metadata rows: {len(rows)}")
	print(f"[OK] wrote metadata: {meta_out}")
	print(f"[OK] wrote combined fasta: {fasta_out}")
	print(f"[OK] wrote stockholm: {sto_out}")

	return 0


if __name__ == "__main__":
	raise SystemExit(main())