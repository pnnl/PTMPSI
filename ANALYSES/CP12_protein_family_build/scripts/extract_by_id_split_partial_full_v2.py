#!/usr/bin/env python3

"""
Extract sequences by new_gene_name from a protein FASTA.

Inputs
	--filtered, path to a filtered TSV that has column 'new_gene_name'
	--fasta, path to the full protein FASTA for the same sample
	--outdir, optional, where to write outputs, default same directory as --filtered
	--id-col, optional, column name for IDs, default 'new_gene_name'
	--write-missing, optional flag, also write a .missing_ids.txt report

Output
	<filtered_basename>.fasta
	<filtered_basename>.matched_headers.txt
	<filtered_basename>_full_sequences.fasta
	<filtered_basename>_partial_sequences.fasta
	<filtered_basename>_full_sequences.txt      [NEW] annotations for full
	<filtered_basename>_partial_sequences.txt   [NEW] annotations for partial
	The header lines are preserved exactly as in the source FASTA.
"""

import argparse
import csv
import gzip
import os
import sys
from typing import Dict, Iterable, Iterator, List, Set, Tuple

NA = "NA"

def parse_args():
	parser = argparse.ArgumentParser(description="Extract sequences by new_gene_name into a new FASTA.")
	parser.add_argument("--filtered", required=True, help="Filtered TSV with column 'new_gene_name'")
	parser.add_argument("--fasta", required=True, help="Protein FASTA to pull sequences from")
	parser.add_argument("--outdir", default=None, help="Output directory, default same as --filtered")
	parser.add_argument("--id-col", default="new_gene_name", help="ID column to use, default new_gene_name")
	parser.add_argument("--write-missing", action="store_true", help="Write a sidecar list of IDs not found")
	return parser.parse_args()

def open_text(path: str):
	if path.endswith(".gz"):
		return gzip.open(path, "rt")
	return open(path, "r", encoding="utf-8")

def load_ids_from_filtered(path: str, id_col: str) -> List[str]:
	with open(path, "r", encoding="utf-8", newline="") as fh:
		r = csv.DictReader(fh, delimiter="\t")
		if id_col not in (r.fieldnames or []):
			sys.exit(f"[ERR] ID column '{id_col}' not found in {path}")
		ids: List[str] = []
		for row in r:
			val = (row.get(id_col, "") or "").strip()
			if not val or val == NA:
				continue
			ids.append(val)
		return ids

def fasta_iter(handle) -> Iterator[Tuple[str, str, str]]:
	"""
	yield (primary_id, full_header_without_gt, sequence)
	primary_id is the first token of the header, split by whitespace
	"""
	header = None
	seq_chunks: List[str] = []
	for line in handle:
		if not line:
			continue
		if line.startswith(">"):
			if header is not None:
				primary = header.split()[0]
				yield primary, header, "".join(seq_chunks)
			header = line[1:].rstrip("\n").strip()
			seq_chunks = []
		else:
			seq_chunks.append(line.strip())
	if header is not None:
		primary = header.split()[0]
		yield primary, header, "".join(seq_chunks)

def ensure_outdir(outdir: str) -> str:
	path = os.path.abspath(outdir)
	os.makedirs(path, exist_ok=True)
	return path

# helpers to classify full vs partial from prodigal-like metadata in headers
def _parse_meta(meta: str) -> Dict[str, str]:
	m: Dict[str, str] = {}
	for kv in meta.split(";"):
		kv = kv.strip()
		if "=" in kv:
			k, v = kv.split("=", 1)
			m[k] = v
	return m

def _is_partial_from_header(full_header_without_gt: str) -> bool:
	# Example header segment:
	# ... # 1 # 351 # 1 # ID=...;partial=10;start_type=Edge;...
	parts = [p.strip() for p in full_header_without_gt.split("#")]
	meta = parts[4] if len(parts) > 4 else ""
	m = _parse_meta(meta)
	partial = m.get("partial", "00")
	start_type = m.get("start_type", "")
	return (partial in {"10","01","11"}) or (start_type == "Edge")

# NEW split the original annotation table into full vs partial using captured ID sets
def _write_split_annotations(filtered_tsv: str, id_col: str, outdir: str, base: str, full_ids: Set[str], partial_ids: Set[str]) -> None:
	full_out = os.path.join(outdir, base + "_full_sequences.txt")
	part_out = os.path.join(outdir, base + "_partial_sequences.txt")

	with open(filtered_tsv, "r", encoding="utf-8") as fh:
		lines = fh.read().splitlines()

	if not lines:
		print(f"[INFO] {os.path.basename(filtered_tsv)} is empty, wrote headers only for split annotation files")
		with open(full_out, "w", encoding="utf-8", newline="") as wf, \
		     open(part_out, "w", encoding="utf-8", newline="") as wp:
			pass
		return

	header = lines[0]
	cols = header.split("\t")
	try:
		idx = cols.index(id_col)
	except ValueError:
		sys.exit(f"[ERR] ID column '{id_col}' not found in {filtered_tsv}")

	full_lines: List[str] = []
	part_lines: List[str] = []

	for ln in lines[1:]:
		if not ln:
			continue
		fields = ln.split("\t")
		if idx >= len(fields):
			continue
		_id = fields[idx].strip()
		if _id in full_ids:
			full_lines.append(ln)
		elif _id in partial_ids:
			part_lines.append(ln)

	# always write files, even if empty, write header only
	with open(full_out, "w", encoding="utf-8", newline="") as wf:
		wf.write(header + "\n")
		if full_lines:
			wf.write("\n".join(full_lines) + "\n")
	print(f"[OK] wrote {full_out} with {len(full_lines)} data rows")

	with open(part_out, "w", encoding="utf-8", newline="") as wp:
		wp.write(header + "\n")
		if part_lines:
			wp.write("\n".join(part_lines) + "\n")
	print(f"[OK] wrote {part_out} with {len(part_lines)} data rows")

def main():
	args = parse_args()
	filtered = os.path.abspath(args.filtered)
	fasta = os.path.abspath(args.fasta)
	outdir = os.path.dirname(filtered) if args.outdir is None else ensure_outdir(args.outdir)

	assert os.path.isfile(filtered), f"[ERR] filtered file not found, {filtered}"
	assert os.path.isfile(fasta), f"[ERR] fasta file not found, {fasta}"

	target_ids_list = load_ids_from_filtered(filtered, args.id_col)
	if not target_ids_list:
		print(f"[INFO] {os.path.basename(filtered)} has no valid IDs in column '{args.id_col}', nothing to do")
		return

	target_ids: Set[str] = set(target_ids_list)
	base = os.path.basename(filtered)
	out_path = os.path.join(outdir, base + ".fasta")
	headers_path = os.path.join(outdir, base + ".matched_headers.txt")
	out_full_path = os.path.join(outdir, base + "_full_sequences.fasta")
	out_part_path = os.path.join(outdir, base + "_partial_sequences.fasta")

	found = 0
	total = len(target_ids)
	wrote_any = False
	missing: Set[str] = set(target_ids)
	matched_headers: List[str] = []
	full_ids: Set[str] = set()     # NEW
	partial_ids: Set[str] = set()  # NEW

	with open(out_path, "w", encoding="utf-8", newline="") as out_fh, \
	     open(out_full_path, "w", encoding="utf-8", newline="") as wf_full, \
	     open(out_part_path, "w", encoding="utf-8", newline="") as wf_part:
		with open_text(fasta) as fh:
			for pid, header, seq in fasta_iter(fh):
				if pid in target_ids:
					if not wrote_any:
						print(f"[INFO] extracting {total} targets from {os.path.basename(fasta)}")
						wrote_any = True
					# write the main matched FASTA
					out_fh.write(f">{header}\n")
					for i in range(0, len(seq), 60):
						out_fh.write(seq[i:i+60] + "\n")
					# write matched header to sidecar
					matched_headers.append(f">{header}")
					# classify and write to full or partial FASTA, track IDs
					if _is_partial_from_header(header):
						wf_part.write(f">{header}\n")
						for i in range(0, len(seq), 60):
							wf_part.write(seq[i:i+60] + "\n")
						partial_ids.add(pid)  # NEW
					else:
						wf_full.write(f">{header}\n")
						for i in range(0, len(seq), 60):
							wf_full.write(seq[i:i+60] + "\n")
						full_ids.add(pid)     # NEW
					found += 1
					if pid in missing:
						missing.remove(pid)

	if found == 0:
		try:
			os.remove(out_path)
			os.remove(out_full_path)
			os.remove(out_part_path)
		except OSError:
			pass
		print(f"[INFO] 0 of {total} IDs found, no FASTA written")
	else:
		print(f"[OK] wrote {out_path}, found {found} of {total} IDs")
		print(f"[OK] wrote {out_full_path} and {out_part_path}")

	# write sidecar with exact matched headers
	if matched_headers:
		with open(headers_path, "w", encoding="utf-8") as hh:
			hh.write("\n".join(matched_headers) + "\n")
		print(f"[OK] wrote {headers_path} with {len(matched_headers)} matched header lines")

	# NEW write paired annotation tables for full and partial, preserving original order and formatting
	if full_ids or partial_ids:
		_write_split_annotations(filtered, args.id_col, outdir, base, full_ids, partial_ids)

	if args.write_missing and missing:
		miss_path = os.path.join(outdir, base + ".missing_ids.txt")
		with open(miss_path, "w", encoding="utf-8") as mfh:
			mfh.write("\n".join(sorted(missing)) + "\n")
		print(f"[INFO] missing IDs list written to {miss_path}")

if __name__ == "__main__":
	main()