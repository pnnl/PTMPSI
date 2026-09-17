#!/usr/bin/env python3
# tabs are used for indentation in this file

"""
Filter a merged kofam + pfam table by KO and PF.
The merged table may already include cat_gtdb:*, cat_nr:*, and virus:* columns. These are preserved.

Inputs
	--merged, merged table produced by clean_join.py, tab delimited
	--ko, a Kxxxxx code or NA to disable KO filtering
	--pf, a PFxxxxx code or NA to disable PF filtering
	--outdir, optional output directory

Output
	<sample>-filtered_KO-<KO or NA>_PF-<PF or NA>.txt
	Only written when at least one row matches.
"""

import argparse
import csv
import os
import re
from typing import Dict, List, Tuple

NA = "NA"

def parse_args():
	parser = argparse.ArgumentParser(description="Filter merged kofam + pfam table by KO and PF.")
	parser.add_argument("--merged", required=True, help="Merged table from clean_join.py")
	parser.add_argument("--ko", required=True, help="Kxxxxx or NA")
	parser.add_argument("--pf", required=True, help="PFxxxxx or NA")
	parser.add_argument("--outdir", default=".", help="Output directory")
	return parser.parse_args()

def ensure_outdir(d: str) -> str:
	path = os.path.abspath(d)
	os.makedirs(path, exist_ok=True)
	return path

def guess_sample_from_merged_filename(path: str) -> str:
	base = os.path.basename(path)
	m = re.match(r"(.+?)-merged_kofam_pfam_", base)
	if m:
		return m.group(1)
	name = os.path.splitext(base)[0]
	return name.split("-merged_kofam_pfam_")[0]

def pf_strip_version(pf: str) -> str:
	if not pf or pf == NA:
		return NA
	pf = pf.strip().upper()
	m = re.match(r"(PF\d+)", pf)
	return m.group(1) if m else pf

def ko_normalize(ko: str) -> str:
	if not ko or ko == NA:
		return NA
	ko = ko.strip().upper()
	m = re.match(r"(K\d{5})", ko)
	return m.group(1) if m else ko

def read_merged(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
	with open(path, "r", encoding="utf-8", newline="") as fh:
		r = csv.DictReader(fh, delimiter="\t")
		header = r.fieldnames or []
		rows: List[Dict[str, str]] = []
		for row in r:
			for k, v in list(row.items()):
				if v is None or v == "":
					row[k] = NA
			rows.append(row)
	return header, rows

def find_ko_col(rows: List[Dict[str, str]]) -> str:
	cand = ["kofam:KO", "KO"]
	if not rows:
		return cand[0]
	keys = list(rows[0].keys())
	for c in cand:
		if c in keys:
			return c
	for c in keys:
		cl = c.lower()
		if cl.endswith(":ko") or cl == "ko":
			return c
	raise SystemExit("KO column not found. Expected 'kofam:KO' or a column ending with ':KO'.")

def find_pf_col(rows: List[Dict[str, str]]) -> str:
	cand = ["pfamsearch:target accession"]
	if not rows:
		return cand[0]
	keys = list(rows[0].keys())
	for c in cand:
		if c in keys:
			return c
	for c in keys:
		cl = c.lower()
		if "pfam" in cl and "accession" in cl:
			return c
	raise SystemExit("PF accession column not found. Expected 'pfamsearch:target accession' or a column containing 'pfam' and 'accession'.")

def filter_rows(rows: List[Dict[str, str]], want_ko: str, want_pf: str, ko_col: str, pf_col: str) -> List[Dict[str, str]]:
	out: List[Dict[str, str]] = []
	for row in rows:
		ok_ko = True
		ok_pf = True
		if want_ko != NA:
			ok_ko = ko_normalize(row.get(ko_col, NA)) == want_ko
		if want_pf != NA:
			pfval = pf_strip_version(row.get(pf_col, NA))
			ok_pf = pfval == want_pf
		if ok_ko and ok_pf:
			out.append(row)
	return out

def write_tsv(path: str, header: List[str], rows: List[Dict[str, str]]) -> None:
	with open(path, "w", encoding="utf-8", newline="") as out:
		w = csv.DictWriter(out, fieldnames=header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
		w.writeheader()
		for r in rows:
			for h in header:
				r.setdefault(h, NA)
			w.writerow(r)

def main():
	args = parse_args()
	outdir = ensure_outdir(args.outdir)

	want_ko = ko_normalize(args.ko if args.ko.strip().lower() != "na" else NA)
	want_pf = pf_strip_version(args.pf if args.pf.strip().lower() != "na" else NA)

	header, rows = read_merged(args.merged)
	if not rows:
		raise SystemExit("No rows found in merged file")

	ko_col = find_ko_col(rows) if want_ko != NA else ""
	pf_col = find_pf_col(rows) if want_pf != NA else ""
	sample = guess_sample_from_merged_filename(args.merged)

	filt_rows = filter_rows(rows, want_ko, want_pf, ko_col, pf_col)
	total = len(rows)
	matched = len(filt_rows)
	print(f"[INFO] sample={sample}, total_rows={total}, matched_rows={matched}, KO={want_ko}, PF={want_pf}")

	if matched == 0:
		print(f"[INFO] sample={sample}, no rows matched, skipping outputs")
		return

	ko_tag = want_ko if want_ko != NA else "NA"
	pf_tag = want_pf if want_pf != NA else "NA"
	out_path = os.path.join(outdir, f"{sample}-filtered_KO-{ko_tag}_PF-{pf_tag}.txt")
	write_tsv(out_path, header, filt_rows)
	print(f"[OK] Wrote {out_path}")

if __name__ == "__main__":
	main()