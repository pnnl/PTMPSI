#!/usr/bin/env python3
# tabs are used for indentation

import argparse
import csv
import os
import re
from typing import List, Dict, Tuple

NA = "NA"
CATEGORIES = ["Viruses","Bacteria","Archaea","Fungi","Viridiplantae","Metazoa","SAR"]

def _split_lineage(s: str) -> list:
	if not s or s == 'NA':
		return []
	return [t.strip() for t in s.split(';') if t.strip()]

def _longest_common_prefix(a: list, b: list) -> list:
	out = []
	for x, y in zip(a, b):
		if x == y:
			out.append(x)
		else:
			break
	return out

def _first_disagreement_index(a: list, b: list) -> int:
	n = min(len(a), len(b))
	for i in range(n):
		if a[i] != b[i]:
			return i
	return n

def _domain_to_category(domain_token: str) -> str:
	if not domain_token:
		return 'Unclassified'
	t = domain_token.strip()
	if t.startswith('d__Bacteria'):
		return 'Bacteria'
	if t.startswith('d__Archaea'):
		return 'Archaea'
	if t.startswith('d__Eukaryota'):
		return 'Unclassified'
	if 'Viruses' in t or t.startswith('d__Viruses'):
		return 'Unclassified'
	return 'Unclassified'

def resolve_tax_fields(record: dict) -> dict:
	virus_tax = record.get('virus:taxonomy', 'NA')
	gtdb_lin = record.get('cat_gtdb:lineage', 'NA')
	nr_lin = record.get('cat_nr:lineage', 'NA')

	out = {}
	if virus_tax and virus_tax != 'NA':
		out['Category'] = 'Viruses'
		out['SourceUsed'] = 'genomad'
		out['MatchBasis-greedy'] = virus_tax
		out['MatchBasis-conservative'] = virus_tax
		return out

	a = _split_lineage(gtdb_lin)
	b = _split_lineage(nr_lin)

	if a or b:
		if a and b:
			common = _longest_common_prefix(a, b)
			dis_idx = _first_disagreement_index(a, b)
			longer = a if len(a) >= len(b) else b
			m = min(len(a), len(b))
			mismatch = any(a[i] != b[i] for i in range(m))
			if mismatch:
				greedy = longer[:dis_idx]
			else:
				greedy = longer[:]
			cons = common
			cons_str = ';'.join(cons) if cons else 'NA'
			greedy_str = ';'.join(greedy) if greedy else 'NA'
			out['MatchBasis-conservative'] = cons_str
			out['MatchBasis-greedy'] = greedy_str
			domain_token = cons[1] if len(cons) > 1 else (greedy[1] if len(greedy) > 1 else '')
			out['Category'] = _domain_to_category(domain_token)
			out['SourceUsed'] = 'gtdb+nr'
			return out
		only = a if a else b
		lin_str = ';'.join(only) if only else 'NA'
		out['MatchBasis-conservative'] = lin_str if lin_str else 'NA'
		out['MatchBasis-greedy'] = lin_str if lin_str else 'NA'
		domain_token = only[1] if len(only) > 1 else ''
		out['Category'] = _domain_to_category(domain_token)
		out['SourceUsed'] = 'gtdb' if a else 'nr'
		return out

	out['Category'] = 'Unclassified'
	out['SourceUsed'] = 'NA'
	out['MatchBasis-greedy'] = 'NA'
	out['MatchBasis-conservative'] = 'NA'
	return out

def parse_args():
	p = argparse.ArgumentParser(description="Append taxonomy fields and optionally filter by targets")
	p.add_argument("--merged", required=True, help="Input table, function filtered or merged")
	p.add_argument("--sample", required=True, help="Sample ID for output naming")
	p.add_argument("--targets", default=",".join(CATEGORIES), help="Comma separated categories to keep")
	p.add_argument("--strict-taxid-only", action="store_true", help="Reserved, keep for compatibility")
	p.add_argument("--use_greedy_filter", choices=["yes","no"], default="no", help="If yes, match targets against MatchBasis-greedy lineage")
	p.add_argument("--skip_taxa_filter", choices=["yes","no"], default="no", help="If yes, do not drop rows by taxonomy, append columns only")
	p.add_argument("--outdir", default=".", help="Output directory")
	return p.parse_args()

def ensure_outdir(d: str) -> str:
	path = os.path.abspath(d)
	os.makedirs(path, exist_ok=True)
	return path

def read_table(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
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

def _row_matches_targets(row: dict, targets: list, use_greedy: bool) -> bool:
	cat = (row.get('Category') or '').strip()
	if cat == 'Viruses':
		return any(t.strip() == 'Viruses' for t in targets)
	if not use_greedy:
		return any(t.strip() == cat for t in targets)
	greedy = (row.get('MatchBasis-greedy') or '').strip()
	if not greedy and cat:
		return any(t.strip() == cat for t in targets)
	for t in targets:
		ts = t.strip()
		if not ts:
			continue
		if ts == cat:
			return True
		if ts in greedy:
			return True
		if ts.lower() in greedy.lower():
			return True
	return False

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
	targets = [t.strip() for t in args.targets.split(",") if t.strip()]
	targets = [t for t in targets if t in CATEGORIES]

	header, rows = read_table(args.merged)
	app_cols = ["Category","SourceUsed","MatchBasis-greedy","MatchBasis-conservative"]
	out_header = list(header)
	for c in app_cols:
		if c not in out_header:
			out_header.append(c)

	kept: List[Dict[str,str]] = []
	use_greedy = (args.use_greedy_filter == "yes")

	for r in rows:
		r2 = dict(r)
		nf = resolve_tax_fields(r2)
		r2.update(nf)

		if args.skip_taxa_filter == "yes":
			kept.append(r2)
			continue

		if _row_matches_targets(r2, targets, use_greedy):
			kept.append(r2)

	total = len(rows)
	matched = len(kept)
	mode = "strict" if args.strict_taxid_only else "lenient"
	skp = args.skip_taxa_filter
	grd = args.use_greedy_filter
	print(f"[INFO] sample={args.sample}, mode={mode}, skip_taxa={skp}, greedy={grd}, total_rows={total}, kept_rows={matched}, targets={'+'.join(targets)}")

	target_tag = "+".join(targets)
	out_path = os.path.join(outdir, f"{args.sample}-filtered_taxa-{target_tag}.txt")

	# Write policy:
	# skip_taxa=yes, always write, even if zero kept, so chain mode has a file
	# skip_taxa=no, write only if matched > 0
	if args.skip_taxa_filter == "yes":
		write_tsv(out_path, out_header, kept)
		print(f"[OK] Wrote {out_path}")
		return

	if matched == 0:
		print(f"[INFO] sample={args.sample}, no rows kept, skipping output")
		return

	write_tsv(out_path, out_header, kept)
	print(f"[OK] Wrote {out_path}")

if __name__ == "__main__":
	main()