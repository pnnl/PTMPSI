#!/usr/bin/env python3
"""
Super merge v10
1) Reformat KOFAM and PFAM to top hit per gene and standardize IDs
2) Merge by new_gene_name with outer or inner join
3) Append CAT GTDB and CAT NR contig2classification by mapping new_gene_name -> contig id
4) Optionally append geNomad viral summary by contig id
5) Write one merged TSV with a descriptive filename
6) Append full input file paths as provenance columns at the end of the table

Rules
- new_gene_name
	- If original gene contains 'PNNL', keep as is
	- Else prefix sample_id + '_PNNL_' + gene
	- sample_id defaults to the token before first underscore in the KOFAM file name, or you can pass --sample to override
- sample_name column is added as the first column
- Headers from sources are prefixed with source tags
	- kofam:*
	- pfamsearch:*
	- cat_gtdb:*, cat_nr:*
	- virus:*
- Missing fields are filled with 'NA'
- Alignment fix
	- When a data line begins with a leading tab and a real header is present, shift fields so columns line up
"""

import argparse
import csv
import os
import re
from typing import Dict, List, Tuple, Optional

NA = "NA"

def ensure_outdir(d: str) -> str:
	path = os.path.abspath(d)
	os.makedirs(path, exist_ok=True)
	return path

def basename_noext(p: str) -> str:
	return os.path.splitext(os.path.basename(p))[0]

def infer_sample_from_fname(path: str) -> str:
	"""
	Default sample id is the first underscore delimited token of the input filename.
	Example: renamed_prodigal_out_20240911-PNNL-MPF-Seq129_megahit... -> 20240911-PNNL-MPF-Seq129
	"""
	base = os.path.basename(path)
	tok = base.split("_")[0]
	return tok

def new_gene_name_for(sample_id: str, gene: str) -> str:
	"""
	If gene contains 'PNNL', return gene, else prefix sample_id + '_PNNL_' + gene
	"""
	if "PNNL" in gene:
		return gene
	return f"{sample_id}_PNNL_{gene}"

def contig_from_new_gene(new_gene: str) -> str:
	"""
	Extract contig id like k141_123456 from new_gene_name using '(k[0-9]+_[0-9]+)'
	"""
	if not new_gene:
		return ""
	m = re.search(r"(k[0-9]+_[0-9]+)", new_gene)
	return m.group(1) if m else ""

def write_tsv(path: str, header: List[str], rows: List[Dict[str, str]]) -> None:
	with open(path, "w", encoding="utf-8", newline="") as out:
		w = csv.DictWriter(out, fieldnames=header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
		w.writeheader()
		for r in rows:
			for h in header:
				r.setdefault(h, NA)
			w.writerow(r)

def _align_to_header(parts: List[str], header_cols: List[str]) -> List[str]:
	"""
	Fix rows that have one extra leading cell caused by a leading tab
	or a bare KOfam significance marker. Then pad or trim to header length.
	"""
	if not header_cols:
		return parts
	h = len(header_cols)
	if len(parts) == h + 1 and parts:
		first = parts[0]
		# drop empty first cell or a bare marker run such as *, !!, --, +, ?
		if first == "" or re.fullmatch(r"[*!+\-?]+", first or ""):
			parts = parts[1:]
	if len(parts) < h:
		parts = parts + [NA] * (h - len(parts))
	elif len(parts) > h:
		parts = parts[:h]
	return parts


def _clean_field_name(name: str) -> str:
	"""
	Remove quotes and trim spaces from a field name.
	"""
	return name.replace('"', "").strip()

def read_kofam_top1(path: str, forced_sample: Optional[str]) -> Tuple[List[str], List[Dict[str, str]], str]:
	"""
	Read KOFAM TSV, keep the first hit per 'gene name'.
	Return header, rows with source prefixed columns and standardized new_gene_name, and sample_id used.
	"""
	sample_id = forced_sample or infer_sample_from_fname(path)
	keep: Dict[str, Dict[str, str]] = {}
	with open(path, "r", encoding="utf-8", errors="ignore") as fh:
		header_cols: List[str] = []
		for raw in fh:
			if not raw.strip():
				continue
			if raw.startswith("#"):
				h = raw.lstrip("#").strip()
				if "\t" in h and "gene name" in h:
					header_cols = [_clean_field_name(c) for c in h.split("\t")]
				continue
			# normalize leading significance markers, single or repeated, optional space or tab
			raw = re.sub(r'^[*!+\-?]+[ \t]?', '', raw)
			parts = [p.rstrip() for p in raw.split('\t')]
			if not parts or len(parts) < 1:
				continue
			if header_cols and "gene name" in header_cols:
				parts = _align_to_header(parts, header_cols)
				row_map: Dict[str, str] = {}
				for i, c in enumerate(header_cols):
					row_map[c] = parts[i] if i < len(parts) else NA
				gene = (row_map.get("gene name", "") or "").strip()
				if not gene:
					continue
				if gene in keep:
					continue
				new_id = new_gene_name_for(sample_id, gene)
				out = {"new_gene_name": new_id, "kofam:gene_name": gene}
				for k, v in row_map.items():
					if k == "gene name":
						continue
					out[f'kofam:{_clean_field_name(k)}'] = v if v != "" else NA
				out["sample_name"] = sample_id
				keep[gene] = out
			else:
				# fallback positions: 0 gene, 1 KO, 2 thrshld, 3 score, 4 E-value, 5 KO definition
				gene = re.sub(r'^[*!+\-?]+', '', parts[0]).strip()
				if not gene:
					continue
				if gene in keep:
					continue
				out = {
					"new_gene_name": new_gene_name_for(sample_id, gene),
					"kofam:gene_name": gene,
					"kofam:KO": parts[1] if len(parts) > 1 else NA,
					"kofam:thrshld": parts[2] if len(parts) > 2 else NA,
					"kofam:score": parts[3] if len(parts) > 3 else NA,
					"kofam:E-value": parts[4] if len(parts) > 4 else NA,
					"kofam:KO definition": parts[5] if len(parts) > 5 else NA,
					"sample_name": sample_id
				}
				keep[gene] = out
	# header
	if keep:
		extra_cols = sorted([c for c in next(iter(keep.values())).keys() if c not in {"sample_name","new_gene_name"}])
		header = ["sample_name", "new_gene_name"] + extra_cols
	else:
		header = ["sample_name", "new_gene_name"]
	rows = list(keep.values())
	return header, rows, sample_id


def read_pfam_top1(path: str, forced_sample: Optional[str], sample_id_hint: Optional[str]) -> Tuple[List[str], List[Dict[str, str]], str]:
	"""
	Read PFAM table, keep the first hit per query name.
	Accepts tab-delimited or space-aligned HMMER tables.
	Returns header, rows with source-prefixed columns and standardized new_gene_name, and sample_id used.
	"""
	sample_id = forced_sample or sample_id_hint or infer_sample_from_fname(path)
	keep: Dict[str, Dict[str, str]] = {}
	with open(path, "r", encoding="utf-8", errors="ignore") as fh:
		header_cols: List[str] = []
		header_is_tab = False
		for raw in fh:
			line = raw.rstrip("\n")
			if not line.strip():
				continue
			if line.startswith("#"):
				# Detect a true tab header only, otherwise use positional fallback
				h = line.lstrip("#").strip()
				if "\t" in h and "query name" in h:
					header_cols = [_clean_field_name(c) for c in h.split("\t")]
					header_is_tab = True
				continue

			# choose splitting strategy for data rows
			if header_is_tab:
				parts = [p.rstrip() for p in line.split("\t")]
				if header_cols:
					parts = _align_to_header(parts, header_cols)
			else:
				# space-aligned, collapse runs of whitespace
				parts = line.strip().split()

			if not parts:
				continue

			if header_is_tab and header_cols and "query name" in header_cols:
				row_map: Dict[str, str] = {}
				for i, c in enumerate(header_cols):
					row_map[c] = parts[i] if i < len(parts) else NA
				gene = (row_map.get("query name", "") or "").strip()
				if not gene or gene in keep:
					continue
				new_id = new_gene_name_for(sample_id, gene)
				out = {"new_gene_name": new_id, "pfamsearch:query_name": gene}
				for k, v in row_map.items():
					if k == "query name":
						continue
					out[f'pfamsearch:{_clean_field_name(k)}'] = v if v != "" else NA
				out["sample_name"] = sample_id
				keep[gene] = out
			else:
				# positional fallback for space-aligned PFAM/HMMER tables
				# expected columns: target name, accession, query name, accession, E-value, score, ...
				if len(parts) < 6:
					continue
				gene = parts[2].strip()
				if not gene or gene in keep:
					continue
				out = {
					"new_gene_name": new_gene_name_for(sample_id, gene),
					"pfamsearch:query_name": gene,
					"pfamsearch:target_name": parts[0],
					"pfamsearch:target_accession": parts[1],
					"pfamsearch:E-value": parts[4],
					"pfamsearch:score": parts[5],
					"sample_name": sample_id
				}
				keep[gene] = out

	# header
	if keep:
		extra_cols = sorted([c for c in next(iter(keep.values())).keys() if c not in {"sample_name", "new_gene_name"}])
		header = ["sample_name", "new_gene_name"] + extra_cols
	else:
		header = ["sample_name", "new_gene_name"]
	rows = list(keep.values())
	return header, rows, sample_id

def merge_kofam_pfam(kh: List[str], kr: List[Dict[str,str]], ph: List[str], pr: List[Dict[str,str]], how: str) -> Tuple[List[str], List[Dict[str,str]]]:
	"""
	Join on new_gene_name, how in {'outer','inner'}
	Fill NA for missing
	"""
	kmap = {r["new_gene_name"]: r for r in kr}
	pmap = {r["new_gene_name"]: r for r in pr}
	if how == "inner":
		keys = sorted(set(kmap).intersection(pmap))
	else:
		keys = sorted(set(kmap).union(pmap))
	cols = sorted({c for c in (kh + ph) if c not in {"sample_name","new_gene_name"}})
	out_header = ["sample_name","new_gene_name"] + cols
	out_rows: List[Dict[str,str]] = []
	for key in keys:
		row = {"sample_name": kmap.get(key, pmap.get(key, {})).get("sample_name", NA), "new_gene_name": key}
		if key in kmap:
			for c, v in kmap[key].items():
				if c in {"sample_name","new_gene_name"}:
					continue
				row[c] = v
		else:
			for c in cols:
				if c.startswith("kofam:"):
					row.setdefault(c, NA)
		if key in pmap:
			for c, v in pmap[key].items():
				if c in {"sample_name","new_gene_name"}:
					continue
				row[c] = v
		else:
			for c in cols:
				if c.startswith("pfamsearch:"):
					row.setdefault(c, NA)
		out_rows.append(row)
	return out_header, out_rows


def read_cat_c2c(cat_path: str, label: str) -> Dict[str, Dict[str, str]]:
	"""
	Read contig2classification into map contig -> dict with prefixed keys
	"""
	cmap: Dict[str, Dict[str, str]] = {}
	with open(cat_path, "r", encoding="utf-8", errors="ignore") as fh:
		for raw in fh:
			if not raw.strip():
				continue
			if raw.startswith("#"):
				continue
			parts = raw.rstrip("\n").split("\t")
			if not parts:
				continue
			contig = parts[0].strip()
			if not contig:
				continue
			cls = parts[1].strip() if len(parts) > 1 and parts[1] else NA
			reason = parts[2].strip() if len(parts) > 2 and parts[2] else NA
			lineage = parts[3].strip() if len(parts) > 3 and parts[3] else NA
			linscores = parts[4].strip() if len(parts) > 4 and parts[4] else NA
			cmap[contig] = {
				f"cat_{label}:classification": cls,
				f"cat_{label}:reason": reason,
				f"cat_{label}:lineage": lineage,
				f"cat_{label}:lineage_scores": linscores
			}
	return cmap

def append_cat_by_contig(header: List[str], rows: List[Dict[str,str]], cat_path: Optional[str], label: str) -> Tuple[List[str], List[Dict[str,str]]]:
	if not cat_path:
		return header, rows
	cmap = read_cat_c2c(cat_path, label)
	cat_cols = [
		f"cat_{label}:classification",
		f"cat_{label}:reason",
		f"cat_{label}:lineage",
		f"cat_{label}:lineage_scores"
	]
	out_header = list(header)
	for c in cat_cols:
		if c not in out_header:
			out_header.append(c)
	out_rows: List[Dict[str,str]] = []
	for r in rows:
		contig = contig_from_new_gene(r.get("new_gene_name",""))
		extra = cmap.get(contig, {c: NA for c in cat_cols})
		nr = dict(r)
		nr.update(extra)
		out_rows.append(nr)
	return out_header, out_rows

VIRUS_KEEP = ["seq_name","length","topology","coordinates","n_genes","genetic_code","virus_score","fdr","n_hallmarks","marker_enrichment","taxonomy"]

def read_viral_summary(vpath: str) -> Dict[str, Dict[str, str]]:
	vmap: Dict[str, Dict[str, str]] = {}
	if not vpath:
		return vmap
	with open(vpath, "r", encoding="utf-8", errors="ignore") as fh:
		r = csv.DictReader(fh, delimiter="\t")
		if not r.fieldnames:
			return vmap
		for row in r:
			seq = (row.get("seq_name","") or "").strip()
			if not seq:
				continue
			out: Dict[str, str] = {}
			for c in VIRUS_KEEP:
				val = row.get(c, NA)
				out[f"virus:{c}"] = val if val not in [None,""] else NA
			vmap[seq] = out
	return vmap

def append_virus_by_contig(header: List[str], rows: List[Dict[str,str]], vpath: Optional[str]) -> Tuple[List[str], List[Dict[str,str]]]:
	if not vpath:
		return header, rows
	vmap = read_viral_summary(vpath)
	vcols = [f"virus:{c}" for c in VIRUS_KEEP]
	out_header = list(header)
	for c in vcols:
		if c not in out_header:
			out_header.append(c)
	out_rows: List[Dict[str,str]] = []
	for r in rows:
		contig = contig_from_new_gene(r.get("new_gene_name",""))
		extra = vmap.get(contig, {c: NA for c in vcols})
		nr = dict(r)
		nr.update(extra)
		out_rows.append(nr)
	return out_header, out_rows

# -------------- main --------------

def parse_args():
	ap = argparse.ArgumentParser(description="Merge KOFAM, PFAM, optional CAT GTDB, CAT NR, virus summary")
	ap.add_argument("--kofam", required=True, help="KOFAM output TSV")
	ap.add_argument("--pfam", required=True, help="PFAM output TSV")
	ap.add_argument("--join", choices=["outer","inner"], default="outer", help="Join type for KOFAM and PFAM")
	ap.add_argument("--sample", default=None, help="Override sample_id for naming and new_gene_name prefixing")
	ap.add_argument("--outdir", default=".", help="Output directory")
	ap.add_argument("--cat_gtdb", default=None, help="CAT GTDB contig2classification file")
	ap.add_argument("--cat_nr", default=None, help="CAT NR contig2classification file")
	ap.add_argument("--virus", default=None, help="geNomad viral summary TSV")
	return ap.parse_args()

def main():
	args = parse_args()
	outdir = ensure_outdir(args.outdir)

	# absolute input paths for provenance
	kofam_path = os.path.abspath(args.kofam)
	pfam_path = os.path.abspath(args.pfam)
	cat_gtdb_path = os.path.abspath(args.cat_gtdb) if args.cat_gtdb else NA
	cat_nr_path = os.path.abspath(args.cat_nr) if args.cat_nr else NA
	virus_path = os.path.abspath(args.virus) if args.virus else NA

	# read kofam and pfam
	kh, kr, k_sample = read_kofam_top1(args.kofam, args.sample)
	ph, pr, p_sample = read_pfam_top1(args.pfam, args.sample, k_sample)

	# sample consistency and choice for filenames
	sample = args.sample or k_sample or p_sample
	if k_sample != p_sample:
		print(f"[WARN] sample names differ, kofam={k_sample}, pfam={p_sample}. Using {sample} for outputs.")

	# merge
	mh, mr = merge_kofam_pfam(kh, kr, ph, pr, args.join)

	# append CAT by contig, contig2classification only
	tags: List[str] = []
	if args.cat_gtdb:
		mh, mr = append_cat_by_contig(mh, mr, args.cat_gtdb, "gtdb")
		tags.append("gtdb")
	if args.cat_nr:
		mh, mr = append_cat_by_contig(mh, mr, args.cat_nr, "nr")
		tags.append("nr")

	# append viral
	if args.virus:
		mh, mr = append_virus_by_contig(mh, mr, args.virus)
		tags.append("viral")

	# provenance columns at the end of the table
	prov_cols = [
		"input_kofam_path",
		"input_pfam_path",
		"input_cat_gtdb_path",
		"input_cat_nr_path",
		"input_virus_path"
	]

	for c in prov_cols:
		if c not in mh:
			mh.append(c)

	for r in mr:
		r["input_kofam_path"] = kofam_path
		r["input_pfam_path"] = pfam_path
		r["input_cat_gtdb_path"] = cat_gtdb_path
		r["input_cat_nr_path"] = cat_nr_path
		r["input_virus_path"] = virus_path

	# filename tag
	join_tag = "union" if args.join == "outer" else "intersection"
	extra_tag = "_".join(tags) if tags else ""
	name_parts = [f"{sample}-merged_kofam_pfam_{join_tag}"]
	if extra_tag:
		name_parts.append(extra_tag)
	out_name = "_".join(name_parts) + ".txt"
	out_path = os.path.join(outdir, out_name)

	# write
	write_tsv(out_path, mh, mr)
	print(f"[OK] Wrote {out_path}")
	print(f"[INFO] rows={len(mr)} cols={len(mh)}")

if __name__ == "__main__":
	main()