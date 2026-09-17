#!/usr/bin/env python3

import argparse
import shutil
import subprocess
from Bio.PDB import PDBParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1

def parse_args():
	p = argparse.ArgumentParser(
		description="Phase 4.1 + 4.2: map CP12 PDB residues to PF02672 HMM columns"
	)
	p.add_argument("--pdb", required=True, help="PDB file with CP12 chains")
	p.add_argument("--chains", required=True, help="Comma-separated CP12 chain IDs, e.g. C,H,M,P")
	p.add_argument("--hmm", required=True, help="PF02672 HMM file")
	p.add_argument("--system", required=True, help="System name, e.g. MED4_0003")
	p.add_argument("--out-prefix", required=True, help="Output prefix")
	p.add_argument(
		"--hmmalign-bin",
		default="hmmalign",
		help="Path to HMMER hmmalign executable. Default: hmmalign on PATH."
	)
	return p.parse_args()

def extract_chain_sequences(pdb_file, chains, system):
	parser = PDBParser(QUIET=True)
	structure = parser.get_structure("struct", pdb_file)

	records = []
	mappings = []

	for chain_id in chains:
		seq = []
		res_map = []

		for model in structure:
			for chain in model:
				if chain.id != chain_id:
					continue
				for res in chain:
					if res.id[0] != " ":
						continue
					if "CA" not in res:
						continue

					aa3 = res.get_resname()
					try:
						aa1 = seq1(aa3)
					except Exception:
						continue

					seq.append(aa1)
					res_map.append(res.id[1])

		if not seq:
			continue

		record = SeqRecord(
			Seq("".join(seq)),
			id=f"{system}|CP12|chain_{chain_id}",
			description=""
		)
		records.append(record)

		for i, res_id in enumerate(res_map, start=1):
			mappings.append({
				"System": system,
				"CP12_Chain": chain_id,
				"ResID": res_id,
				"SeqPos": i,
				"AA": seq[i - 1]
			})

	return records, mappings

def run_hmmalign(hmm, fasta, out_sto, hmmalign_bin):
	if shutil.which(hmmalign_bin) is None:
		raise SystemExit(
			"Cannot find hmmalign. Install HMMER or pass --hmmalign-bin /path/to/hmmalign. "
			"For conda/mamba: mamba install -c bioconda hmmer. For Homebrew: brew install hmmer."
		)
	cmd = [hmmalign_bin, hmm, fasta]
	with open(out_sto, "w") as fh:
		subprocess.run(cmd, stdout=fh, check=True)

def parse_stockholm(sto_file):
	hmm_cols = {}
	with open(sto_file) as fh:
		for line in fh:
			if line.startswith("#") or line.strip() == "" or line.startswith("//"):
				continue
			fields = line.strip().split()
			if len(fields) != 2:
				continue

			seq_id, aln = fields
			col = 0
			seq_pos = 0
			pos_map = {}

			for c in aln:
				if c.isupper() or c == "-":
					col += 1
				if c.isupper():
					seq_pos += 1
					pos_map[seq_pos] = col

			hmm_cols[seq_id] = pos_map
	return hmm_cols

def main():
	args = parse_args()
	chains = args.chains.split(",")

	records, res_map = extract_chain_sequences(args.pdb, chains, args.system)

	fasta_file = f"{args.out_prefix}_cp12_chains.faa"
	SeqIO.write(records, fasta_file, "fasta")

	sto_file = f"{args.out_prefix}_cp12_hmmalign.sto"
	run_hmmalign(args.hmm, fasta_file, sto_file, args.hmmalign_bin)

	hmm_map = parse_stockholm(sto_file)

	out_tsv = f"{args.out_prefix}_cp12_resid_to_hmmcol.tsv"
	with open(out_tsv, "w") as out:
		out.write("System\tCP12_Chain\tResID\tSeqPos\tAA\tHMM_col\n")
		for row in res_map:
			seq_id = f"{row['System']}|CP12|chain_{row['CP12_Chain']}"
			hmm_col = hmm_map.get(seq_id, {}).get(row["SeqPos"], "NA")
			out.write(
				f"{row['System']}\t{row['CP12_Chain']}\t{row['ResID']}\t"
				f"{row['SeqPos']}\t{row['AA']}\t{hmm_col}\n"
			)

	print(f"Wrote residue-to-HMM mapping: {out_tsv}")

if __name__ == "__main__":
	main()