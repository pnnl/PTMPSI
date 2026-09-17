import concurrent.futures
from itertools import combinations
import ptmpsi
from ptmpsi.protein import Protein
import os

ptm_dict = {"SNO": "nitrosylation", "SOH": "sulfenylation", "SSG": "glutathionylation", "None": "reduction"}


gap2 = [("None", "None", "None", "None")]
for ptm in ["SNO", "SOH", "SSG"]:
    for ituple in range(1,5):
        combo = list(combinations([0, 1, 2, 3], r=ituple))
        for icombo in combo:
            temp = ["None", "None", "None", "None"]
            for ipos in icombo:
                temp[ipos] = ptm
            gap2.append(tuple(temp))

prk = [("None", "None"), ("SS", "SS")]
for ptm in ["SNO", "SOH", "SSG"]:
    for ituple in range(1,3):
        combo = list(combinations([0, 1], r=ituple))
        for icombo in combo:
            temp = ["None", "None"]
            for ipos in icombo:
                temp[ipos] = ptm
            prk.append(tuple(temp))

cp12 = [ ("None", "None"), ("SS", "SS") ]

gap2_chains = ["A", "D", "E", "F", "I", "K", "N", "O"]
prk_chains = ["B", "L", "G", "J"]
cp12_chains = ["C", "H", "M", "P"]
gap2_cys = ["CYS20", "CYS78", "CYS157", "CYS161"]
prk_cys = ["CYX201", "CYX225"]
cp12_cys = ["CYX56", "CYX65"]

full_combo_list = []
for igap in gap2:
    for iprk in prk:
        for icp12 in cp12:
            case = []
            for chain in gap2_chains:
                for i, ptm in enumerate(igap):
                    if "None" in ptm: continue
                    case.append(tuple((f"{chain}:{gap2_cys[i]}", f"{ptm_dict[ptm]}")))
            for chain in prk_chains:
                for i, ptm in enumerate(iprk):
                    if ptm == "SS": continue
                    case.append(tuple((f"{chain}:{prk_cys[i]}", f"{ptm_dict[ptm]}")))
            for chain in cp12_chains:
                for i, ptm in enumerate(icp12):
                    if "SS" in ptm: continue
                    case.append(tuple((f"{chain}:{cp12_cys[i]}", f"{ptm_dict[ptm]}")))
            full_combo_list.append(tuple(case.copy()))

machine = 'frontier'
account = 'bip258'
long_partition='extended'
submit_command = 'sbatch'

protein = Protein(filename="hm2_dark_complex_modeller_CYX_250101.pdb")

# Parallelized
gromacs_directory = '/path/to/ptmpsi-production-runs'
gromacs_path = os.path.join(gromacs_directory, "gromacs_sims/HM2_NAD")
gromacs_options = ptmpsi.slurm.machines[machine].partitions[long_partition].options["gromacs"]
gromacs_options['ntasks'] = 8

def process_full_combo_list():
    return protein.gen_ptm_files(
        [full_combo_list],
        bundling=True,
        auto_submit_ti=False,
        lenmd=100,
        lennvt=500,
        lennpt=500,
        path=gromacs_path,
        machine=machine,
        account=account,
        **gromacs_options,
        partition=long_partition,
        time=24,
	cofactor="NAD_HM2.pdb"
    )

# Use ProcessPoolExecutor for parallel execution
with concurrent.futures.ProcessPoolExecutor() as executor:
    future = executor.submit(process_full_combo_list)
    uid = future.result()
