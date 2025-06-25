EAC_frg = """$EAC 
   19    1    1    0
EAC   
    1 N    N         1    1    0    1    1   -0.542610    0.000000
    2 H    H         0    0    0    1    1    0.326920    0.000000
    3 CE   CT        0    0    0    1    1    0.049700    0.000000
    42HE   H1        0    0    0    1    1    0.081380    0.000000
    53HE   H1        0    0    0    1    1    0.081380    0.000000
    6 CD   CT        0    0    0    1    1   -0.158800    0.000000
    72HD   HC        0    0    0    1    1    0.084150    0.000000
    83HD   HC        0    0    0    1    1    0.084150    0.000000
    9 CG   CT        0    0    0    1    1   -0.116460    0.000000
   102HG   HC        0    0    0    1    1    0.080440    0.000000
   113HG   HC        0    0    0    1    1    0.080440    0.000000
   12 CB   CT        0    0    0    1    1   -0.096240    0.000000
   132HB   HC        0    0    0    1    1    0.077710    0.000000
   143HB   HC        0    0    0    1    1    0.077710    0.000000
   15 CA   CT        0    0    0    1    1   -0.289820    0.000000
   162HA   HC        0    0    0    1    1    0.120240    0.000000
   173HA   HC        0    0    0    1    1    0.120240    0.000000 
   18 C    C         2    1    0    1    1    0.553070    0.000000
   19 O    O         0    0    0    1    1   -0.613600    0.000000
    2    1
    3    1
    4    3
    5    3
    6    3
    7    6
    8    6
    9    6
   10    9
   11    9
   12    9
   13   12
   14   12
   15   12
   16   15
   17   15
   18   15
   19   18
"""

prepare = """start {name}
memory total {memory} mb

prepare
  system {name}_mm
  source ./{complex}
  new_top new_seq
  new_rst
  {center}
  {orient}
  solvate box {size} {size} {size}
  {counter}
  update lists
  ignore
  write {name}_mm.rst
  write {name}_mm_initial.pdb
end

task prepare

md
  system {name}_mm
  noshake solute
  cutoff 1.2
  pme grid 32 alpha 1e-6 order 6 fft 2
  sd 1000 min 0.0001
  record rest 10
end

task md optimize

prepare
  read {name}_mm.qrs
  write {name}_mm_optimized.pdb
end

task prepare
  
prepare
  system {name}2_qmmm
  source ./{name}_mm_optimized.pdb
  new_top new_seq
  new_rst
  {modify}
  update lists
  ignore
  write {name}2_qmmm.rst
  write {name}2_qmmm_initial.pdb
end

task prepare

md
  system {name}2_qmmm
  noshake solute
  cutoff 1.2 qmmm 1.2
  #pme grid 32 alpha 1e-6 order 6 fft 2 
  record rest 1 scoor 1
end


basis "ao basis" spherical
 * library {aobasis}
end

basis "cd basis" spherical bse
 * library {cdbasis}
end

charge {charge}

dft
 noio
 {disp}
 maxiter {nscf}
 xc {xcfun}
 grid {grid} nodisk
 convergence lshift {lshift} energy 1d-7
 noprint "final vectors analysis"
end

qmmm
  region qmlink mm
  maxiter 10 500
  ncycles {nopt}
  bqzone {bqzone}
  density espfit
  xyz qmregion
end

set driver:linopt 0

task qmmm dft optimize

prepare
  read {name}2_qmmm.qrs
  write {name}2_qmmm_optimized.pdb
end

task prepare

analysis
  system {name}2_qmmm
  reference {name}2_qmmm.rst
  file {name}2_qmmm.trj
  copy {name}2_qmmm.xyz
end
task analysis
"""

amber2nwchem = """#!/bin/bash

sed -i \
  -e '/ ALA /{s/ HB1/2HB /;s/ HB2/3HB /;s/ HB3/4HB /}' \
  -e '/ ARG /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /;s/ HD2/2HD /;s/ HD3/3HD /;s/HH11/2HH1/;s/HH12/3HH1/;s/HH21/2HH2/;s/HH22/3HH2/}' \
  -e '/ ASH /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ ASN /{s/ HB2/2HB /;s/ HB3/3HB /;s/HD21/2HD2/;s/HD22/3HD2/}' \
  -e '/ ASP /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ CYM /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ CYS /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ CYX /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ GLH /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /}' \
  -e '/ GLN /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /;s/HE21/2HE2/;s/HE22/3HE2/}' \
  -e '/ GLU /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /}' \
  -e '/ GLY /{s/ HA2/2HA /;s/ HA3/3HA /}' \
  -e '/ HID /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ HIE /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ HIP /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ ILE /{s/HG21/2HG2/;s/HG22/3HG2/;s/HG23/4HG2/;s/HG12/2HG1/;s/HG13/3HG1/;s/HD11/2HD /;s/HD12/3HD /;s/HD13/4HD /}' \
  -e '/ LEU /{s/ HB2/2HB /;s/ HB3/3HB /;s/HD11/2HD1/;s/HD12/3HD1/;s/HD13/4HD1/;s/HD21/2HD2/;s/HD22/3HD2/;s/HD23/4HD2/}' \
  -e '/ LYN /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG/;s/ HG3/3HG /;s/ HD2/2HD /;s/ HD3/3HD /;s/ HE2/2HE /;s/ HE3/3HE /;s/ HZ2/2HZ /;s/ HZ3/3HZ /}' \
  -e '/ LYS /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /;s/ HD2/2HD /;s/ HD3/3HD /;s/ HE2/2HE /;s/ HE3/3HE /;s/ HZ1/2HZ /;s/ HZ2/3HZ /;s/ HZ3/4HZ /}' \
  -e '/ MET /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /;s/ HE1/2HE /;s/ HE2/3HE /;s/ HE3/4HE /}' \
  -e '/ PHE /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ PRO /{s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /;s/ HD2/2HD /;s/ HD3/3HD /}' \
  -e '/ SER /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ THR /{s/HG21/2HG2/;s/HG22/3HG2/;s/HG23/4HG2/}' \
  -e '/ TRP /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ TYR /{s/ HB2/2HB /;s/ HB3/3HB /}' \
  -e '/ VAL /{s/HG21/2HG2/;s/HG22/3HG2/;s/HG23/4HG2/;s/HG11/2HG1/;s/HG12/3HG1/;s/HG13/4HG1/}' \
  -e '/ ACE /{s/HH31/2HH3/;s/HH32/3HH3/;s/HH33/4HH3/}' \
  -e '/ NME /{s/HH31/2HH3/;s/HH32/3HH3/;s/HH33/4HH3/}' \
  -e '/ EAC /{s/ HA2/2HA /;s/ HA3/3HA /;s/ HB2/2HB /;s/ HB3/3HB /;s/ HG2/2HG /;s/ HG3/3HG /;s/ HD2/2HD /;s/ HD3/3HD /;s/ HE2/2HE /;s/ HE3/3HE /}' \
  -e '{s/ H1/2H /;s/ H2/3H /;s/ H3/4H /}' \
  $1"""
