#!/bin/bash
flux resource list
while IFS= read -r line; do
  dir=$(dirname $line)
  base=$(basename $line)
  flux batch -N 1 -t 23.8h -l --gpus-per-slot=8 --cores-per-slot=56 -x --cwd=$dir --output=$base.flux.out --error=$base.flux.err $line
done < md_fluxjobs.txt
flux queue drain
