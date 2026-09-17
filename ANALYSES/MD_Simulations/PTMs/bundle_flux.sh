#!/bin/bash
flux resource list
while IFS= read -r line; do
  dir=$(dirname $line)
  base=$(basename $line)
  flux batch -N 1 -t 23.8h -l --gpus-per-slot=8 --cores-per-slot=56 -x --cwd=$dir --output=$base.flux.out --error=$base.flux.err $line
done < fluxjobs.txt
flux queue drain

maxestruns=0
while IFS= read -r line; do
  dir=$(dirname $line)
  estruns=$(cat ${dir}/estruns.txt)
  if [ $estruns -gt $maxestruns ]; then
    maxestruns=$estruns
  fi
done < md_fluxjobs.txt
echo $maxestruns > maxestruns.txt
echo "Number of additional runs needed: $maxestruns"

# Obtain the jobid of the first run
jobid=$(cat md.jobid)

for i in $(seq 1 $maxestruns); do
  sleep 1s
  echo "Batching run $i"
  jobid=$(sbatch --dependency=afterok:$jobid md_bundle.sbatch | sed 's/Submitted batch job //')
  echo "Submitted job $jobid"
done

echo $jobid > $md.jobid
