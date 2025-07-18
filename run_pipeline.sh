#!/bin/bash
##############################################3

export TMPDIR="/scratch/tmp/"

CMD="snakemake \
-s Snakefile \
-j 999 \
-r \
--latency-wait=180 \
--local-cores 8 \
--cluster-config ./lib/threeprime_cluster.yaml \
--rerun-incomplete \
--max-jobs-per-second=50 \
--printshellcmds"

while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
  -u | --unlock)
    CMD="$CMD --unlock"
    shift
    ;;
  -d | --dryrun)
    CMD="$CMD --dryrun"
    shift
    ;;
  --cluster)
    CMD="$CMD --cluster $2"
    shift
    shift
    ;;
  --configfile)
    CMD="$CMD --configfile $2"
    shift
    shift
    ;;
  esac
done

echo "EXECUTING COMMAND"
echo $CMD
$CMD
