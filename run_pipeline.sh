#!/bin/bash
##############################################3


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/grenkocm/conda/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/grenkocm/conda/etc/profile.d/conda.sh" ]; then
        . "/home/grenkocm/conda/etc/profile.d/conda.sh"
    else
        export PATH="/home/grenkocm/conda/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/home/grenkocm/conda/etc/profile.d/mamba.sh" ]; then
    . "/home/grenkocm/conda/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<

export TMPDIR="/scratch/tmp/"

CMD=\
"snakemake \
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
    -u|--unlock)
      CMD="$CMD --unlock"
      shift
      ;; 
    -d|--dryrun)
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

