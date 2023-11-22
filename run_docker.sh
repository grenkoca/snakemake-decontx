# Get repo to pass to Docker
SNK_REPO="$(pwd)"

# Run command inside Docker image -- snakemake only uses singularity natively
mkdir -p results # make results dir to mount 
docker run --name islet_decontamination \
     --mount type=bind,source="${SNK_REPO}",target=/home/container_user/predict \
     --mount type=bind,source="${PWD}",target=/home/container_user/analysis \
     --rm -t \
     letaylor/sc_decontx:latest \
     /bin/bash -c "cd /home/container_user/analysis && \
          source activate base && \
          chmod +777 * && \
	  snakemake -s Snakefile -j 999 -r --latency-wait=180 --local-cores 8 --cluster-config ./lib/threeprime_cluster.yaml --rerun-incomplete --max-jobs-per-second=50 --printshellcmds --configfile workflow/src/threeprime.yaml
          "
