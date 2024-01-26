#!/bin/bash
#SBATCH --export=ALL
#SBATCH -J SNAKEMASTER
#SBATCH -o __LOG_DIRECTORY__/workflow_output_%j.log
#SBATCH --ntasks 1
#SBATCH --time 240:00:00
#SBATCH --mem=8G

#-----------------------------------------------------------------------------------------------------------------------
# STORMqc SLURM Submit Template
#
# Config Creation Data:
#     CREATE_DATE
#
# Template Created by:
#     Jacob Morrison
#
# Template Creation Date:
#     January 2024
#
# Notes:
#     Jan 2024 -
#         Initial creation
#-----------------------------------------------------------------------------------------------------------------------

cd ${SLURM_SUBMIT_DIR}

CONFIG_FILE=__CONFIG_YAML__
LOG_DIR=__LOG_DIRECTORY__

# Add snakemake to PATH here
if [[ `which snakemake 2>&1 /dev/null` ]]; then
    snakemake_module="bbc2/snakemake/snakemake-7.32.3"

    module load $snakemake_module
fi

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --configfile ${CONFIG_FILE} --dry-run         > ${LOG_DIR}/workflow_${TIME}.txt
snakemake --configfile ${CONFIG_FILE} --dag | dot -Tpng > ${LOG_DIR}/workflow_${TIME}.png

# Default to using conda, if using environment modules, then replace --use-conda with --use-envmodules
# Note, this requires downloading mamba (conda install -n base -c conda-forge mamba)
snakemake \
    --printshellcmds \
    --latency-wait 20 \
    --keep-going \
    --use-conda \
    --jobs 20 \
    --configfile ${CONFIG_FILE} \
    --cluster "mkdir -p ${LOG_DIR}/{rule}; sbatch \
        --export=ALL \
        --ntasks {threads} \
        --mem={resources.mem_gb}G \
        -t {resources.time} \
        -o ${LOG_DIR}/{rule}/{rule}-%j.log"
