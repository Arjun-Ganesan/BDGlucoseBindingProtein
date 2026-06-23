#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --job-name 6471_1_ligmpnn62_BGC   ##name that will show up in the queue
#SBATCH --output /home/qzhong/AlphaFold3/Outputs/05122025/6471_1_ligmpnn62_BGC-%j.out
#SBATCH --error /home/qzhong/AlphaFold3/Outputs/05122025/6471_1_ligmpnn62_BGC-%j.err
#SBATCH --nodes 1  ##number of nodes to use
#SBATCH --ntasks 1  ##number of tasks (analyses) to run
#SBATCH --time 1-00:00:00  ##time for analysis (day-hour:min:sec)
#SBATCH --cpus-per-task 48  ##the number of threads the code will use
#SBATCH --partition contrib-gpuq  ##the partition to run in [options: normal, gpu, debug]
#SBATCH --qos gpu
#SBATCH --gres gpu:A100.80gb:1  ##GPU(s) allocated
#SBATCH --mem 32GB


## Load modules, insert code, and run your programs here.
#load the singularity module:

module load singularity

#set environments for the databases and model parameters:

DB_DIR=/datasets/alphafold3/databases
MODEL_PARAMETERS_DIR=/datasets/alphafold3/model_parameters
AF3_SCRIPTS=/containers/dgx/Containers/alphafold3/v3.0.1

#set the environment commands to run the image

AF3_IMAGE=/containers/dgx/Containers/alphafold3/v3.0.1/singularity/alphafold3.sif

INPUT_DIR=/home/qzhong/AlphaFold3/Inputs/05122025
OUTPUT_DIR=/home/qzhong/AlphaFold3/Outputs/05122025

SINGULARITY_RUN="singularity exec --nv \
    --bind $MODEL_PARAMETERS_DIR:/root/models \
    --bind $DB_DIR:/root/public_databases \
    --bind $AF3_SCRIPTS:/root/scripts \
    --bind $INPUT_DIR:/root/inputdir \
    --bind $OUTPUT_DIR:/root/outputdir"

#run alphafold3 scripts:

${SINGULARITY_RUN} ${AF3_IMAGE} python3 /root/scripts/run_alphafold.py \
  --json_path=/root/inputdir/6471_1_ligmpnn62_BGC.json \
  --model_dir=/root/models \
  --db_dir=/root/public_databases \
  --output_dir=/root/outputdir/6471_1_ligmpnn62_BGC