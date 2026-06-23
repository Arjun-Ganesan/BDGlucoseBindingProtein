#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --job-name 17_53_105_Fixed_ligmpnn_ATPs   ##name that will show up in the queue
#SBATCH --output /home/qzhong/AlphaFold3/Outputs/07032025/17_53_105_Fixed/17_53_105_Fixed-%j.out
#SBATCH --error /home/qzhong/AlphaFold3/Outputs/07032025/17_53_105_Fixed/17_53_105_Fixed-%j.err
#SBATCH --nodes 1  ##number of nodes to use
#SBATCH --ntasks 1  ##number of tasks (analyses) to run
#SBATCH --time 1-00:00:00  ##time for analysis (day-hour:min:sec)
#SBATCH --cpus-per-task 64  ##the number of threads the code will use
#SBATCH --partition contrib-gpuq  ##the partition to run in [options: normal, gpu, debug]
#SBATCH --qos gpu
#SBATCH --gres gpu:A100.80gb:1  ##GPU(s) allocated
#SBATCH --mem 320GB

## Load modules, insert code, and run your programs here.
#load the singularity module:

module load singularity

#set environments for the databases and model parameters:

DB_DIR=/datasets/alphafold3/databases
MODEL_PARAMETERS_DIR=/datasets/alphafold3/model_parameters
AF3_SCRIPTS=/containers/dgx/Containers/alphafold3/v3.0.1

#set the environment commands to run the image

AF3_IMAGE=/containers/dgx/Containers/alphafold3/v3.0.1/singularity/alphafold3.sif

GROUP=17_53_105_Fixed
INPUT_DIR=/home/qzhong/AlphaFold3/Inputs/07032025/${GROUP}
OUTPUT_DIR=/home/qzhong/AlphaFold3/Outputs/07032025/${GROUP}

SINGULARITY_RUN="singularity exec --nv \
    --bind $MODEL_PARAMETERS_DIR:/root/models \
    --bind $DB_DIR:/root/public_databases \
    --bind $AF3_SCRIPTS:/root/scripts \
    --bind $INPUT_DIR:/root/inputdir \
    --bind $OUTPUT_DIR:/root/outputdir"

#run alphafold3 scripts:

${SINGULARITY_RUN} ${AF3_IMAGE} bash -c '
for json in /root/inputdir/*.json; do
  seq=$(basename "$json" .json)
  echo ">>> Running $seq"
  python3 /root/scripts/run_alphafold.py \
      --json_path="$json" \
      --model_dir=/root/models \
      --db_dir=/root/public_databases \
      --output_dir=/root/outputdir/"$seq"
done
'