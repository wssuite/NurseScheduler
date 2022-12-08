#!/bin/bash
#SBATCH --account=def-foo   # some account
#SBATCH --time=0-4:01        # specify time limit (D-HH:MM)
#SBATCH --cpus-per-task=8     # specify number threads
#SBATCH --mem=32G              # specify total memory
#SBATCH --nodes=1             # number of nodes
#SBATCH --array=0-19
#SBATCH --output=slurm/%A_%a.out

module load StdEnv/2020       # for versions > 9.0.2
module load gurobi/9.5.0
module load openblas

./run-benchmark.py -b benchmark/w4.yml -t ${SLURM_ARRAY_TASK_ID}-${SLURM_ARRAY_TASK_ID} --prefix ${SLURM_ARRAY_TASK_ID}_