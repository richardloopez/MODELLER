#!/bin/bash
#SBATCH --job-name=P_MODELLER_JOB
#SBATCH --output=salida_slurm_%j.out
#SBATCH --error=error_slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48  ###################

module load Python/3.11.3-GCCcore-12.3.0

export MODELLER_CORES=$SLURM_CPUS_PER_TASK

srun --cpus-per-task=$SLURM_CPUS_PER_TASK python3 controller.py > salida.out

