#!/bin/bash
#SBATCH --job-name=Modeller_Job
#SBATCH --output=salida_slurm_%j.out
#SBATCH --error=error_slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16  ###################

module load python/3.8.12

export MODELLER_CORES=$SLURM_CPUS_PER_TASK

srun --cpus-per-task=$SLURM_CPUS_PER_TASK python3 controller.py > salida.out

