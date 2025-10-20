#!/bin/bash
#SBATCH --job-name=Modeller_Job      # Nombre del trabajo
#SBATCH --output=salida_slurm_%j.out # Nombre del archivo de salida
#SBATCH --error=error_slurm_%j.err   # Nombre del archivo de errores
#SBATCH -n 288
##SBATCH --ntasks=1                   # Solo 1 tarea principal
##SBATCH --cpus-per-task=4           # Solicitar procesadores (CPUs)

# 1. Cargar el entorno Python
module load python/3.8.12

# 2. Exportar la variable de Slurm (número real de CPUs asignados)
# La variable $SLURM_CPUS_PER_TASK contendrá 20 en este ejemplo.
export MODELLER_CORES=$SLURM_CPUS_PER_TASK

# 3. Ejecutar el script principal
srun --cpus-per-task=288 python3 controller.py > salida.out





