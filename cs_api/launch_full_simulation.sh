#!/bin/bash
#-------------------------------------------------------------------------------
#
# Batch options for SLURM (Simple Linux Utility for Resource Management)
# =======================
#SBATCH --nodes=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=bm
#SBATCH --wckey=P120R:CODE_SATURNE
#SBATCH --output=farm.out.log
#SBATCH --error=farm.err.log
#SBATCH --job-name=farm

mesh_jobid=$(sbatch --parsable launch_mesh.sh)
echo $mesh_jobid

#Run  wind farm simulation
sbatch --dependency=afterok:$mesh_jobid launch_code_saturne.sh
