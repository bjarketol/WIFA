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
#SBATCH --output=mesh.out.log
#SBATCH --error=mesh.err.log
#SBATCH --job-name=mesh

/software/rd/salome/logiciels/salome/V9_10_0/salome -t python3 generate_salome_mesh.py
