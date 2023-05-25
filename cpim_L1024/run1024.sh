#!/bin/bash
#SBATCH -J CPIM1024
#SBATCH -p slims
#SBATCH -n 1
#SBATCH --output=output_%j.out
#SBATCH --error=errors_%j.err
#SBATCH --mail-user=ajlhomme@uc.cl
#SBATCH --mail-type=ALL

make
./CPIM "CPIM_L=1024_2.txt"