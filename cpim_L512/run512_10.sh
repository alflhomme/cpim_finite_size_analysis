#!/bin/bash
#SBATCH -J CPIM512_10
#SBATCH -p slims
#SBATCH -n 1
#SBATCH --output=output_%j.out
#SBATCH --error=errors_%j.err
#SBATCH --mail-user=ajlhomme@uc.cl
#SBATCH --mail-type=ALL

./CPIM "CPIM_L=512_11.txt"