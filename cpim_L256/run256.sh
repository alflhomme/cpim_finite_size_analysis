#!/bin/bash
#SBATCH -J CPIM256
#SBATCH -p slims
#SBATCH -n 1
#SBATCH --output=output_%j.out
#SBATCH --error=errors_%j.err
#SBATCH --mail-user=ajlhomme@uc.cl
#SBATCH --mail-type=ALL

make
./CPIM "CPIM_L=256.txt"