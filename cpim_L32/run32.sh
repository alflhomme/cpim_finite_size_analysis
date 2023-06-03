#!/bin/bash
#SBATCH -J CPIM32
#SBATCH -n 1
#SBATCH --output=output_%j.out
#SBATCH --error=errors_%j.err
#SBATCH --mail-user=ajlhomme@uc.cl
#SBATCH --mail-type=ALL

make
./CPIM "CPIM_L=32.txt"