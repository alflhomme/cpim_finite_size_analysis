#!/bin/bash
#SBATCH -J CPIM
#SBATCH -p slims
#SBATCH -n 1
#SBATCH --output=output_%j.out
#SBATCH --error=errors_%j.err
#SBATCH --mail-user=ajlhomme@uc.cl
#SBATCH --mail-type=ALL

cd cpim_L1024
make
./CPIM "CPIM_L=1024.txt"
cd ..
#
cd cpim_L512
make
./CPIM "CPIM_L=512.txt"
cd ..
#
cd cpim_L256
make
./CPIM "CPIM_L=256.txt"
cd ..
#
cd cpim_L128
make
./CPIM "CPIM_L=128.txt"
cd ..
#
cd cpim_L64
make
./CPIM "CPIM_L=64.txt"
cd ..
#
cd cpim_L32
make
./CPIM "CPIM_L=32.txt"
cd ..