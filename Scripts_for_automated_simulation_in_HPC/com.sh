#!/bin/bash
#SBATCH --account=MY_GROUP
#SBATCH -N 1
#SBATCH -p std
##SBATCH -C IVY,IB,INTEL
##SBATCH -C SKYLAKE,OPA,INTEL
##SBATCH -C BROADWELL,OPA,INTEL
##SBATCH -C CASCADELAKE,OPA,INTEL
##SBATCH --requeue
#SBATCH -C NOPREEMPT
#SBATCH --exclude=cna[30-38]
#SBATCH -J minimize_shear
#SBATCH -n 32
#SBATCH --ntasks-per-node=32
#SBATCH -t 00-10:30:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=gbizana3@gatech.edu
#SBATCH --output zip_750.output
zip  -9 -r 250_arch.zip  /home/uda69/lcz56/Gen_GB_min_Final_109/250
