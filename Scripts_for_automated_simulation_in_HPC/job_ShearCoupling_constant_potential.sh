#!/bin/bash
#SBATCH --account=MY_GROUP
#SBATCH -N 3
#SBATCH -p std
##SBATCH -C IVY,IB,INTEL
##SBATCH -C SKYLAKE,OPA,INTEL
##SBATCH -C BROADWELL
##SBATCH -C CASCADELAKE,OPA,INTEL
##SBATCH --requeue
#SBATCH -C NOPREEMPT
#SBATCH --exclude=cna[30-38]
#SBATCH -J constant_potential
#SBATCH -n 96
#SBATCH --ntasks-per-node=32
#SBATCH -t 00-03:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=gbizana3@gatech.edu
#SBATCH --output min.output
GB="Sigma_89_3130"
T=50
tau=0.6
TN="$T"
latP=4.0966
Nrun=1000000
mkdir $SLURM_SUBMIT_DIR/Syst/ $SLURM_SUBMIT_DIR/Syst/$TN/  $SLURM_SUBMIT_DIR/Syst/$TN/$tau $SLURM_SUBMIT_DIR/Syst/$TN/$tau/constant_potential   $SLURM_SUBMIT_DIR/Syst/$TN/$tau/constant_potential/${GB} 
#SBATCH --output $SLURM_SUBMIT_DIR/$GB/$TN/CS/$tau
cd $SLURM_SUBMIT_DIR/Syst/$TN/$tau/constant_potential/${GB}
#module purge
#module load intel/2018.5 intelmpi/2018.5.274 lammps-2018/lammps_2018_INTELMPI
module purge
module load intel/2017.4.196 intelmpi/2017.4.196 #lammps-2018/lammps_2018_INTELMPI
#module load intel/2017.4.196 intelmpi/2017.4.196 lammps/11Aug17/intel_cpu_intelmpi
unset I_MPI_PMI_LIBRARY
srun --mpi=pmi2 ~/lammps-27Oct2021/src/lmp_mpi -v gbn $GB -v T_eq $T  -v tau $tau -v Lp ${latP} -v NR $Nrun <$SLURM_SUBMIT_DIR/ShearCoupling_constant_potential.input
#cp $SLURM_SUBMIT_DIR/${SLURM_JOB_ID}.out $SLURM_SUBMIT_DIR/$TN/$tau/constant_potential/${GB}/
