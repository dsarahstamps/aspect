#!/bin/bash 
#PBS -N g5_vtu
#PBS -o g5_vtu 
#PBS -l nodes=1:ppn=10
#PBS -q stamps_geo_lab
#PBS -l walltime=48:00:00 
#PBS -W group_list=stamps_geo_lab
#PBS -A mantLith

module purge
module load gcc openmpi/1.10.0 glm lua python boost hdf5 phdf5 atlas trilinos p4est dealii
cd $WORK/aspect
mpirun -np $PBS_NP ./build/aspect g5_traction_vector_rerun_v4.prm | tee jobout_g5_vtu.$PBS_JOBID
