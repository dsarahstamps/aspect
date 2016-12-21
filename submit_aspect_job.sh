#!/bin/bash 
#PBS -N g3_point_t1
#PBS -o g3_point_t1
#PBS -l nodes=1:ppn=20
#PBS -q stamps_geo_lab
#PBS -l walltime=4:00:00 
#PBS -W group_list=stamps_geo_lab
#PBS -A mantLith

module purge
module load gcc openmpi/1.10.0 glm lua python boost hdf5 phdf5 atlas trilinos p4est dealii
cd $WORK/aspect
mpirun -np $PBS_NP ./build/aspect g3_traction_vector_t1.prm | tee jobout_g3_point_t1.$PBS_JOBID
