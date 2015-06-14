#!/bin/bash
#SBATCH -J g5_run_africa_litho1.0_crust1.0_smooth_harmAvg_tractionVector
#SBATCH -o g5_TV_africa_litho1.0_crust1.0_smooth_harmAvg_tractionVector.o%j
#SBATCH -n 64 
#SBATCH -p normal 
#SBATCH -t 20:00:00
ibrun $HOME/packages/aspect/build/aspect africa_litho1.0_crust1.0_smooth_harmAvg.prm
