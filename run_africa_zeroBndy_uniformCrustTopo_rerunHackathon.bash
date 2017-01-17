#!/bin/bash
#SBATCH -J run_africa_rerunHackathon
#SBATCH -o africa_zeroBndy_uniformCrustTopo_rerunHackathon.o%j
#SBATCH -n 5 
#SBATCH -p normal 
#SBATCH -t 00:30:00
ibrun $HOME/packages/aspect/build/aspect africa_zeroBndy_uniformCrustTopo_rerunHackathon.prm
