#!/bin/bash

#SBATCH -J sim_vQTL           # Job name
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p long          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 96:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=mcb4548@uncw.edu
#SBATCH --mail-type=all    # Send email at begin and end of job


# Other commands must follow all #SBATCH directives...

mkdir output
Rscript --vanilla --verbose ./SimulatedvQTL.R > ./output.Rout
