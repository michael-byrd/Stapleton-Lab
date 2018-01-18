#!/bin/bash 

#SBATCH -J ManchingMultiple 

#SBATCH -o ManchingMultiple%j.txt

#SBATCH -N 1

#SBATCH -n 68

#SBATCH -p normal

#SBATCH -t 03:00:00 

#SBATCH -A TG-DMS170019

#SBATCH --mail-user=trb9259@uncw.edu

#SBATCH --mail-type=end 

sadfasdf

cd $HOME/MultipleStress # Change to home directory first, just to make sure SLURM reads from the proper directory

module load Rstats # Load the R module along with some popular packages so it will run R files

Rscript  famrandfinal.R # You may need to change this depending on what exactly With.R is, but if it's a standard R script that can be run from command line, this would suffice

### 4.4min for core=276, n.boot=100
