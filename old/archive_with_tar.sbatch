#!/usr/bin/env bash
#SBATCH -J archive
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 1-00:00
#SBATCH --array=1-1
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

tar -czvf BarcodersGuide_10Nov2021.tar.gz BarcodersGuide