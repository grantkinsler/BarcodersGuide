#!/usr/bin/env bash
#SBATCH -J copy_data
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-1
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

ml system
ml rclone

rclone copy teamdrive:yeast_fitness_remeasurements/1BigBatch_AllRawData AllRawData

rclone copy teamdrive:yeast_fitness_remeasurements/FRM_Data_From_Yuping AllRawData
