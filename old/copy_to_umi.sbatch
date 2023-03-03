#!/usr/bin/env bash
#SBATCH -J copy_umi
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

rsync -avz --include="*_barcode.fastq" --exclude="*unmappedReads*"   --include="*_R1.fastq" --include="*_R2.fastq" --include="*allCRs*" --exclude="*"  BarcodeCounts_technical_replicates/  BarcodeCounts_umi_analysis/