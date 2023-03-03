#!/usr/bin/env bash
#SBATCH -J techRep_demultiplex
#SBATCH -p hns,dpetrov
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -t 7-00:00
#SBATCH --array=1,96,138,148,267,268,269,270,271
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

sample_sheets=$(sed -n "$SLURM_ARRAY_TASK_ID"p techRep_demultiplex.inp)

module load python/3.6.1
module load biology
module load ncbi-blast+/2.7.1

python3 

python3 /home/groups/dpetrov/SOFTWARE/BarcodeCounter2-master_11162021update/barcodeCounter.py -fastqDir ../AllRawData/ -outputDir ../BarcodeCounts_fixed/ -templateSeq ../Template_doubleBC_Read1First.txt -sample $sample_sheets -multiBCFasta ../Primers_noP.fasta -pairedEnd -useUMI -numThreads 32 -demultiplexOnly -useBowtie2 -bowtie2Path /home/groups/dpetrov/SOFTWARE/bowtie2-2.2.6_new/ -readLength 108