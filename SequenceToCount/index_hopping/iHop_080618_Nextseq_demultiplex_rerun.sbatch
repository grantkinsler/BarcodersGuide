#!/usr/bin/env bash
#SBATCH -J iHop_080618_Nextseq_demultiplex.sbatch
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=31
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

sample_sheets=$(sed -n "$SLURM_ARRAY_TASK_ID"p iHop_080618_Nextseq_demultiplex.inp)

module load python/3.6.1
module load biology
module load ncbi-blast+/2.7.1
python3 /home/groups/dpetrov/SOFTWARE/BarcodeCounter2-master_06242019update/barcodeCounter.py -fastqDir ../AllRawData/ -outputDir ../BarcodeCounts_index_hopping/080618_Nextseq/ -templateSeq ../Template_doubleBC_Read1First.txt -sample $sample_sheets -multiBCFasta ../Primers_noP.fasta -pairedEnd -useUMI -numThreads 16 -demultiplexOnly -useBowtie2 -bowtie2Path /home/groups/dpetrov/SOFTWARE/bowtie2-2.2.6_new/ -readLength 108
