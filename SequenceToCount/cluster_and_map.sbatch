#!/usr/bin/env bash
#SBATCH -J FitnessMeasurements
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

module load python/3.6.1
module load biology
module load ncbi-blast+/2.7.1
python3 /home/groups/dpetrov/SOFTWARE/BarcodeCounter2-master_06242019update/barcodeCounter.py -fastqDir AllRawData/ -outputDir BarcodeCounts/ -templateSeq Template_double_longbc.txt -sample All_samples.txt -multiBCFasta Primers_noP.fasta -pairedEnd -useUMI -numThreads 16 -skipSplitFastq -DNAclustPath /home/groups/dpetrov/SOFTWARE/dnaclust_linux_release3/ -useBowtie2 -bowtie2Path /home/groups/dpetrov/SOFTWARE/bowtie2-2.2.6_new/ -readLength 108
