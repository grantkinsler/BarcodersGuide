#!/usr/bin/env bash
#SBATCH -J umi_comparison
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-271
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

sample_sheets=$(sed -n "$SLURM_ARRAY_TASK_ID"p techRep_demultiplex.inp)

module load python/3.6.1
module load biology
module load ncbi-blast+/2.7.1
python3 umi_comparison.py -fastqDir ../AllRawData/ -outputDir ../BarcodeCounts_technical_replicates/ -templateSeq ../Template_doubleBC_Read1First.txt -sample $sample_sheets -multiBCFasta ../Primers_noP.fasta -pairedEnd -useUMI -numThreads 16 -skipSplitFastq -barcodeList ../500pool_noconstant_bothBCs_withSpikeIns_withAncestor.fasta -useBowtie2 -bowtie2Path /home/groups/dpetrov/SOFTWARE/bowtie2-2.2.6_new/ -readLength 108 -remapBarcodes
