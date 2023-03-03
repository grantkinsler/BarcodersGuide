#!/usr/bin/env bash
#SBATCH -J iHop_082818_Hiseq_map_barcodes.sbatch
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
python3 /home/groups/dpetrov/SOFTWARE/BarcodeCounter2-master_07272021update/barcodeCounter.py -fastqDir ../AllRawData/ -outputDir ../BarcodeCounts_index_hopping/082818_Hiseq/ -templateSeq ../Template_doubleBC_Read1First.txt -sample iHop_082818_Hiseq_All_Samples.txt -multiBCFasta ../Primers_noP.fasta -pairedEnd -useUMI -numThreads 16 -skipSplitFastq -barcodeList ../500pool_noconstant_bothBCs_withSpikeIns.fasta -useBowtie2 -bowtie2Path /home/groups/dpetrov/SOFTWARE/bowtie2-2.2.6_new/ -readLength 108 -remapBarcodes
