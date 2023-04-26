Barcoder's Guide README

This is the Github Repository for:

Kinsler et al. 2023. Extreme sensitivity of fitness to environmental conditions; lessons from #1BigBatch. 

The "Processing Sequencing Data" section contains details about how sequencing data was processed to produce barcode counts.

The "Analysizing Fitness Data" section contains details on how to find code for analysis and figures in the paper.

----------------------------
Processing Sequencing Data.
----------------------------

All processing steps are provided in the "SequenceToCount" directory.

Sequencing data was processed using BarcodeCounter.py

First, run demultiplex_bc_files.sbatch. This step only uses the initial steps of BarcodeCounter.py which demuliplexes samples and extracts barcode/UMI sequences from each read.
Second, run cluster_and_map.sbatch. This step maps barcodes sequences to the provided list of barcodes in parallel.
Finally, run count_barcodes.sbatch. This step compiles all counts of mapped barcodes across all samples into a single table.

These steps were initially run all samples.

For analysis involving index-misassignment and template switching, additional sample files were created with all possible combinations of N,S,F,R primers (96x96=9216 total samples) rather than the true samples. For this analysis, running_sbatch.py was used.

Relevant lanes:
	070418_diagonal	- HiSeqX run with indices on the diagonal (just 8 samples)
	080618_NextSeq 	- NextSeq run with 96 samples (same exact library as 082818_HiSeq)
	082818_HiSeq	- HiSeqX run with 96 samples (same exact library as 080618_NextSeq)

----------------------------
Analyzing Fitness Data
----------------------------

Code corresponding to the analysis and the creation of each figure is contained in the code directory.

code/figure1+4 contains R scripts used to create figures 1 and 4.
code/figure2 contains R scripts used to create figure 3.
code/figure3 contains python notebook used to analyze index mis-assignment and generate figure 3


