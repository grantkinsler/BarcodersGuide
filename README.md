Barcoder's Guide README

----------------------------
Processing Sequencing Data.
----------------------------

Sequencing data was processed using BarcodeCounter.py

First, run demultiplex_b

For analysis involving index hopping, sample files were created with all possible combinations of N,S,F,R primers (96x96=9216 total samples). 

SequenceToCount/running_sbatch.py 

Relevant lanes:
	070418_diagonal	- HiSeqX run with indices on the diagonal
	080618_NextSeq 	- NextSeq run with 96 samples (same library as 082818)
	082818_HiSeq	- HiSeqX run with 96 samples (same library as 080618)

Created job array for each lane - submits 96 jobs (for each NS primer combo), each of which has a sample sheet with 96 samples (the F,R primers)


----------------------------
Analyzing Fitness Data
----------------------------

Code corresponding to the creation of each figure is contained in the code directory. 
