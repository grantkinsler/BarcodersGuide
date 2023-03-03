README.txt


Processing counts for Barcoder's Guide. 
_____________________________________________________
SUMMARY
(1) Demultiplex the samples
    sbatch techRep_demultiplex.sbatch
(2) Build the bowtie index for the barcodes. 
    /home/groups/dpetrov/SOFTWARE/bowtie2-2.2.6_new/bowtie2-build ../500pool_noconstant_bothBCs_withSpikeIns_withAncestor.fasta ../500pool_noconstant_bothBCs_withSpikeIns_withAncestor.fasta
(3) sbatch techRep_map_barcodes.sbatch
(4) sbatch techRep_count_barcodes.sbatch

For a comparison of using different UMI metrics for each sample, also run:
(5) sbatch umi_comparision.sbatch
_____________________________________________________

(1) Run techRep_demuliplex.sbatch 

	This demultiplexes the samples based on the Nextera indicies (N/S primers) (already demultiplexed by sequencing company into separate files) as well as the inline indicies (F/R primers).

    This step also extracts the relevant regions of the read (index sequences, UMIs, and barcodes) and writes them to corresponding files.

	During demuliplexing, I received errors from Groups 96, 138, 148. Reports of non-gzipped files and/or fastq records not starting with @.
	Corresponds with following samples/sequencing files.

	EE1-DE1-PCRa        090118_Nextseq_N716S518
	U3-DE1-PCRa        090118_Nextseq_N726S521
	W3-DE1-PCRa        120118_NextseqHTP_N728S521

    This doesn't matter for this study, so I'm ignorning it for now - but if you need these files, you may want to double check they are the correct file and/or not corrupted.

(2) Build the bowtie index for the barcodes.

    This builds the bowtie index. Note that if you are only using one "group" of samples, you can remove the -skipBowtieBuild option in the code. 

    When processing the data with many groups of samples, it's best to just manually run this, since otherwise there will be lots of attempts to delete/rewrite the bowtie index which will cause jobs to fail

(3) Map barcodes

    Maps the barcode reads to the bowtie index, finding only barcodes that are supposed to be in the list.

    


(4) Count barcodes

    Aggregates all of the counts from all of the separate groups into a single count file.

